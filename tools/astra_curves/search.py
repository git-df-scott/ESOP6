#!/usr/bin/env python3
"""Bounded curve-coefficient searches. Floating endpoints are NEVER certificates."""
import argparse, json, math, time
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.optimize import least_squares

OUT=Path('results/astra_curves_2026_09_05')
def mons(n,d):
    if n==1: return [(d,)]
    return [(i,)+m for i in range(d,-1,-1) for m in mons(n-1,d-i)]
def rank(a):
    s=np.linalg.svd(a,compute_uv=False)
    return {str(t):int(np.count_nonzero(s>t*max(s[0],1))) for t in (1e-6,1e-8,1e-10)}
def power(p,k):
    r=np.array([1.],dtype=p.dtype)
    for _ in range(k): r=np.convolve(r,p)
    return r

class RationalQuartic:
    def __init__(self,double=False):
        self.double=double
        # p_i(s,t)=sum coefficient[k] s^(4-k)t^k.
        # Common scalar p1[0]=1. Single symmetry uses t scaling p1[1]=1.
        if not double:
            self.n=14; self.base=np.zeros((6,5));self.base[0,0]=self.base[1,0]=1
            self.base[0,1]=1;self.base[1,1]=-1
            self.A=np.zeros((6,5,self.n)); j=0
            for k in (2,3,4):
                self.A[0,k,j]=1;self.A[1,k,j]=(-1)**k;j+=1
            for k in range(5):
                self.A[2,k,j]=1;self.A[3,k,j]=(-1)**k;j+=1
            for i in (4,5):
                for k in (0,2,4):self.A[i,k,j]=1;j+=1
            self.rows=np.arange(0,25,2)
        else:
            self.n=8;self.base=np.zeros((6,5));self.A=np.zeros((6,5,self.n))
            self.base[0,0]=self.base[1,0]=self.base[2,4]=self.base[3,4]=1
            for k in range(1,5):
                self.A[0,k,k-1]=1;self.A[1,k,k-1]=(-1)**k
                self.A[2,4-k,k-1]=1;self.A[3,4-k,k-1]=(-1)**k
            for i,j in [(4,4),(5,6)]:
                self.A[i,0,j]=self.A[i,4,j]=1;self.A[i,2,j+1]=1
            self.rows=np.arange(0,13,2)
        self.m=len(self.rows);self.sign=np.array([1,1,1,1,1,-1])
        self.weights=np.sqrt(np.array([math.comb(24,int(k)) for k in self.rows]))
    def coefficients(self,x): return self.base+np.einsum('ikj,j->ik',self.A,x)
    def both(self,x):
        ps=self.coefficients(x);v=np.zeros(25,dtype=x.dtype);J=np.zeros((25,self.n),dtype=x.dtype)
        for i,p in enumerate(ps):
            p5=power(p,5);v+=self.sign[i]*np.convolve(p5,p)
            for k in range(5):J[k:k+21]+=6*self.sign[i]*p5[:,None]*self.A[i,k][None,:]
        return v[self.rows]/self.weights,J[self.rows]/self.weights[:,None]
    def seed(self,rng,complex_mode=False):
        x=rng.normal(0,0.6,self.n)
        if self.double: x[3]=rng.uniform(-1,1);x[6]=1.5;x[7]=2
        else: x[-3:]=[1.5,2,1.5]
        if complex_mode:x=x+1j*rng.normal(0,.35,self.n)
        return x
    def diagnose(self,x):
        ps=self.coefficients(x);s=np.linalg.svd(ps,compute_uv=False)
        real=bool(np.max(np.abs(ps.imag))<1e-7) if np.iscomplexobj(ps) else True
        out={'coefficient_rank':int(sum(s>1e-7*max(s[0],1))), 'real_coefficients_heuristic':real}
        if real:
            def definite(p):
                p=p.real
                if abs(p[-1])<1e-7 or abs(p[0])<1e-7:return False
                rr=np.roots(p[::-1]);return bool(all(abs(z.imag)>1e-6 for z in rr))
            out['rhs_definite_heuristic']=definite(ps[5]);out['all_definite_heuristic']=all(definite(p) for p in ps)
            out['no_identically_zero_coordinate_heuristic']=bool(all(np.linalg.norm(p)>1e-6 for p in ps))
        out['coordinate_coefficients']=[[[float(z.real),float(z.imag)] for z in p] for p in ps]
        return out

class EllipticQuartic:
    def __init__(self):
        self.n=self.m=84;self.mm={d:mons(4,d) for d in (1,2,4,5,6)}
        self.ind={d:{a:i for i,a in enumerate(m)} for d,m in self.mm.items()}
        self.free=[i for i,a in enumerate(self.mm[2]) if a not in [(2,0,0,0),(0,2,0,0)]]
        self.bfree=[i for i,a in enumerate(self.mm[4]) if a[0]<2]
        self.pivots=[self.ind[2][(2,0,0,0)],self.ind[2][(0,2,0,0)]]
        self.mul=np.array([[self.ind[6][tuple(x+y for x,y in zip(a,b))] for b in self.mm[4]] for a in self.mm[2]])
        self.weights=np.array([math.factorial(6)/math.prod(math.factorial(k) for k in a) for a in self.mm[6]])
        self.lpweights={d:np.array([math.factorial(d)/math.prod(math.factorial(k) for k in a) for a in self.mm[d]]) for d in (5,6)}
        self.exps={d:np.array(self.mm[d]) for d in (5,6)}
        self.constant=np.zeros(84)
        for i in range(4): a=[0]*4;a[i]=6;self.constant[self.ind[6][tuple(a)]]=1
    def lp(self,l,d):return self.lpweights[d]*np.prod(l[None,:]**self.exps[d],axis=1)
    def unpack(self,x):
        q=np.zeros((2,10),dtype=x.dtype);q[0,self.pivots[0]]=q[1,self.pivots[1]]=1
        q[0,self.free]=x[8:16];q[1,self.free]=x[16:24]
        a=x[24:59];b=np.zeros(35,dtype=x.dtype);b[self.bfree]=x[59:84]
        return x[:4],x[4:8],q,a,b
    def both(self,x):
        l,m,q,a,b=self.unpack(x);v=self.constant.astype(x.dtype)+self.lp(l,6)-self.lp(m,6)
        J=np.zeros((84,84),dtype=x.dtype)
        for off,z,sgn in [(0,l,1),(4,m,-1)]:
            z5=self.lp(z,5)
            for k in range(4):
                ids=[]
                for e in self.mm[5]:ee=list(e);ee[k]+=1;ids.append(self.ind[6][tuple(ee)])
                J[ids,off+k]=sgn*6*z5
        for iq,h,off in [(0,a,8),(1,b,16)]:
            np.add.at(v,self.mul.ravel(),-(q[iq,:,None]*h[None,:]).ravel())
            for j,k in enumerate(self.free):np.add.at(J[:,off+j],self.mul[k],-h)
        for k in range(35):np.add.at(J[:,24+k],self.mul[:,k],-q[0])
        for j,k in enumerate(self.bfree):np.add.at(J[:,59+j],self.mul[:,k],-q[1])
        return v/self.weights,J/self.weights[:,None]
    def seed(self,rng,complex_mode=False):
        x=rng.normal(0,.35,84);x[4:8]=rng.normal(0,.8,4)
        if complex_mode:x=x+1j*rng.normal(0,.3,84)
        # Least-squares initialize the cofactors for this plane and pencil.
        x[24:]=0;v,J=self.both(x);x[24:]=np.linalg.lstsq(J[:,24:],-v,rcond=1e-11)[0]
        return x
    def diagnose(self,x):
        l,m,q,a,b=self.unpack(x);mat=[]
        for qq in q:
            M=np.zeros((4,4),dtype=x.dtype)
            for c,e in zip(qq,self.mm[2]):
                ii=[i for i,k in enumerate(e) for _ in range(k)]
                if ii[0]==ii[1]:M[ii[0],ii[1]]=c
                else:M[ii[0],ii[1]]=M[ii[1],ii[0]]=c/2
            mat.append(M)
        ts=np.array([0,1,-1,2,-2.]);ds=np.array([np.linalg.det(mat[0]+t*mat[1]) for t in ts])
        co=np.linalg.solve(np.vander(ts,5,increasing=True),ds)
        rr=np.roots(co[::-1]);sep=min([abs(rr[i]-rr[j]) for i in range(len(rr)) for j in range(i)]+[0] if len(rr)<4 else [abs(rr[i]-rr[j]) for i in range(4) for j in range(i)])
        regular=bool(abs(co[-1])>1e-7 and sep>1e-4 and np.linalg.norm(co)>1e-6)
        return {'real_coefficients_heuristic':bool(np.max(np.abs(x.imag))<1e-7) if np.iscomplexobj(x) else True,
                'smooth_pencil_heuristic':regular,'pencil_root_separation':float(sep),'det_pencil_coefficients':[[float(z.real),float(z.imag)] for z in co]}

def run(args):
    sys=EllipticQuartic() if args.lane=='elliptic' else RationalQuartic(args.lane=='double')
    rng=np.random.default_rng(args.seed);start=time.time();records=[];seen=[];endpoints=0
    for i in range(args.starts):
        cm=i>=args.real_starts;x0=sys.seed(rng,cm)
        if cm:
            def fun(z):v,j=sys.both(z[:sys.n]+1j*z[sys.n:]);return np.r_[v.real,v.imag]
            def jac(z):v,j=sys.both(z[:sys.n]+1j*z[sys.n:]);return np.block([[j.real,-j.imag],[j.imag,j.real]])
            z0=np.r_[x0.real,x0.imag]
        else:
            def fun(z):return sys.both(z)[0]
            def jac(z):return sys.both(z)[1]
            z0=x0
        sol=least_squares(fun,z0,jac=jac,max_nfev=args.max_nfev,ftol=1e-12,xtol=1e-12,gtol=1e-12)
        x=sol.x[:sys.n]+1j*sol.x[sys.n:] if cm else sol.x
        v,j=sys.both(x);res=float(np.max(np.abs(v)));accepted=res<1e-9
        rec={'start':i,'complex_start':cm,'nfev':sol.nfev,'status':int(sol.status),'residual_inf_weighted':res,'jacobian_ranks':rank(j),'small_residual_endpoint':accepted}
        if accepted:
            endpoints+=1;new=all(np.linalg.norm(x-y)>1e-4*(1+np.linalg.norm(x)) for y in seen)
            rec['distinct_coefficient_endpoint_heuristic']=new
            if new:seen.append(x)
            rec.update(sys.diagnose(x));rec['parameters']=[[float(z.real),float(z.imag)] for z in x]
        records.append(rec)
        print(json.dumps({k:v for k,v in rec.items() if k not in ('parameters','coordinate_coefficients','det_pencil_coefficients')}),flush=True)
        payload={'lane':args.lane,'seed':args.seed,'unknowns_complex':sys.n,'equations_complex':sys.m,'starts_completed':i+1,'requested_starts':args.starts,'max_nfev':args.max_nfev,'elapsed_s':time.time()-start,'small_residual_endpoints':endpoints,'distinct_coefficient_endpoints_heuristic':len(seen),'certified_curve_count':0,'complete_complex_solution_count':None,'records':records}
        (OUT/(args.lane+'.json')).write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps({k:v for k,v in payload.items() if k!='records'}),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('lane',choices=['elliptic','single','double']);p.add_argument('--starts',type=int,default=32);p.add_argument('--real-starts',type=int,default=24);p.add_argument('--max-nfev',type=int,default=300);p.add_argument('--seed',type=int,default=20260905);run(p.parse_args())
