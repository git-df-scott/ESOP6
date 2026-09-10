"""Family I sweep. For integer A=a*E1+c*E2 and Omega=w0 E1^2+w1 E1E2+w2 E2^2 compute
 H = E1^6 - E2^6 - 2*Phi(A,1,Omega), Phi = A^6+15A^4 Om+15A^2 Om^2+Om^3,
and test whether H has Waring rank <= 2 (4x4 catalecticant rank <= 2), i.e. H = lam*M^6 + mu*N^6.
Exact integer arithmetic via numpy int64 on vectorised boxes (coefficients stay < 2^63 for this box)."""
import numpy as np, itertools, sys, json
from fractions import Fraction
R=int(sys.argv[1]) if len(sys.argv)>1 else 6
W=int(sys.argv[2]) if len(sys.argv)>2 else 8
def poly_pow(P,k):  # P: (...,deg+1) int arrays, full-precision python ints would be slow; use object only for verification
    out=np.ones(P.shape[:-1]+(1,),dtype=object)
    for _ in range(k): out=poly_mul(out,P)
    return out
def poly_mul(P,Q):
    n=P.shape[-1]+Q.shape[-1]-1
    out=np.zeros(P.shape[:-1]+(n,),dtype=object)
    for i in range(P.shape[-1]):
        for j in range(Q.shape[-1]):
            out[...,i+j]+=P[...,i]*Q[...,j]
    return out
hits=[]
binom6=[1,6,15,20,15,6,1]
# sextic coefficients ordered by power of E2 : coef[k] multiplies E1^(6-k) E2^k
count=0
for a in range(-R,R+1):
  for c in range(-R,R+1):
    if a==0 and c==0: continue
    A=np.array([[a,c]],dtype=object)
    A2=poly_pow(A,2); A4=poly_pow(A,4); A6=poly_pow(A,6)
    ws=np.array(list(itertools.product(range(-W,W+1),repeat=3)),dtype=object)
    ws=ws[(ws!=0).any(axis=1)]
    Om=ws
    Om2=poly_pow(Om,2); Om3=poly_pow(Om,3)
    def padto(P,n=7):
        out=np.zeros(P.shape[:-1]+(n,),dtype=object); out[...,:P.shape[-1]]=P; return out
    Phi=padto(np.broadcast_to(A6,(len(ws),7)))+15*padto(poly_mul(np.broadcast_to(A4,(len(ws),5)),Om))+15*padto(poly_mul(np.broadcast_to(A2,(len(ws),3)),Om2))+padto(Om3)
    H=-2*Phi
    H[:,0]+=1; H[:,6]-=1
    # catalecticant: h_k = coef_k / binom(6,k); Cat[i][j]=h_{i+j}, i,j=0..3
    h=np.array([[Fraction(int(x),b) for x,b in zip(row,binom6)] for row in H],dtype=object)
    # rank<=2 iff all 3x3 minors of the 4x4 Hankel vanish; test via determinant of full 4x4 and its leading 3x3 minors is not sufficient; do exact rank
    for idx in range(len(ws)):
        hk=h[idx]
        M=[[hk[i+j] for j in range(4)] for i in range(4)]
        # exact rank by fraction gaussian elimination
        Mt=[row[:] for row in M]; rank=0
        for col in range(4):
            piv=None
            for r in range(rank,4):
                if Mt[r][col]!=0: piv=r;break
            if piv is None: continue
            Mt[rank],Mt[piv]=Mt[piv],Mt[rank]
            for r in range(4):
                if r!=rank and Mt[r][col]!=0:
                    f=Mt[r][col]/Mt[rank][col]
                    Mt[r]=[x-f*y for x,y in zip(Mt[r],Mt[rank])]
            rank+=1
        count+=1
        w0,w1,w2=[int(x) for x in ws[idx]]
        if rank<=2 and w1*w1-4*w0*w2!=0 and any(int(x)!=0 for x in H[idx][1:6]):
            hits.append(dict(a=a,c=c,omega=[int(x) for x in ws[idx]],rank=rank,H=[int(x) for x in H[idx]]))
print("tested",count,"hits",len(hits))
json.dump(hits,open("/tmp/claude-0/-home-user/15becd99-5cb2-5637-beeb-15d2ba0dfed3/scratchpad/famI_hits.json","w"))
for hh in hits[:40]: print(hh)
