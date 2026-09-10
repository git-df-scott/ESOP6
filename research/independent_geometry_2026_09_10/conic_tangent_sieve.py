#!/usr/bin/env python3
"""Exhaustive F7 necessary tangent identity sieve, not a rational lift search."""
import itertools, json, time
from collections import Counter
from pathlib import Path
import numpy as np
p=7
# Coefficients ascending in t, homogeneous degree 2 understood.
forms=[(0,0,0)]+[(1,b,c) for b in range(p) for c in range(p)]+[(0,1,c) for c in range(p)]+[(0,0,1)]
assert len(forms)==58 and len(set(forms))==58

def mul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b): out[i+j]=(out[i+j]+x*y)%p
 return out

def power(a,n):
 out=[1]
 for _ in range(n):out=mul(out,a)
 return out
weights=np.array([p**i for i in range(13)],dtype=np.int64)
def encode(a):return int(np.asarray(a,dtype=np.int64)@weights)
powers=np.array([power(a,6) for a in forms],dtype=np.int64)
pairs=list(itertools.combinations_with_replacement(range(58),2))
pairpowers=np.array([(powers[a]+powers[b])%p for a,b in pairs])
targets={}
for wi,W in enumerate(forms[1:],1):
 w5=power(W,5)
 for H in itertools.product(range(p),repeat=3):
  if H==(0,0,0):continue
  rhs=[6*x%p for x in mul(w5,H)]
  targets.setdefault(encode(rhs),[]).append((wi,H))
nontrivial=set(); trivial=set(); examples=[]; started=time.time(); hits=0
for j in range(len(pairs)):
 sums=(pairpowers[j:]+pairpowers[j])%p
 keys=sums@weights
 for off,key in enumerate(keys):
  matches=targets.get(int(key))
  if not matches:continue
  ids=tuple(sorted(pairs[j]+pairs[j+off]))
  for wi,H in matches:
   hits+=1
   item=(ids,wi,H)
   if any(i not in (0,wi) for i in ids):
    nontrivial.add(item)
    if len(examples)<12:
     examples.append({'U':[forms[i] for i in ids], 'W':forms[wi], 'H':H})
   else:trivial.add(item)
def rank7(rows):
 a=[list(row) for row in rows]; rank=0
 for col in range(len(a[0])):
  pivot=next((i for i in range(rank,len(a)) if a[i][col]%7),None)
  if pivot is None:continue
  a[rank],a[pivot]=a[pivot],a[rank]
  inv=pow(int(a[rank][col]),-1,7)
  a[rank]=[v*inv%7 for v in a[rank]]
  for i in range(len(a)):
   if i!=rank:
    c=a[i][col]
    a[i]=[(v-c*w)%7 for v,w in zip(a[i],a[rank])]
  rank+=1
 return rank
def solve7(jac,rhs):
 a=[[int(v)%7 for v in row]+[int(b)%7] for row,b in zip(jac,rhs)]; pivots=[]; r=0
 for c in range(len(jac[0])):
  pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r]; inv=pow(a[r][c],-1,7)
  a[r]=[v*inv%7 for v in a[r]]
  for i in range(len(a)):
   if i!=r:
    z=a[i][c]; a[i]=[(v-z*w)%7 for v,w in zip(a[i],a[r])]
  pivots.append(c);r+=1
 assert not any(not any(row[:-1]) and row[-1] for row in a)
 sol=[0]*len(jac[0])
 for i,c in enumerate(pivots):sol[c]=a[i][-1]
 return sol
seeds=[]
for ids,wi,H in sorted(nontrivial):
 U=[forms[i] for i in ids]; W=forms[wi]
 lhs=[sum(power(u,6)[k] for u in U)%7 for k in range(13)]
 rhs=[6*v%7 for v in mul(power(W,5),H)]
 assert lhs==rhs
 columns=[]
 bases=[[6*v%7 for v in power(u,5)] for u in U]
 bases += [[-30*v%7 for v in mul(power(W,4),H)],[-6*v%7 for v in power(W,5)]]
 for base in bases:
  for j in range(3):columns.append(([0]*j+base+[0]*(2-j)))
 jac=list(zip(*columns)); rank=rank7(jac)
 def mul_int(a,b):
  out=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):out[i+j]+=x*y
  return out
 def pow_int(a,n):
  out=[1]
  for _ in range(n):out=mul_int(out,a)
  return out
 exact=[sum(pow_int(u,6)[k] for u in U)-6*mul_int(pow_int(W,5),H)[k] for k in range(13)]
 assert all(v%7==0 for v in exact)
 compatible=rank7([list(row)+[-exact[k]//7%7] for k,row in enumerate(jac)])==rank
 lift49=None; lift343=None; next_compatible=False
 if compatible:
  delta=solve7(jac,[-v//7 for v in exact])
  lift49=[[x+7*delta[3*i+j] for j,x in enumerate(form)] for i,form in enumerate(U+[W,H])]
  uu=lift49[:4];ww,hh=lift49[4:]
  ex49=[sum(pow_int(u,6)[k] for u in uu)-6*mul_int(pow_int(ww,5),hh)[k] for k in range(13)]
  assert all(v%49==0 for v in ex49)
  next_compatible=rank7([list(row)+[-ex49[k]//49%7] for k,row in enumerate(jac)])==rank
  if next_compatible:
   delta2=solve7(jac,[-v//49 for v in ex49])
   lift343=[[x+49*delta2[3*i+j] for j,x in enumerate(form)] for i,form in enumerate(lift49)]
   uu=lift343[:4];ww,hh=lift343[4:]
   assert all((sum(pow_int(u,6)[k] for u in uu)-6*mul_int(pow_int(ww,5),hh)[k])%343==0 for k in range(13))
 seed={'U':U,'W':W,'H':H,'W_discriminant':(W[1]**2-4*W[0]*W[2])%7,
       'H_proportional_W':any(tuple(c*x%7 for x in W)==H for c in range(7)),
       'jacobian_rank_mod7':rank,'lifts_to_mod49':compatible,'canonical_mod49_lift':lift49,
       'canonical_mod49_branch_lifts_to_mod343':next_compatible,'canonical_mod343_lift':lift343}
 seeds.append(seed)
Path(__file__).with_name('conic_tangent_seeds.json').write_text(json.dumps(seeds,indent=2)+'\n')
result={'field':7,'projective_quadratic_forms_including_zero':58,
 'unordered_pair_count':len(pairs),'unordered_pair_pair_iterations':len(pairs)*(len(pairs)+1)//2,
 'unique_rhs_polynomials':len(targets),'raw_matching_presentations':hits,
 'unique_trivial_solutions':len(trivial),'unique_nonproportional_solutions':len(nontrivial),
 'jacobian_rank_histogram':dict(Counter(x['jacobian_rank_mod7'] for x in seeds)),
 'W_discriminant_histogram':dict(Counter(x['W_discriminant'] for x in seeds)),
 'canonical_mod49_branches_lifting_to_mod343_count':sum(x['canonical_mod49_branch_lifts_to_mod343'] for x in seeds),
 'lifts_to_mod49_count':sum(x['lifts_to_mod49'] for x in seeds),
 'H_proportional_W_count':sum(x['H_proportional_W'] for x in seeds),
 'nonproportional_examples':examples[:1],'elapsed_seconds':time.time()-started,
 'scope':'Exhaustive sum of four sixth powers = 6 H W^5 over F7 with quadratic forms; modular solutions are not Q7 or rational lifts.'}
Path(__file__).with_name('conic_tangent_sieve.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
