"""Two square ratios by construction on a cubic tangent section.

B=(1,-1,1,-1,2,2), D1+D2=D3+D4=-1, D5=0,D6=-1/2.
Impose C6/C5=(a/b)^2. The remaining conic is X^2+Y^2=2K,
K=-a^4+3a^2*b^2-b^4. Rotate an integer norm representation rationally.
"""
import argparse,json,math,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--r-bound',type=int,default=60);p.add_argument('--rotation-bound',type=int,default=60);p.add_argument('--unequal-seed',action='store_true');arg=p.parse_args()
out=Path('results/fourfold_routes_2026_09_05');start=time.monotonic()
sa,sb,sc=(1,2,3) if arg.unequal_seed else (1,1,2)
h=sa*sa+sb*sb
seeds=[];conic_tests=0
for b in range(1,arg.r_bound+1):
 for a in range(b+1,2*b):
  if math.gcd(a,b)!=1:continue
  K=(36*h**3*a*a*b*b+(12*h**3-6*sc**6)*(a*a-b*b)**2) if arg.unequal_seed else 2*(-a**4+3*a*a*b*b-b**4)
  if K<=0:continue
  conic_tests+=1
  for X in range(math.isqrt(K)+1):
   yy=K-X*X;Y=math.isqrt(yy)
   if Y*Y==yy:
    seeds.append((a,b,X,Y));break
counts={'r_values_tested':conic_tests,'r_values_with_integer_norm_representation':len(seeds),
 'rotations':0,'positive':0,'exact_cubic_checks':0,'forced_square_ratio_checks':0,'all_remaining_squares':0}
best=0;bestC=[];sample=[];sign=(1,1,1,1,1,-1);maxC=0
for a,b,X,Y in seeds:
 d=a*a-b*b
 for n in range(1,arg.rotation_bound+1):
  for m in range(-arg.rotation_bound,arg.rotation_bound+1):
   if math.gcd(m,n)!=1:continue
   counts['rotations']+=1
   E=m*m+n*n
   xx=X*(m*m-n*n)-2*Y*m*n; yy=2*X*m*n+Y*(m*m-n*n)
   # V=(b^2*E+xx)/(2*d*E), similarly R. Clear all denominators in D.
   if arg.unequal_seed:
    V=6*h*sa*b*b*E+xx;R=6*h*sb*b*b*E+yy;den=3*sc**3*d*E
    D=(sc*sc*(-den+V),sc*sc*(-den-V),sc*sc*(-den+R),sc*sc*(-den-R),0,-2*h*den)
   else:
    V=b*b*E+xx;R=b*b*E+yy;den=2*d*E
    D=(-den+V,-den-V,-den+R,-den-R,0,-den)
   gd=math.gcd(*D);D=tuple(v//gd for v in D)
   B=(sa,-sa,sb,-sb,sc,sc)
   f=sum(e*v**3 for e,v in zip(sign,D));k=sum(e*u*v*v for e,u,v in zip(sign,B,D))
   C=tuple(f*u-3*k*v for u,v in zip(B,D))
   if C[-1]<0:C=tuple(-v for v in C)
   if min(C)<=0:continue
   counts['positive']+=1;gc=math.gcd(*C);C=tuple(v//gc for v in C)
   assert sum(v**3 for v in C[:-1])==C[-1]**3;counts['exact_cubic_checks']+=1
   assert C[-1]*b*b==C[-2]*a*a;counts['forced_square_ratio_checks']+=1
   score=sum(math.isqrt(v*C[-2])**2==v*C[-2] for v in C[:4])
   if score>best:best=score;bestC=C
   maxC=max(maxC,C[-1])
   if len(sample)<64:sample.append({'r':[a,b],'C':[str(v) for v in C]})
   if score==4:
    assert all(math.isqrt(v)**2==v for v in C)
    sol=[math.isqrt(v) for v in C]
    assert sum(v**6 for v in sol[:-1])==sol[-1]**6
    counts['all_remaining_squares']+=1
    (out/'CANDIDATE.json').write_text(json.dumps(sol)+'\n')
    print('EXACT HIT',sol,flush=True);break
report={'seed':[sa,sb,sc],'r_bound':arg.r_bound,'rotation_bound':arg.rotation_bound,'counts':counts,
 'rational_conic_seeds':seeds,'best_of_four_remaining_square_ratios':best,
 'best_primitive_cube_bases':[str(x) for x in bestC],'largest_primitive_cube_base':str(maxC),
 'samples':sample,'elapsed_seconds':time.monotonic()-start}
(out/('two_squares_unequal.json' if arg.unequal_seed else 'two_squares.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('samples','rational_conic_seeds')},indent=2))
