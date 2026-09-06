#!/usr/bin/env python3
"""Exact controls, F7 degree-four certificate, and plane-section ideals."""
import itertools,json,math
from pathlib import Path
import sympy as S
OUT=Path('results/astra_curves_2026_09_05');OUT.mkdir(exist_ok=True)
u,v,w,s,t=S.symbols('u v w s t');a,d=S.symbols('alpha delta')

def conic_controls():
    A=u*u+v*v
    L=[a*u,a*v,u+v,u-v,w,d*w]
    F=S.expand(sum(x**6 for x in L[:5])-L[5]**6)
    F=S.rem(S.Poly(F,a),S.Poly(a*a-2,a)).as_expr()
    F=S.rem(S.Poly(F,d),S.Poly(d**6-11,d)).as_expr()
    Q=A-w*w;R=10*(A*A+A*w*w+w**4)
    assert S.expand(F-Q*R)==0
    P=[s*s-t*t,2*s*t,s*s+2*s*t-t*t,s*s-2*s*t-t*t,s*s+t*t,s*s+t*t]
    weights=[8,8,1,1,1,-11]
    assert S.expand(sum(c*p**6 for c,p in zip(weights,P)))==0
    columns=[]
    for c,p in zip(weights,P):
        for mon in [s*s,s*t,t*t]:
            h=S.Poly(6*c*p**5*mon,s,t)
            columns.append([h.coeff_monomial(s**(12-k)*t**k) for k in range(13)])
    J=S.Matrix.hstack(*[S.Matrix(c) for c in columns])
    jac_rank=J.rank()
    dcolumns=[]
    for i,(c,p) in enumerate(zip(weights,P)):
        mm=[s*s,s*t,t*t]
        if i==0:mm=[s*t,t*t]
        if i==4:mm=[s*s,t*t]
        if i==5:mm=[s*s+t*t]
        for mon in mm:
            h=S.Poly(6*c*p**5*mon,s,t)
            dcolumns.append(S.Matrix([h.coeff_monomial(s**(12-k)*t**k) for k in range(13)]))
    d_rank=S.Matrix.hstack(*dcolumns).rank()
    # The factorization chart after the exact diagonal field extension.
    cs=S.symbols('l0:9');qs=S.symbols('q0:5');rs=S.symbols('r0:15')
    lm=[cs[3*i]*u+cs[3*i+1]*v+cs[3*i+2]*w for i in range(3)]
    qmons=[u*u,u*v,u*w,v*v,v*w];rmons=[u**i*v**j*w**(4-i-j) for i in range(4,-1,-1) for j in range(4-i,-1,-1)]
    q=sum(c*m for c,m in zip(qs,qmons))-w*w;r=sum(c*m for c,m in zip(rs,rmons))
    ff=8*u**6+8*v**6+sum(z**6 for z in lm)-11*w**6
    eq=S.Poly(ff-q*r,u,v,w)
    mon6=[u**i*v**j*w**(6-i-j) for i in range(6,-1,-1) for j in range(6-i,-1,-1)]
    eqs=[eq.coeff_monomial(m) for m in mon6]
    vs=list(cs)+list(qs)+list(rs)
    point=list(map(S.Integer,[1,1,0,1,-1,0,0,0,1]))+list(map(S.Integer,[1,0,0,1,0]))+[S.Poly(R,u,v,w).coeff_monomial(m) for m in rmons]
    sub=dict(zip(vs,point));assert all(e.subs(sub)==0 for e in eqs)
    jr=S.Matrix(eqs).jacobian(vs).subs(sub).rank()
    # Original X: first, second and sixth linear coordinates fixed u,v,w.
    original=S.Poly(u**6+v**6+sum(z**6 for z in lm)-w**6-q*r,u,v,w)
    (OUT/'conic_plane_ideal.json').write_text(json.dumps({'chart':'x1=u,x2=v,x6=w; q[w^2]=-1','variables':[str(z) for z in vs],'equations':[str(original.coeff_monomial(m)) for m in mon6],'unknowns':29,'equations_count':28,'conditions':'rank of plane embedding is 3 automatically; det(q)!=0; all Grassmann charts and conic coefficient charts required for global coverage'},indent=2)+'\n')
    return {'exact_factorization':str(S.factor(F)),'parameter_identity_zero':True,'field_relations':['alpha^2=2, alpha>0','delta^6=11, delta>0'],'coordinate_forms':[str(z) for z in P],'coordinate_scales':['alpha','alpha','1','1','1','delta'],'conic_parameter_jacobian_rank':jac_rank,'conic_parameter_unknowns':18,'conic_parameter_equations':13,'normalized_D_coordinate_scales':['1','1','1/alpha','1/alpha','1/alpha','delta/alpha'],'normalized_D_unknowns':14,'normalized_D_equations':13,'normalized_D_jacobian_rank':d_rank,'plane_factorization_jacobian_rank_over_field':jr,'plane_factorization_unknowns':29,'plane_factorization_equations':28,'p6_definite':True,'all_coordinates_definite':False,'positive_real_specialization':'s=1,t=1/4','defined_over_Q':False,'rational_point_claim':False}

def mul(p,q,mod=7):
    r=[0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]=(r[i+j]+a*b)%mod
    return r
def pow6(p):
    r=[1]
    for _ in range(6):r=mul(r,p)
    return tuple(r)
def quartic_F7():
    # Homogeneous coefficients ordered s^4,...,t^4. One projective representative each.
    forms=[]
    for first in range(5):
        for tail in itertools.product(range(7),repeat=4-first):forms.append((0,)*first+(1,)+tail)
    assert len(forms)==2801
    sixth={pow6(p):p for p in forms};assert len(sixth)==2801
    points=[(1,k) for k in range(7)]+[(0,1)]
    def support(p):return tuple(i for i,(s,t) in enumerate(points) if sum(p[k]*pow(s,4-k)*pow(t,k) for k in range(5))%7)
    by_support={support(p):p for p in forms if len(support(p))==4}
    assert len(by_support)==70
    records=[]
    for supp,p in sorted(by_support.items()):
        comp=tuple(i for i in range(8) if i not in supp)
        if supp>comp:continue
        q=by_support[comp];summ=tuple((a+b)%7 for a,b in zip(pow6(p),pow6(q)))
        assert summ not in sixth
        records.append({'support':supp,'complement':comp,'p':p,'q':q,'sixth_sum':summ,'is_sixth_of_quartic':False})
    assert len(records)==35
    return {'field':7,'degree':4,'nonzero_projective_forms':len(forms),'four_point_support_forms':70,'complementary_support_pairs':35,'matching_rhs_sixth_powers':0,'coefficient_identity_solutions':'one nonzero LHS polynomial proportional to RHS; all others zero, or all zero','proof':'Supports of nonzero quartics have size >=4. Distinct nonzero left forms have disjoint supports because sixth powers are 0 or 1 and five terms cannot sum to 7. At most two survive. If two, each support has size 4 and they partition P1(F7); the 35 exhaustive complement pairs have no sixth-power RHS. If one, equality of sixth powers in F7(s/t) implies constant proportionality.','records':records}

def paired_planes():
    # Symbolic representative valid over Q(i); roots of z^6=-1 provide all six choices.
    h=S.symbols('h');F=u**6+(h*u)**6+v**6+(h*v)**6+w**6-w**6
    assert S.rem(S.Poly(F,h),S.Poly(h**6+1,h)).as_expr()==0
    triples=set(itertools.product(range(1,12,2),range(1,12,2),range(0,12,2)))
    orbits=[]
    while triples:
        z=min(triples);orb={tuple(k*a%12 for a in z) for k in (1,5,7,11)}
        assert orb<=triples;triples-=orb;orbits.append(sorted(orb))
    degrees=[len(o) for o in orbits]
    assert degrees.count(2)==4 and degrees.count(4)==52
    return {'pairings':15,'root_choices_each_pair':6,'distinct_standard_complex_planes':15*6**3,'conics_in_each_plane_dimension':5,'real_standard_planes':0,'representative':['u','i*u','v','i*v','w','w'],'plane_parameter_ideal_one_pairing':['a^6+1','b^6+1','c^6-1'],'Galois_orbits_one_pairing':orbits,'degree_2_closed_plane_points_total':15*4,'degree_4_closed_plane_points_total':15*52,'residue_fields':['Q(i) for degree 2','Q(zeta_12) for degree 4'],'reason_not_real':'two pairs among positive slots require sixth root of -1, with no real choices','complete_plane_classification_claim':False,'complete_conic_component_classification_claim':False}

if __name__=='__main__':
    c=conic_controls();f=quartic_F7();p=paired_planes()
    (OUT/'quartic_F7_certificate.json').write_text(json.dumps(f,indent=2)+'\n')
    data={'conic_control':c,'paired_planes':p,'quartic_F7':{k:v for k,v in f.items() if k!='records'},'plane_cubic':{'unknowns':28,'equations':28,'real_plane_cubic_curves':0,'rational_plane_cubic_curves':0,'distinct_complex_solutions':None,'reason':'The hyperplane x6=0 has no real points on X. An odd-degree projective real curve intersects any hyperplane not containing it in a conjugation-invariant effective divisor of odd degree, hence with a real point. If contained, it has real points by odd-degree real linear section.','method':'proof; no numerical root search'},'elliptic_quartic_dimension_count':{'P3_in_P5':8,'quadric_pencil_Gr2_10':16,'restrictions_h0_OE6':24,'expected_dimension':0}}
    (OUT/'exact_results.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data,indent=2))
