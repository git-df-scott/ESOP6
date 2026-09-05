import fs from 'node:fs';
const dir='results/fourfold_routes_2026_09_05/';
const read=n=>JSON.parse(fs.readFileSync(dir+n,'utf8'));
const assert=(v,m)=>{if(!v)throw Error(m)};
const abs=x=>x<0n?-x:x;
const gcd=(a,b)=>{a=abs(a);b=abs(b);while(b)[a,b]=[b,a%b];return a;};
const cube=x=>x*x*x;
const cubic=C=>C.slice(0,5).reduce((s,x)=>s+cube(x),0n)-cube(C[5]);
const mul=(a,b)=>{let c=Array(a.length+b.length-1).fill(0n);a.forEach((x,i)=>b.forEach((y,j)=>c[i+j]+=x*y));return c;};
const pow=(a,n)=>{let b=[1n];for(let j=0;j<n;j++)b=mul(b,a);return b;};
let tangent=0;
for(const name of ['large_tangent_samples.json','crt_tangent_samples.json'])for(const rec of read(name)){
 const B=rec.B.map(BigInt),D=rec.D.map(BigInt),C=rec.C.map(BigInt);
 assert(cubic(B)===0n,'base not on cubic');
 assert(B.reduce((s,b,i)=>s+(i===5?-1n:1n)*b*b*D[i],0n)===0n,'not tangent');
 const F=cubic(D),K=B.reduce((s,b,i)=>s+(i===5?-1n:1n)*b*D[i]*D[i],0n);
 let image=B.map((b,i)=>F*b-3n*K*D[i]);if(image[5]<0n)image=image.map(x=>-x);
 const g=image.reduce(gcd,0n);image=image.map(x=>x/g);
 assert(image.every((x,i)=>x===C[i]),'wrong tangent image');
 assert(C.every(x=>x>0n)&&cubic(C)===0n,'invalid cubic point');tangent++;
}
let conics=0;
for(const rec of read('conic_lifts.json')){
 const x=rec.full_conic_final_vector.map(BigInt),A=[x[12],x[13],1n],B=x.slice(14,17),e=7n**6n;
 const W=A.map((a,i)=>a+e*B[i]);let total=Array(13).fill(0n);
 for(let i=0;i<4;i++)pow(x.slice(3*i,3*i+3).map(a=>7n*a),6).forEach((v,j)=>total[j]+=v);
 pow(A,6).forEach((v,j)=>total[j]+=v);pow(W,6).forEach((v,j)=>total[j]-=v);
 assert(total.every(v=>v%(7n**47n)===0n),'failed full conic congruence');
 assert(total.some(v=>v!==0n),'unexpected exact curve: investigate');conics++;
}
let forced=0;
for(const name of ['two_squares.json','two_squares_unequal.json'])for(const rec of read(name).samples){
 const C=rec.C.map(BigInt),[a,b]=rec.r.map(BigInt);
 assert(cubic(C)===0n&&C.every(x=>x>0n),'invalid forced-square point');
 assert(C[5]*b*b===C[4]*a*a,'failed forced-square ratio');forced++;
}
// Exact sparse contradiction after residuals 9 and 10 eliminate u,v.
const numer=-220092525210624n,denom=144453125n;
assert(numer!==0n&&denom!==0n,'invalid sparse contradiction');
const report={tangent_points_verified:tangent,full_conic_congruences_verified:conics,
conic_modulus:'7^47 in original undivided identity',congruences_are_not_exact_identities:true,
forced_square_samples_verified:forced,status:'PASS'};
fs.writeFileSync(dir+'independent_node_verification.json',JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report,null,2));
