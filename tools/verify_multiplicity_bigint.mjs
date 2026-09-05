import fs from 'node:fs';
import assert from 'node:assert/strict';
const dir='results/multiplicity_square_2026_09_05/';
const read=p=>JSON.parse(fs.readFileSync(dir+p,'utf8'));
function mul(a,b){let c=Array(a.length+b.length-1).fill(0n);a.forEach((x,i)=>b.forEach((y,j)=>c[i+j]+=x*y));return c;}
function pow(a,n){let c=[1n];for(let k=0;k<n;k++)c=mul(c,a);return c;}
let prefixes=0;
for(const row of read('conic_crt.json').rows){
 const x=row.combined_coefficients.map(BigInt),A=[x[12],x[13],1n],B=x.slice(14,17);
 const rhs=pow(A.map((v,i)=>v+7n**6n*B[i]),6);
 const lhs=pow(A,6);
 for(let i=0;i<4;i++){const s=pow(x.slice(3*i,3*i+3),6);s.forEach((v,j)=>lhs[j]+=7n**6n*v);}
 const residual=lhs.map((v,j)=>v-rhs[j]);assert(residual.some(v=>v!==0n));
 for(const m of [8n,9n,5n,7n**47n])assert(residual.every(v=>v%m===0n));
 prefixes++;
}
let cubics=0;
for(const row of read('unequal_square_cube.json').samples){
 const C=row.C.map(BigInt),a=BigInt(row.a),b=BigInt(row.b);
 assert(C.every(x=>x>0n));assert.equal(C.slice(0,5).reduce((s,x)=>s+x**3n,0n),C[5]**3n);
 assert.equal(C[5]*b*b,C[4]*a*a);assert.equal(182n*(C[0]+C[1]),C[2]+C[3]);
 cubics++;
}
const result={status:'PASS',combined_conic_prefixes:prefixes,positive_cubic_samples:cubics,scope:'Node BigInt independent replay of retained samples and all CRT conic prefixes'};
fs.writeFileSync(dir+'node_verification.json',JSON.stringify(result,null,2)+'\n');console.log(result);
