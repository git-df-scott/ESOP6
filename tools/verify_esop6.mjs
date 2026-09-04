// Independent standalone verifier. No packages; JavaScript BigInt only.
const args = process.argv.slice(2);
if (args.length !== 6 || args.some(x=>!/^[-+]?\d+$/.test(x))) {
  throw new Error('Usage: node verify_esop6.mjs a b c d e f');
}
const values = args.map(BigInt);
function pow6(x) { let y=1n; for(let i=0;i<6;i++) y*=x; return y; }
function gcd(a,b) { a=a<0n?-a:a; b=b<0n?-b:b;
  while(b) {const r=a%b;a=b;b=r;} return a; }
const powers=values.map(pow6), lhs=powers.slice(0,5).reduce((a,b)=>a+b,0n);
const rhs=powers[5], g=values.reduce(gcd,0n);
let normalized=values.slice(0,5).sort((a,b)=>a<b?-1:a>b?1:0).concat(values[5]);
if(g) normalized=normalized.map(x=>x/g);
const np=normalized.map(pow6), nl=np.slice(0,5).reduce((a,b)=>a+b,0n);
const positive=values.every(x=>x>0n);
const result={integers:values.map(String),positivity:positive,
  sixth_powers:powers.map(String),lhs:String(lhs),rhs:String(rhs),equality:lhs===rhs,
  gcd:String(g),sorted_normalized_tuple:normalized.map(String),
  normalized_sixth_powers:np.map(String),normalized_lhs:String(nl),
  normalized_rhs:String(np[5]),normalized_equality:nl===np[5],
  solution:positive&&lhs===rhs};
process.stdout.write(JSON.stringify(result,null,2)+'\n');
process.exitCode=result.solution?0:1;
