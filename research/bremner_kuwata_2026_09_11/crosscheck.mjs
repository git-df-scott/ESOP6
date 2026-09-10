// Independent BigInt re-verification of the Bremner/Kuwata mechanism.
const q = (a,b) => a*a + a*b - b*b;
let fails = 0;
const chk = (name, cond) => { if(!cond){fails++; console.log("FAIL "+name);} else console.log("PASS "+name); };

// 1. the engine identity  q(s,t)^3 - q(t,s)^3 = 2(s^6 - t^6)
let bad = 0;
for (let i=0;i<4000;i++){
  const s = BigInt(Math.floor(Math.random()*20001)-10000);
  const t = BigInt(Math.floor(Math.random()*20001)-10000);
  if (q(s,t)**3n - q(t,s)**3n !== 2n*(s**6n - t**6n)) bad++;
}
chk("q(s,t)^3 - q(t,s)^3 = 2(s^6-t^6)  (4000 random integer pairs)", bad===0);

// 2. the same-sign version needs i:  q(s,it)^3 - q(it,s)^3 = 2(s^6+t^6).
//    Over Z we verify the equivalent real identity
//    (s^2+t^2)^3 ... use the explicit expansion of q(s,it)=s^2+i s t+t^2:
//    (s^2+st i+t^2)^3 + (s^2-st i+t^2)^3 = 2(s^6+t^6)  <=>
//    2*[ (s^2+t^2)^3 - 3 s^2 t^2 (s^2+t^2) ] = 2(s^6+t^6)
bad = 0;
for (let i=0;i<4000;i++){
  const s = BigInt(Math.floor(Math.random()*20001)-10000);
  const t = BigInt(Math.floor(Math.random()*20001)-10000);
  const A = (s*s+t*t)**3n - 3n*s*s*t*t*(s*s+t*t);
  if (2n*A !== 2n*(s**6n + t**6n)) bad++;
}
chk("real part of the Q(i) same-sign device is correct (4000 pairs)", bad===0);

// 3. (1.3) => (1.1) and (1.2) on the 20 census points found on K_B
const PTS = [[3,19,22,-23,10,-15],[15,52,65,37,67,36],[23,54,73,-74,47,-33],
 [3,55,80,-32,81,43],[11,65,78,50,81,37],[40,125,129,-113,-51,-136],
 [1,132,133,147,71,-92],[26,169,225,111,230,121],[14,163,243,75,245,142],
 [11,188,243,-148,249,103],[29,197,261,131,267,139],[36,179,275,65,276,169],
 [113,241,282,186,293,173],[37,199,309,-99,311,173],[23,282,311,-326,107,-243],
 [92,311,317,-277,-124,-337],[27,317,356,271,372,127],[148,299,507,-177,508,281],
 [93,409,512,293,528,271],[1,500,515,556,197,-409]];
let n=0;
for (const P of PTS){
  const [x,y,z,u,v,w] = P.map(BigInt);
  const e1 = q(x,u)-q(w,z), e2 = q(y,v)-q(u,x), e3 = q(z,w)-q(v,y);
  const s6 = x**6n+y**6n+z**6n-u**6n-v**6n-w**6n;
  const s2 = x*x+y*y+z*z-u*u-v*v-w*w;
  if (e1===0n && e2===0n && e3===0n && s6===0n && s2===0n) n++;
  else console.log("  offending", P, [e1,e2,e3,s6,s2].map(String));
}
chk(`all ${PTS.length} census points satisfy (1.3), (1.2), (1.1)`, n===PTS.length);

// 4. the (2,1,5) multigrade has no positive solutions:  sum a^6 < (sum a^2)^3
bad = 0;
for (let i=0;i<200000;i++){
  const a = Array.from({length:5},()=>BigInt(1+Math.floor(Math.random()*100000)));
  const s2 = a.reduce((A,t)=>A+t*t,0n);
  const s6 = a.reduce((A,t)=>A+t**6n,0n);
  if (s6 >= s2**3n) bad++;
}
chk("sum a_i^6 < (sum a_i^2)^3 for 200000 random positive 5-tuples", bad===0);
process.exit(fails?1:0);
