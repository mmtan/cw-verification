// Exact verification for Character and Multiplier Obstructions for
// Circulant Weighing Matrices. Reconstructed from the stated finite systems.
// Integer arithmetic only. No original repository implementation is used.
// See README.md for the supported cases, integer bounds, and commands.
#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <limits>
using namespace std;
using U64=uint64_t;
static_assert(numeric_limits<int>::digits >= 31, "At least 32-bit int is required.");
static_assert(numeric_limits<U64>::digits == 64, "A 64-bit counter type is required.");
void check(bool v,const string& msg){if(!v)throw runtime_error(msg);}
template<class T> void pv(const vector<T>&a){for(auto x:a)cout<<x<<' ';cout<<'\n';}
struct Orbits{
 int n,t; vector<int> reps,w,id;
 Orbits(int n_,int t_):n(n_),t(t_),id(n_,-1){
  check(gcd(n,t)==1,"multiplier not invertible");
  for(int a=0;a<n;a++)if(id[a]<0){int k=w.size(),x=a,c=0;reps.push_back(a);
   do{id[x]=k;c++;x=(x*t)%n;}while(x!=a);w.push_back(c);
  }
 }
};
// For integers, |a|<=a^2. Thus an energy bound N also bounds |sum| by N.
// dp[i][e][s+N] exactly records suffix feasibility, not a heuristic bound.
struct EnumReal{
 const vector<int>& w; int b,N,S; vector<vector<vector<uint8_t>>> dp;
 EnumReal(const vector<int>&ww,int bb,int NN):w(ww),b(bb),N(NN),S(2*NN+1){
  int k=w.size();dp.assign(k+1,vector<vector<uint8_t>>(N+1,vector<uint8_t>(S)));
  dp[k][0][N]=1;
  for(int i=k-1;i>=0;i--)for(int e=0;e<=N;e++)for(int s=-N;s<=N;s++)if(dp[i+1][e][s+N])
   for(int v=-b;v<=b;v++){int ee=e+w[i]*v*v,ss=s+w[i]*v;if(ee<=N&&abs(ss)<=N)dp[i][ee][ss+N]=1;}
 }
 bool feasible(int i,int e,int s)const{return e>=0&&e<=N&&abs(s)<=N&&dp[i][e][s+N];}
 void run(int sum,int energy,const function<void(const vector<int>&)>& fn){
  vector<int>a(w.size());
  function<void(int,int,int)> dfs=[&](int i,int e,int s){
   if(!feasible(i,e,s))return;
   if(i==(int)w.size()){fn(a);return;}
   for(int v=-b;v<=b;v++){a[i]=v;dfs(i+1,e-w[i]*v*v,s-w[i]*v);}
  };dfs(0,energy,sum);
 }
};
struct Term{int a,b,m;};
vector<Term> realterms(const Orbits&o,int h){
 map<pair<int,int>,int> mm;
 for(int j=0;j<o.n;j++){int a=o.id[j],b=o.id[(j+h)%o.n];if(a>b)swap(a,b);mm[{a,b}]++;}
 vector<Term> tt;for(auto kv:mm)tt.push_back({kv.first.first,kv.first.second,kv.second});return tt;
}
int corr(const vector<Term>&tt,const vector<int>&a){int s=0;for(auto t:tt)s+=t.m*a[t.a]*a[t.b];return s;}
void realtest(int n,int t,int b,int row,int energy,vector<int> shifts,vector<U64> expected={},bool printlast=false){
 Orbits o(n,t);cout<<"REAL n="<<n<<" t="<<t<<" bound="<<b<<" row="<<row<<" energy="<<energy<<'\n';
 cout<<"reps ";pv(o.reps);cout<<"sizes ";pv(o.w);
 vector<vector<Term>>ts;for(int h:shifts)ts.push_back(realterms(o,h));
 vector<U64>counts(shifts.size()+1);vector<vector<int>>last,full;
 EnumReal en(o.w,b,energy);
 en.run(row,energy,[&](const vector<int>&a){counts[0]++;
  for(int i=0;i<(int)shifts.size();i++){
   if(printlast&&i==(int)shifts.size()-1)last.push_back(a);
   if(corr(ts[i],a)) return;
   counts[i+1]++;
  }if(full.size()<10000)full.push_back(a);
 });
 cout<<"shifts ";pv(shifts);cout<<"counts ";pv(counts);
 if(!expected.empty())check(counts==expected,"real filtering mismatch");
 if(printlast){cout<<"Vectors before final shift and final correlation:\n";for(auto a:last){for(int v:a)cout<<v<<' ';cout<<" | "<<corr(ts.back(),a)<<'\n';}}
 if(counts.back()<=30){cout<<"Full survivors\n";for(auto a:full)pv(a);}
 cout<<'\n';
}
void gausstest(int n,int energy=64,int b=2,int row=8,vector<int>shifts={},vector<U64>expected={}){
 Orbits o(n,2);cout<<"GAUSSIAN n="<<n<<" bound="<<b<<" row="<<row<<" energy="<<energy<<'\n';
 cout<<"reps ";pv(o.reps);cout<<"sizes ";pv(o.w);
 if(shifts.empty())for(int h=1;h<=n/2;h++)shifts.push_back(h);
 vector<vector<Term>>ts;
 for(int h:shifts){map<pair<int,int>,int>mm;for(int j=0;j<n;j++)mm[{o.id[(j+h)%n],o.id[j]}]++;
  vector<Term>terms;for(auto p:mm)terms.push_back({p.first.first,p.first.second,p.second});ts.push_back(terms);}
 EnumReal en(o.w,b,energy);vector<U64> counts(shifts.size()+1);
 vector<vector<pair<int,int>>>full;
 for(int e=0;e<=energy;e++){
  if(!en.feasible(0,e,row)||!en.feasible(0,energy-e,0))continue;
  vector<vector<int>>rr,ii;
  en.run(row,e,[&](const vector<int>&a){rr.push_back(a);});
  en.run(0,energy-e,[&](const vector<int>&a){ii.push_back(a);});
  for(auto&r:rr)for(auto&s:ii){counts[0]++;bool ok=true;
   for(int h=0;h<(int)shifts.size();h++){int re=0,im=0;
    for(auto z:ts[h]){re+=z.m*(r[z.a]*r[z.b]+s[z.a]*s[z.b]);im+=z.m*(s[z.a]*r[z.b]-r[z.a]*s[z.b]);}
    if(re||im){ok=false;break;}counts[h+1]++;
   }
   if(ok&&full.size()<10000){vector<pair<int,int>>a;for(int i=0;i<(int)r.size();i++)a.push_back({r[i],s[i]});full.push_back(a);}
  }
 }
 cout<<"shifts ";pv(shifts);cout<<"counts ";pv(counts);
 if(!expected.empty())check(counts==expected,"Gaussian filtering mismatch");
 if(counts.back()<=30){cout<<"Full survivors\n";for(auto a:full){for(auto v:a)cout<<'('<<v.first<<','<<v.second<<") ";cout<<'\n';}}
 cout<<'\n';
}
vector<int> mul3(const vector<int>&a,const vector<int>&b){vector<int>c(a.size()+b.size()-1);for(int i=0;i<(int)a.size();i++)for(int j=0;j<(int)b.size();j++)c[i+j]=(c[i+j]+a[i]*b[j]+9)%3;return c;}
string canon(const vector<int>&a){string best(35,'3');for(int sign=1;sign<=2;sign++)for(int h=0;h<35;h++){string w;for(int j=0;j<35;j++)w+=char('0'+a[(j+h)%35]*sign%3);best=min(best,w);}return best;}
void eisenstein(){
 vector<int>f={1,-1,-1,1,-1,1,0,1,-1,0,1,0,1};vector<int>fr(f.rbegin(),f.rend());
 vector<int>g=mul3(mul3(mul3({-1,1},{1,1,1,1,1}),{1,1,1,1,1,1,1}),f);
 auto prod=mul3(g,fr);for(int j=0;j<36;j++)check(prod[j]==(j==0?2:j==35?1:0),"F3 factorization mismatch");
 map<int,U64> dist;map<string,int>classes;
 for(int v=0;v<531441;v++){int q=v;vector<int>a(35);for(int i=0;i<12;i++){int c=q%3;q/=3;for(int j=0;j<24;j++)a[i+j]=(a[i+j]+c*g[j])%3;}
  int wt=0;for(int c:a)wt+=c!=0;dist[wt]++;if(wt==12)classes[canon(a)]++;
 }
 cout<<"TERNARY DISTRIBUTION\n";for(auto [w,n]:dist)cout<<w<<' '<<n<<'\n';
 check(dist==map<int,U64>{{0,1},{12,420},{15,2520},{18,37590},{21,158550},{24,218610},{27,102620},{30,11130}},"code weight distribution mismatch");
 check(classes.size()==6,"code class count mismatch");
 // Appendix A representatives, retaining the order in the manuscript.
 vector<vector<int>> signedpositions={{9,13,-14,-15,19,20,22,-23,-27,28,-33,-34},{7,10,12,-14,16,-17,-19,24,-27,-30,-31,34},{6,-13,14,18,-19,-21,24,26,27,-31,-32,-34},{6,8,12,13,-18,-22,24,-27,-31,32,-33,-34},{6,-8,-12,-14,19,21,22,24,-26,-31,33,-34},{6,7,-10,-11,-14,-16,-20,25,-27,30,31,34}};
 set<string> seen;
 vector<int>lastshift={6,6,4,8,5,8},lastcount={6,45,324,2,72,2};
 for(int cl=0;cl<6;cl++){
  vector<int>a(35),pos,sign;for(int sp:signedpositions[cl]){int p=abs(sp);a[p]=sp>0?1:2;pos.push_back(p);sign.push_back(sp>0?1:-1);}
  string key=canon(a);check(classes.count(key)&&classes[key]==70,"Appendix A class not found");seen.insert(key);
  vector<vector<pair<int,int>>> pairs(18);
  for(int u=0;u<12;u++)for(int v=0;v<12;v++){int h=(pos[u]-pos[v]+35)%35;if(h>=1&&h<=17)pairs[h].push_back({u,v});}
  vector<U64>count(18);int phases[12];int roots[3][2]={{1,0},{0,1},{-1,-1}};
  for(int code=0;code<177147;code++){int q=code;phases[0]=0;for(int j=1;j<12;j++){phases[j]=q%3;q/=3;}count[0]++;
   for(int h=1;h<=17;h++){int aa=0,bb=0;for(auto [u,v]:pairs[h]){int power=(phases[u]-phases[v]+3)%3;aa+=sign[u]*sign[v]*roots[power][0];bb+=sign[u]*sign[v]*roots[power][1];}if(aa||bb)break;count[h]++;}
  }
  cout<<"EISENSTEIN CLASS "<<cl+1<<" stages ";pv(count);
  check(count[lastshift[cl]]==U64(lastcount[cl])&&count[lastshift[cl]+1]==0,"Eisenstein terminal filtering mismatch");
 }
 check(seen.size()==6,"Appendix representatives not distinct");
}
void localfour(){int tot=0,n3=0;for(int a=-1;a<=1;a++)for(int b=-1;b<=1;b++)for(int c=-1;c<=1;c++)for(int d=-1;d<=1;d++){tot++;if(abs(a+b+c+d)==3){n3++;check(abs(a-b+c-d)==1,"C4 local character failure");}}cout<<"C4 local fibers tested "<<tot<<"; |principal sum|=3 fibers "<<n3<<"\n";}
int main(int argc,char**argv){try{
 if(argc>1){string mode=argv[1];
  check((mode=="real"&&argc==7)||(mode=="gauss"&&argc==6),
        "Usage: verify_cw [real n multiplier bound sum energy | gauss n energy bound sum]");
  if(mode=="real"){int n=stoi(argv[2]),t=stoi(argv[3]),b=stoi(argv[4]),s=stoi(argv[5]),e=stoi(argv[6]);vector<int>h;for(int j=1;j<=n/2;j++)h.push_back(j);realtest(n,t,b,s,e,h);return 0;}if(mode=="gauss"){gausstest(stoi(argv[2]),stoi(argv[3]),stoi(argv[4]),stoi(argv[5]));return 0;}}
 vector<int>h35;for(int h=1;h<=17;h++)h35.push_back(h);
 realtest(35,4,3,6,36,h35,{1434,163,14,14,14,6,6,2,2,2,2,2,2,2,2,2,2,2});
 realtest(35,4,4,6,36,h35,{1600,189,20,20,20,8,8,2,2,2,2,2,2,2,2,2,2,2});
 eisenstein();localfour();
 gausstest(35,64,2,8,{1,2,3,4,5},{1152,4,4,4,4,0});
 gausstest(45,64,2,8,{1,2,3},{58188,1242,1242,0});
 gausstest(49,64,2,8,{1},{32,0});
 realtest(20,7,6,7,49,{5,10},{716,76,0});
 realtest(116,7,1,7,49,{1,2,3,4},{6088,624,40,16,0},true);
 realtest(64,7,3,7,49,{1,2,3,4,5,6,7,8,11,12,16,24},{22880810,5952866,455924,413640,14670,14670,7836,7836,338,338,22,4,0},true);
 cout<<"ALL ASSERTED INDEPENDENT CHECKS PASSED\n";
}catch(exception&e){cerr<<e.what()<<'\n';return 1;}}
