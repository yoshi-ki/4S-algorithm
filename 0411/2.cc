#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <functional>
#include <cmath>
#include <cstring>
#include <numeric>
#include <random>
#include <map>
#include <unordered_map>
#include <queue>
#include <regex>
#include <complex>
#include <list>
#include <cassert>
#include <deque>
#include <iomanip>


#define addf(T) [](T a, T b){return (a + b);}
#define minf(T) [](T a, T b){return min(a, b);}
#define maxf(T) [](T a, T b){return max(a, b);}
#define rep(i,n) for(int i = 0; i < (n); ++i)
#define all(x) (x).begin(), (x).end()
#define INF 2000000000 //2e9
#define LLINF 2000000000000000000ll //2e18
#define sq(x) ((x)*(x))

using namespace std;
using ll = long long;



//普通のライブラリ

//最大公約数
ll gcd(ll x, ll y) { return y ? gcd(y,x%y) : x;}



/*
// --- MOD用start ---
template <int MOD>
struct ModInt {
  int val;

  ModInt (ll v = 0) : val(v % MOD){
    if(val < 0) val += MOD;
  };

  //unary operator
  ModInt operator+() const {return ModInt(val);}
  ModInt operator-() const {return ModInt(MOD-val);}
  //invはfermatの小定理を用いて
  ModInt inv() const { return this -> pow(MOD - 2);}

  //arithmetic operator
  //thisは自分自身を表し、それに対して+xのような演算がされるとして考えると良い。
  ModInt operator+(const ModInt& x) const { return ModInt(*this) += x;}
  ModInt operator-(const ModInt& x) const { return ModInt(*this) -= x;}
  ModInt operator*(const ModInt& x) const { return ModInt(*this) *= x;}
  ModInt operator/(const ModInt& x) const { return ModInt(*this) /= x;}

  //省略表現について
  ModInt& operator+=(const ModInt& x){
    if((val += x.val) >= MOD) val -= MOD;
    return *this;
  }
  ModInt& operator-=(const ModInt& x){
    if((val -= x.val) <0) val += MOD;
    return *this;
  }
  ModInt& operator*=(const ModInt& x){
    val = (ll(val) * x.val) % MOD;
    return *this;
  }
  ModInt& operator/=(const ModInt& x){
    return *this *= x.inv();
  }

  // ex) mint(2).pow(n);

  //べき乗の計算のライブラリ化
  ModInt pow(ll n) const{
    //繰り返し二乗法でべき乗の計算を行う
    auto x = ModInt(1);
    auto b = *this;
    while(n>0){
      if(n & 1) x*= b;
      n >>= 1;
      b *= b;
    }
    return x;
  }

  //比較演算
  bool operator==(const ModInt& b) const { return val == b.val; }
  bool operator!=(const ModInt& b) const { return val != b.val; }

  //IO
  friend istream& operator>>(istream& is, ModInt& x) noexcept { return is >> x.val; }
  friend ostream& operator<<(ostream& os, const ModInt& x) noexcept { return os << x.val; }
};
const int MOD = 1e9 + 7; // TODO: 適宜ここを調整
using mint = ModInt<MOD>;
//modの時用のnCrの計算、
struct combination {
  vector<mint> fact, ifact;
  combination(int n):fact(n+1),ifact(n+1) {
    assert(n < MOD);
    fact[0] = 1;
    for (int i = 1; i <= n; ++i) fact[i] = fact[i-1]*i;
    ifact[n] = fact[n].inv();
    for (int i = n; i >= 1; --i) ifact[i-1] = ifact[i]*i;
  }
  mint operator()(int n, int k) {
    if (k < 0 || k > n) return 0;
    return fact[n]*ifact[k]*ifact[n-k];
  }
};
//modの時用のnCrの計算、
mint nCr(int n, int r){
  assert(!(n<r));
  assert(!(n<0 || r<0));
  mint ans = 1;
  for(int i = 1; i <= r; i++){
    ans *= n-i+1;
    ans /= i;
  }
  return ans;
}
// --- MOD用end ---*/

/*
// --- Union Find tree start ---
struct UnionFind {
    // union tree
    // used for tracking a set of elements partitioned into a number of disjoint subsets
    // rankなしでの実装

    vector<int> par;
    // par[i] means i's parent
    // (ex) par[3] = 2 means 2 is a parent of 3

    UnionFind(int N) : par(N) {
      //used for initialization
      for(int i = 0; i < N; i++) par[i] = i;
    }

    int root(int x) {
      // データxが属する木の根を再帰で得る
      // root(x) = {xの木の根}
      if (par[x] == x) return x;
      return par[x] = root(par[x]);
    }

    void unite(int x, int y) {
      // xとyのtreeのmerge
      int rx = root(x); //xの根をrx
      int ry = root(y); //yの根をry
      if (rx == ry) return; //同じ時はそのままに
      par[rx] = ry; //xの根rxをyの根ryにつける
    }

    bool same(int x, int y) {
      // 2つのデータx, yが属する木が同じならtrueを返す
      int rx = root(x);
      int ry = root(y);
      return rx == ry;
    }
};
// --- Union Find tree end ---
*/

/*
// --- segment tree start ---
template<typename T>
struct SegTree {
  //segment tree
  //0-indexed

  //木は配列で表現
  vector<T> seg;
  //木の1/2の大きさ
  int size;
  //単位元
  const T unit;
  //比較関数の型
  using F = function<T(T,T)>;
  //mergeする関数
  const F f;

  //単位元がaで大きさnのsegtreeの初期化を行う
  //大きさnとはすなわち葉の数がnであることをいう
  SegTree(int n, const T a, const F f) : unit(a), f(f) {
    size = 1;
    while(size < n) size <<= 1;
    seg.assign(2*size, unit);
    size--;
  }

  //k番目にtを代入する
  void set(int k, const T& t) {
    seg[k + size] = t;
  }

  //関数fによって木を構築する
  void build(){
    for (int i = size - 1; i >=0; i--){
      seg[i] = f(seg[i*2+1],seg[i*2+2]);
    }
  }

  //i番目の要素を返す
  T operator[](const int i){
    return seg[i + size];
  }

  //k番目の値をaに更新する
  void update(int k, T a){
    k = k + size;
    //必要ならばここを変えても良い
    seg[k] = a;
    //自分の親の情報を再帰的に更新していく
    while(k>0){
      k = ((k-1) >> 1);
      seg[k] = f(seg[k*2+1],seg[k*2+2]);
    }
  }

  //[a,b)についてfをした結果を返す
  //例えば、fがminだったらどうなるかを考えるとわかりやすい
  //unitについては足し算なら0,掛け算なら1,minならinf,どんな値が来てもその値を返すもの
  T query(int a, int b, int k, int l, int r) {
    if (r <= a || b <= l) {
      return unit;
    } else if (a <= l && r <= b) {
      return seg[k];
    } else {
      T lch = query(a, b, k * 2 + 1, l, (l + r) / 2);
      T rch = query(a, b, k * 2 + 2, (l + r) / 2, r);
      return f(lch, rch);
    }
  }

  //[a,b)についてfをした結果を返す
  T query(int a, int b) { return query(a, b, 0, 0, size+1); }
};
// --- segment tree end ---
*/


/*
// --- 素数系のライブラリ ---

//素因数分解をpairで出す関数
vector<pair<ll,int>> factorize(ll n) {
  vector<pair<ll,int>> res;
  for (ll i = 2; i*i <= n; ++i) {
    if (n%i) continue;
    res.emplace_back(i,0);
    while (n%i == 0) {
      n /= i;
      res.back().second++;
    }
  }
  if (n != 1) res.emplace_back(n,1);
  return res;
};

//エラトステネスの篩(連続する整数について素数判定や素因数分解をしたい時)
struct Sieve {
  int n;
  vector<int> f, primes;
  Sieve(int n=1):n(n), f(n+1) {
    f[0] = f[1] = -1;
    for (ll i = 2; i <= n; ++i) {
      if (f[i]) continue;
      //素数になった時は...?primesに素数を列挙しておく
      primes.push_back(i);
      f[i] = i;
      for (ll j = i*i; j <= n; j += i) {
        if (!f[j]) f[j] = i;
      }
    }
  }

  //素数か否かを判定する関数
  bool isPrime(int x) { return f[x] == x;}

  //素因数分解をする関数
  vector<int> factorList(int x) {
    vector<int> res;
    while (x != 1) {
      res.push_back(f[x]);
      x /= f[x];
    }
    return res;
  }
  vector<pair<int,int>> factor(int x) {
    vector<int> fl = factorList(x);
    if (fl.size() == 0) return {};
    vector<pair<int,int>> res(1, pair<int,int>(fl[0], 0));
    for (int p : fl) {
      if (res.back().first == p) {
        res.back().second++;
      } else {
        res.emplace_back(p, 1);
      }
    }
    return res;
  }
};
// --- 素数系のライブラリend ---
*/


//when you want to cout double ...
//cout << fixed << setprecision(14) << (double)(i*i*i) << endl;

/* --- ここからコード --- */



int main() {
  int n,a,b,c;
  cin >> n >> a >> b >> c;
  vector<int> l(n);
  rep(i,n){
    cin >> l[i];
  }
  int ans = INF;
  rep(i,1<<(2*n)){
    int tempa = 0;
    int tempb = 0;
    int tempc = 0;
    int tempmagic = 0;
    rep(j,n){
      int digit1 = (i>>(2*j) & 1);
      int digit2 = (i>>(2*j+1) & 1);
      //cout << i << " " <<  j << " " << digit1 << " " << digit2 << endl;
      if(digit1 == 1 && digit2 == 1){
        //Aにいれる
        tempa += l[j];
        tempmagic += 10;
      }
      else if(digit1 == 1 && digit2 == 0){
        tempb += l[j];
        tempmagic += 10;
      }
      else if(digit1 == 0 && digit2 == 1){
        tempc += l[j];
        tempmagic += 10;
      }
      else{
        continue;
      }
    }
    // if(tempa == 98){
    // cout << tempa <<" " <<   tempb << " " << tempc << endl;
    // }
    if(tempa == 0 || tempb == 0 || tempc == 0) continue;

    tempmagic += abs(tempa-a);
    tempmagic += abs(tempb-b);
    tempmagic += abs(tempc-c);

    // if(tempmagic == 140) cout << tempa << " " << tempb << " " << tempc << endl;

    ans = min(ans,tempmagic);
  }
  cout << ans-30 << endl;
  return 0;
}
