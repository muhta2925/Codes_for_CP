i128 norm(i128 a,i128 m){
    a %= m;
    if(a<0)a+=m;
    return a;
}
i128 egcd(i128 a,i128 b,i128 &x,i128 &y){
    if(b==0){
        x = 1;
        y = 0;
        return a;
    }
    i128 x1,y1;
    i128 g = egcd(b,a%b,x1,y1);
    x = y1;
    y = x1 - (a/b)*y1;
    return g;
}
bool modinv(ll a,ll m,ll &inv){
    if(m<=1)return false;
    i128 x,y;
    ll g = (ll) egcd(norm(a,m),m,x,y);
    if(g != 1)return false;
    inv = (ll) norm(x,m);
    return true;
}
bool solve_linear(ll a,ll b,ll m,ll &x0,ll &g,ll &M){
    i128 x,y;
    g = egcd(a,m,x,y);
    if(b%g!=0)return false;
    ll A = a/g;
    ll B = b/g;
    M = m/g;
    ll inv;
    if(!modinv(A,M,inv)){
        return false;
    }
 
    x0 = norm(B*inv,M);
    // x = xo + (0,1,2,3...g-1*M);
    return true;
}


bool classical_crt(vector<ll>&a,vector<ll>&m,ll &x0,ll &M){
    int n = a.size();
    if(n!=(int)m.size())return false;

    // pairwise coprime check: etay lagbe na
    // for(ll i=0;i<n;i++){
    //     for(ll j=0;j<i;j++){
    //         if(gcd(a[i],a[j])!=1)return false;
    //     }
    // }

    i128 N = 1;
    for(ll i=0;i<(int)m.size();i++){
        if(N > (i128)LLONG_MAX/m[i]){
            return false;
        }
        N *= m[i];
    }
    i128 ans = 0;
    for(ll i=0;i<n;i++){
        if(m[i]==1)continue;
        i128 x,y;
        i128 Mi = N/m[i];
        ll g = egcd(norm(Mi,m[i]),m[i],x,y);
        if(g!=1)return false;
        ll yi = (ll) norm(x,m[i]);
        ll ai = (ll)norm(a[i],m[i]);
        ans = (ans + (ai*yi*Mi)) % N;
    }
    x0 = (ll) ans;
    M = (ll)N;
    return true;
}
bool merge_crt(ll a,ll b,ll m,ll n,ll &x0,ll &M){
    if(m<=0 || n<=0)return false;
    a = norm(a,m);
    b = norm(b,n);
    ll g = gcd(m,n);
    i128 dif = b-a;
    if(dif%g!=0)return false;
    i128 m1 = m/g;
    i128 n1 = n/g;
    i128 d = dif/g;
    i128 k = 0;
    if(n1>1){
        i128 x,y;
        ll m1inv;
        modinv(m1,n1,m1inv);
        m1inv = norm(m1inv,n1);
        k = norm(m1inv*d,n1);
    }
    i128 L = (i128)m * n1;
    if(L > (i128)LLONG_MAX)return false;
    i128 x = norm(a+((i128)m*k),L);
    x0 = (ll)x;
    M = (ll)L;
    return true;
}
bool generalized_crt(vector<ll>&a,vector<ll>&m,ll &x0,ll &M){

    x0 = 0LL;
    M = 1LL;

    if(a.size() != m.size()){
        return false;
    }

    ll n = a.size();
    for(ll i=0; i<n; i++){
        
        if(!merge_crt(a[i],x0,m[i],M,x0,M)){
            return false;
        }
    }

    return true;
}
void solve(int tc) {
   
    vector<ll> a = {7, -2, 4};
    vector<ll> m = {1, 5, 6};
    
    ll x0,M;

    if(generalized_crt(a,m,x0,M)){
        cout << x0 << " " << M << "\n";
    }
    else {
        cout << "Not found\n";
    }
}
