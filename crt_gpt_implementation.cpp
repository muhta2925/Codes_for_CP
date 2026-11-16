
#include <bits/stdc++.h>
using namespace std;

using ll = long long;
using i128 = __int128_t;

// Normalize into [0, m-1].
i128 norm(i128 a, i128 m) {
    a %= m;
    if (a < 0) a += m;
    return a;
}

// Extended Euclidean Algorithm.
i128 extgcd(i128 a, i128 b, i128 &x, i128 &y) {

    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }

    i128 x1, y1;

    i128 g = extgcd(b, a % b, x1, y1);

    x = y1;
    y = x1 - (a / b) * y1;

    return g;
}

enum class MergeStatus {
    OK,
    NO_SOLUTION,
    OVERFLOW,
    INVALID_INPUT
};

// Merge:
// x = a (mod m)
// x = b (mod n)
//
// Success:
// a = new normalized solution
// m = lcm(m,n)
MergeStatus merge_crt(ll &a, ll &m, ll b, ll n) {

    if (m <= 0 || n <= 0) {
        return MergeStatus::INVALID_INPUT;
    }

    i128 A = norm(a, m);
    i128 B = norm(b, n);

    ll g = std::gcd(m, n);

    i128 diff = B - A;

    // Compatibility condition.
    if (diff % g != 0) {
        return MergeStatus::NO_SOLUTION;
    }

    i128 m1 = (i128)m / g;
    i128 n1 = (i128)n / g;
    i128 d = diff / g;

    i128 k = 0;

    if (n1 > 1) {

        i128 inv, y;

        extgcd(m1, n1, inv, y);

        inv = norm(inv, n1);

        k = norm(d * inv, n1);
    }

    // Calculate LCM safely.
    i128 new_mod = (i128)m * n1;

    if (new_mod > LLONG_MAX) {
        return MergeStatus::OVERFLOW;
    }

    i128 new_a = norm(A + (i128)m * k, new_mod);

    // Update only after successful calculation.
    a = (ll)new_a;
    m = (ll)new_mod;

    return MergeStatus::OK;
}

// Iterative Generalized CRT.
MergeStatus generalized_crt(
    const vector<ll> &r,
    const vector<ll> &mod,
    ll &x,
    ll &M
) {

    if (r.size() != mod.size()) {
        return MergeStatus::INVALID_INPUT;
    }

    // Initially all integers are valid.
    x = 0;
    M = 1;

    for (int i = 0; i < (int)r.size(); i++) {

        MergeStatus status =
            merge_crt(x, M, r[i], mod[i]);

        if (status != MergeStatus::OK) {
            return status;
        }
    }

    return MergeStatus::OK;
}

int main() {

    vector<ll> r = {2, 5, 3};
    vector<ll> mod = {6, 9, 7};

    ll x, M;

    MergeStatus status =
        generalized_crt(r, mod, x, M);

    if (status == MergeStatus::OK) {
        cout << x << " " << M << '\n';
    }
    else if (status == MergeStatus::NO_SOLUTION) {
        cout << "No solution\n";
    }
    else if (status == MergeStatus::OVERFLOW) {
        cout << "Overflow\n";
    }
    else {
        cout << "Invalid input\n";
    }

    return 0;
}
