//#include "BigInt.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <random> 
#include <string.h>
#include <assert.h>
#include <time.h>

#include <mpir.h>

using namespace std;

//vector<int> small_primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101 };
//
//BigInt randWithRange(int n_bit, BigInt range) { // within [0, range - 1]
//    _ASSERT(range > 1);
//    if (n_bit < 0) {
//        n_bit = log(2, range - 1) + 1;
//    }
//    BigInt CPOW15 = pow(Integer(2), 15);
//    BigInt s;
//    do {
//        s = Integer(0);
//        while (n_bit >= 15) {
//            s *= CPOW15;
//            s += Integer(rand()) % CPOW15;
//            n_bit -= 15;
//        }
//        if (n_bit) {
//            s *= pow(2, n_bit);
//            s += Integer(rand()) % pow(2, n_bit);
//        }
//    } while (s >= range);
//    return s;
//}
//
//BigInt int_pow_mod(BigInt v, BigInt p, BigInt mod) {   // v^p % m 找餘數
//    BigInt res = Integer(1);
//    //cout << "v: " << v << "p: " << p << "mod: " << mod;
//    while (p > 0) {
//        if (p % 2 == 1) {
//            res *= v;
//            res %= mod;
//        }
//        p /= 2;
//        v = v * v % mod;
//    }
//    return res;
//}
//
//bool is_prime(BigInt x) { // 根據費馬小定理判定質數 a^(x-1) = 1 mod x [if x is prime]
//    if (x <= Integer(101)) {
//        for (int p : small_primes)
//            if (x == Integer(p)) return true;
//        return false;
//    }
//    for (int i = 0; i < 20; ++i) {
//        int n_bit = log(2, x);
//        BigInt rnd = randWithRange(n_bit, x - 2) + 2;
//        if (int_pow_mod(rnd, x - 1, x) != Integer(1)) {
//            return false;
//        }
//    }
//    return true;
//}
//
//BigInt find_inverse(Ext a, Ext b) {
//    Ext b0 = b, t, q;
//    Ext x0 = Integer(0), x1 = Integer(1);
//    if (b == Integer(1)) return Integer(1);
//    while (a > Integer(1)) {
//        q = a / b;
//        t = b, b = a % b, a = t;
//        t = x0, x0 = x1 - q * x0, x1 = t;
//    }
//    if (x1 < Integer(0)) x1 += b0;
//    x1 = x1 % b;
//    if (x1 < 0)
//        x1 = x1 + b;
//    return x1.n;
//}

class RSA {
private:
public:
    int bit_size;
    mpz_t p, q;
    mpz_t N, e, d;

    RSA(int, gmp_randstate_t);
    RSA(mpz_t, mpz_t);
    ~RSA();

    //void gen_ed();
    //bool test(const BigInt&, const BigInt&);

    void find_two_prime(int, gmp_randstate_t);
    //bool find_Factorization(BigInt &, BigInt &);
};

RSA::RSA(int bit_size, gmp_randstate_t rstate) {
    this->bit_size = bit_size;
    mpz_init(p);
    mpz_init(q);
    mpz_init(N);
    mpz_init(e);
    mpz_init(d);

    find_two_prime(bit_size, rstate);
    //this->N = p * q;
    //gen_ed();
}

RSA::~RSA() {
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
    mpz_clear(e);
    mpz_clear(d);
}

int main() {
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, time(0));

    RSA rsa(32, rstate);

    return 0;
}

//int main() {
//    mpz_t a;
//    gmp_randstate_t rstate;
//
//    mpz_init(a);
//    gmp_randinit_default(rstate);
//    gmp_randseed_ui(rstate, time(0));
//
//    mpz_urandomb(a, rstate, 100);
//    printf("100 random bits gives: "); mpz_out_str(stdout, 10, a);
//    cout << "\nA: " << mpz_sizeinbase(a, 2) << endl;
//    putchar('\n');
//
//    return 0;
//}

//RSA::RSA(mpz_t p, mpz_t q) {
//    _ASSERT(p != q);
//    this->N = p * q;
//    if (p > q) {
//        this->q = p;
//        this->p = q;
//    }
//    else {
//        this->p = p;
//        this->q = q;
//    }
//    this->bit_size = log(2, this->p - 1);
//}
//
//void RSA::gen_ed() {
//    BigInt pminus = p - 1, qminus = q - 1;
//    BigInt M = pminus * qminus / gcd(pminus, qminus);
//    //cout << p << q;
//    while (true) {
//        e = randWithRange(-1, M - 2) + 2;
//        if (gcd(e, pminus) == 1 && gcd(e, qminus) == 1)
//            break;
//    }
//    d = find_inverse(Ext(e), Ext(M));
//    //cout << e << "d: " << d << "M: " << M << e * d % M;
//    _ASSERT(e * d  % M == Integer(1));
//}
//
//bool RSA::test(const BigInt& rp, const BigInt& rq) {
//    return p == rp && q == rq;
//}

void RSA::find_two_prime(int bit_size, gmp_randstate_t rstate) {
    mpz_t tmp;
    do {
        mpz_urandomb(p, rstate, bit_size);
    } while (!mpz_probab_prime_p(p, 25));

    mpz_init(tmp);
    mpz_ui_pow_ui(tmp, 2, mpz_sizeinbase(p, 2) - 1);
    do {
        mpz_urandomm(q, rstate, tmp);
        mpz_add(q, q, tmp);
    } while (!mpz_probab_prime_p(q, 25));

    _ASSERT(mpz_cmp(q, p) != 0);
    if (mpz_cmp(q, p) < 0) {
        mpz_set(tmp, q);
        mpz_set(q, p);
        mpz_set(p, tmp);
    }
    mpz_clear(tmp);
    printf("p: "); mpz_out_str(stdout, 10, p);
    cout << "\n[ p ]bit_size: " << mpz_sizeinbase(p, 2) << endl;
    printf("q: "); mpz_out_str(stdout, 10, q);
    cout << "\n[ q ]bit_size: " << mpz_sizeinbase(q, 2) << endl;
    _ASSERT(mpz_cmp(p, q) < 0);
}

//bool RSA::find_Factorization(BigInt& rp, BigInt& rq) {
//    BigInt start = Integer(2);
//    BigInt end = sqrt(N);
//    for (BigInt i = start; i <= end; i+=1) {
//        if (is_prime(i)) {
//            if (Integer(0) == (N % i)) {
//                rp = i;
//                rq = N / i;
//                return true;
//            }
//        }
//    }
//    return false;
//}
//
//ostream& operator << (ostream& out, const Ext& obj) {
//    if (obj.neg)
//        printf("-");
//    Print(obj.n);
//    return out;
//}

//int main() {
//
//    gmp_randstate_t s;
//    gmp_randinit_mt(s);
//
//    clock_t start, end;
//    double cpu_time_used;
//
//    int bit_size, samples = 1;
//    for (bit_size = 8; bit_size < 1024; bit_size *= 2) {
//        start = clock();
//
//        //mpz_t n;
//        //mpz_init(n);
//        //mpz_set_ui(n, 2<<20);
//        //mpz_out_str(stdout, 10, n);
//        //cout << "A: " << mpz_sizeinbase(n, 2) << endl;
//        //mpz_clear(n);
//        // 
//        //for (int i = 0; i < samples; ++i)
//        //{
//        //    RSA rsa = RSA(bit_size);
//        //    //cout << "p:" << rsa.p << "q:" << rsa.q << "N:" << rsa.N << endl;
//        //    //cout << "e:" << rsa.e << "d:" << rsa.d << endl;
//
//        //    BigInt rp, rq;
//        //    _ASSERT(rsa.find_Factorization(rp, rq));
//        //    //cout << "p:" << rsa.p << "q:" << rsa.q;
//        //    _ASSERT(rsa.test(rp, rq));
//        //}
//        end = clock();
//        cpu_time_used = ((double)(end - start) / samples) / CLOCKS_PER_SEC;
//        printf("Avg-Time = %f with bit-size %i \n", cpu_time_used, bit_size);
//    }
//    return 0;
//}