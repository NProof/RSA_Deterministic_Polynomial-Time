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

class RSA {
private:
public:
    int bit_size;
    mpz_t p, q;
    mpz_t N, e, d;

    RSA(int, gmp_randstate_t);
    RSA(mpz_t, mpz_t);
    ~RSA();

    void gen_ed(gmp_randstate_t);
    bool test(const mpz_t, const mpz_t);

    void find_two_prime(int, gmp_randstate_t);
    bool find_Factorization(mpz_t, mpz_t);
    bool factorizationCase1(mpz_t, mpz_t, mpz_t);
};

RSA::RSA(int bit_size, gmp_randstate_t rstate) {
    this->bit_size = bit_size;
    mpz_init(p);
    mpz_init(q);
    mpz_init(N);
    mpz_init(e);
    mpz_init(d);

    find_two_prime(bit_size, rstate);
    mpz_mul(N, p, q);
    //cout << "N: "; mpz_out_str(stdout, 10, N);
    //cout << "\n[ N ]bit_size: " << mpz_sizeinbase(N, 2) << endl;
    gen_ed(rstate);
    //cout << "\ne: "; mpz_out_str(stdout, 10, e);
    //cout << "\nd: "; mpz_out_str(stdout, 10, d);
}

RSA::~RSA() {
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
    mpz_clear(e);
    mpz_clear(d);
}

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

void RSA::gen_ed(gmp_randstate_t rstate) {
    mpz_t pminus , qminus, phi, gcd;
    mpz_init(pminus); mpz_init(qminus);
    mpz_init(phi); mpz_init(gcd);
    mpz_sub_ui(pminus, p, 1);
    mpz_sub_ui(qminus, q, 1);
    mpz_mul(phi, pminus, qminus);
    while (true) {
        mpz_urandomm(e, rstate, phi);
        mpz_gcd(gcd, e, phi);
        if (mpz_cmp_ui(gcd, 1) == 0)
            break;
    }
    mpz_invert(d, e, phi);
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul(tmp, e, d);
    mpz_mod(tmp, tmp, phi);
    _ASSERT(mpz_cmp_ui(tmp, 1) == 0);
    mpz_clear(tmp);
    mpz_clear(gcd);
    mpz_clear(phi);
    mpz_clear(qminus);
    mpz_clear(pminus);
}

bool RSA::test(const mpz_t rp, const mpz_t rq) {
    //cout << "p:";
    //mpz_out_str(stdout, 10, p); putchar('\n');
    //cout << "rp:";
    //mpz_out_str(stdout, 10, rp); putchar('\n');
    return (mpz_cmp(p, rp) == 0) && (mpz_cmp(q, rq) == 0);
}

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
    //printf("p: "); mpz_out_str(stdout, 10, p);
    //cout << "\n[ p ]bit_size: " << mpz_sizeinbase(p, 2) << endl;
    //printf("q: "); mpz_out_str(stdout, 10, q);
    //cout << "\n[ q ]bit_size: " << mpz_sizeinbase(q, 2) << endl;
    _ASSERT(mpz_cmp(p, q) < 0);
}

bool RSA::find_Factorization(mpz_t rp, mpz_t rq) {
    mpz_t ed, n3f2;
    mpz_init(ed);
    mpz_init(n3f2);
    mpz_mul(ed, e, d);
    mpz_pow_ui(n3f2, N, 3);
    mpz_sqrt(n3f2, n3f2);
    if (mpz_cmp(ed, n3f2) <= 0) {
        mpz_clear(n3f2);
        cout << "factorizationCase1" << endl;
        return factorizationCase1(ed, rp, rq);
        mpz_clear(ed);
    }
    cout << "other" << endl;
    mpz_t i, end, r;
    mpz_init(i);
    mpz_init(end);
    mpz_init(r);
    mpz_set_si(i, 2);
    mpz_sqrt(end, N);
    while (mpz_cmp(i, end) <= 0) {
        if (mpz_probab_prime_p(i, 25)) {
            mpz_mod(r, N, i);
            if (mpz_cmp_ui(r, 0) == 0) {
                mpz_set(rp, i);
                mpz_divexact(rq, N, i);
                mpz_clear(i);
                mpz_clear(end);
                mpz_clear(r);
                return true;
            }
        }
        mpz_add_ui(i, i, 1);
    }
    mpz_clear(i);
    mpz_clear(end);
    mpz_clear(r);
    //mpz_clear(n3f2);
    //mpz_clear(ed);
    return false;
}


bool RSA::factorizationCase1(mpz_t ed, mpz_t rp, mpz_t rq) {
    mpz_t k, kp, edtmp;
    mpz_init(k); mpz_init(kp); mpz_init(edtmp);
    mpz_sub_ui(kp, kp, 1);
    mpz_cdiv_q(kp, kp, N);
    mpz_ui_sub(edtmp, 1, ed);
    bool founded = false;
    mpz_t paq, start, end;
    mpz_init(paq); mpz_init(start); mpz_init(end);
    mpz_t mid, midN;
    mpz_init(mid); mpz_init(midN);
    int icase;
    for (int i = 1; i < 6; ++i) {
        mpz_add_ui(k, kp, i);
        mpz_cdiv_q(paq, edtmp, k);
        mpz_add_ui(paq, paq, 1);
        mpz_add(paq, paq, N);
        mpz_set_ui(start, 2);
        mpz_sub_ui(end, paq, 1);
        mpz_fdiv_q_ui(end, end, 2);
        while (mpz_cmp(start, end) <= 0) {
            mpz_sub(mid, end, start);
            mpz_tdiv_q_ui(mid, mid, 2);
            mpz_add(mid, mid, start);
            mpz_sub(midN, paq, mid);
            mpz_mul(midN, midN, mid);
            icase = mpz_cmp(midN, N);
            if (icase < 0)
                mpz_add_ui(start, mid, 1);
            else if (icase > 0)
                mpz_sub_ui(end, mid, 1);
            else {
                founded = true;
                mpz_set(rp, mid);
                mpz_divexact(rq, N, mid);
                break;
            }
        }
        if (founded)
            break;
    }
    mpz_clear(mid); mpz_clear(midN);
    mpz_clear(paq); mpz_clear(start); mpz_clear(end);
    mpz_clear(kp);
    mpz_clear(k);
    return founded;
}

int main() {
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, time(0));

    clock_t start, end;
    double cpu_time_used;

    int bit_size = 21, samples = 100;
    for (int i = 0; i < samples; ++i)
    {
        start = clock();
        RSA rsa = RSA(bit_size, rstate);
        //cout << "p:" << rsa.p << "q:" << rsa.q << "N:" << rsa.N << endl;
        //cout << "e:" << rsa.e << "d:" << rsa.d << endl;

        mpz_t rp, rq;
        mpz_init(rp);
        mpz_init(rq);
        bool ans1 = rsa.find_Factorization(rp, rq);
        _ASSERT(ans1);
        //cout << "p:";
        //mpz_out_str(stdout, 10, rsa.p); putchar('\n');
        //cout << "rp:";
        //mpz_out_str(stdout, 10, rp); putchar('\n');
        bool ans2 = rsa.test(rp, rq);
        _ASSERT(ans2);
        mpz_clear(rp);
        mpz_clear(rq);
        end = clock();
        cpu_time_used = ((double)(end - start) / samples) / CLOCKS_PER_SEC;
        printf("time = %f with bit-size %i \n", cpu_time_used, bit_size);
    }
    return 0;
}
//
//int main() {
//    gmp_randstate_t rstate;
//    gmp_randinit_default(rstate);
//    gmp_randseed_ui(rstate, time(0));
//
//    clock_t start, end;
//    double cpu_time_used;
//
//    int bit_size, samples = 100;
//    for (bit_size = 8; bit_size <= 256; bit_size += 8) {
//        start = clock();
//        for (int i = 0; i < samples; ++i)
//        {
//            RSA rsa = RSA(bit_size, rstate);
//            //cout << "p:" << rsa.p << "q:" << rsa.q << "N:" << rsa.N << endl;
//            //cout << "e:" << rsa.e << "d:" << rsa.d << endl;
//
//            mpz_t rp, rq;
//            mpz_init(rp);
//            mpz_init(rq);
//            _ASSERT(rsa.find_Factorization(rp, rq));
//            //cout << "p:";
//            //mpz_out_str(stdout, 10, rsa.p); putchar('\n');
//            //cout << "q:";
//            //mpz_out_str(stdout, 10, rsa.q); putchar('\n');
//            _ASSERT(rsa.test(rp, rq));
//            mpz_clear(rp);
//            mpz_clear(rq);
//        }
//        end = clock();
//        cpu_time_used = ((double)(end - start) / samples) / CLOCKS_PER_SEC;
//        printf("Avg-Time = %f with bit-size %i \n", cpu_time_used, bit_size);
//    }
//    return 0;
//}