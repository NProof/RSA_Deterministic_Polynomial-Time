#include <bits/stdc++.h>
#include <chrono>
#include <random>

std::vector< int > small_prime = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                              47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101 };

unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 g1(std::chrono::system_clock::now().time_since_epoch().count());

int int_pow_mod( int v, int p, int mod ){
    int res = 1;
    while( p ){
        if( p & 1 ) res = 1LL * res * v % mod;
        p >>= 1;
        v = 1LL * v * v % mod;
    }
    return res;
}

bool is_prime( int x ){
    if( x <= 101 ){
        for( int p : small_prime )
            if( x == p ) return true;
        return false;
    }
    for( int i = 0; i < 20; ++i ){
        int rnd = rand() % ( x - 2 ) + 2; // [ 0, x - 3 ] -> [ 2, x - 1 ]
        if( int_pow_mod( rnd, x - 1, x ) != 1 )
            return false;
    }
    return true;
}

void find_two_prime(uint32_t &x, uint32_t &y){
    uint32_t m = pow(std::numeric_limits<uint32_t>::max(), (double)10/21);
    uint32_t M = pow(std::numeric_limits<uint32_t>::max(), (double)12/21);
//    std::cout << std::hex << " m " << m << std::endl;
    while(true){
        x = llrint( (double)g1() / g1.max() * ( M - m ) + m );
        if(is_prime(x)){
            break;
        }
    }
    uint32_t m2 = pow(std::numeric_limits<uint32_t>::max(), (double)9/21);
    uint32_t M2 = (std::numeric_limits<uint32_t>::max() / x);
    while(true){
        y = llrint( (double)g1() / g1.max() * ( M2 - m2 ) + m2 );
        if(is_prime(y)){
            break;
        }
    }
}

uint32_t gcd(uint32_t &x, uint32_t &y){
    uint32_t tmp = x % y;
    if(tmp == 0)
        return y;
    else
        return gcd(y, tmp);
}

long mul_inv(long a, long b)
{
	long long b0 = b, t, q;
	long long x0 = 0, x1 = 1;
	if (b == 1) return 1;
	while (a > 1) {
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

void gen_ed(uint32_t &p, uint32_t &q, uint32_t &e, long &d, long &M){
    uint32_t pminus = p-1, qminus = q-1;
    std::cout << std::dec << pminus << " , " << qminus;
    M = 1LL * pminus*qminus / gcd(pminus, qminus);
    std::cout << " " << M;
    while(true){
        e = llrint( (double)g1() / g1.max() * ( M - 2 ) + 2 );
        if(gcd(e, pminus) == 1 && gcd(e, qminus) == 1)
            break;
    }
    std::cout << " " << e;
    d = mul_inv(e, M);
    std::cout << " " << d << " ";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie( 0 ), std::cout.tie( 0 );

    uint32_t x = 0, y = 0, e = 0;
    long d = 0, M = 0;
    for(int i = 0; i < 10000; ++i){
        std::cout << " *";
        find_two_prime(x, y);
        std::cout << " " << x << " " << y << " ";
        gen_ed(y, x, e, d, M); // 1LL* e*d % M  = 1
        std::cout << "< ";
        std::cout << std::dec << 1LL* e*d % M;
        std::cout << " > " << std::endl;
//        std::cout << std::hex << e << " & " << d << std::endl;
//        std::cout << std::hex << x * y << " - " << x << ", " << y << std::endl;
    }
    return 0;
}

// * 101723 28607 28606 , 101722 1454929766 813851573 1233904035 < 1 >
//* 130409 23159 23158 , 130408 1509994232 353883097 1341512809 < 1 >

// http://h0rnet.hatenablog.com/entry/2016/07/27/%E9%9A%A8%E6%A9%9F%E7%AE%97%E6%B3%95%E5%BF%AB%E9%80%9F%E5%88%A4%E6%96%B7%E8%B3%AA%E6%95%B8_O%28_lg_V_%29
// https://www.hexadecimaldictionary.com/hexadecimal/0xfcb5e9/
// https://www.mathsisfun.com/prime_numbers.html
// https://www.inf.pucrs.br/~calazans/graduate/TPVLSI_I/RSA-oaep_spec.pdf
// https://cp-algorithms.com/algebra/chinese-remainder-theorem.html
// https://rosettacode.org/wiki/Modular_inverse#Recursive_implementation