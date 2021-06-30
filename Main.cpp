#include <bits/stdc++.h>

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
    std::cout << " " << e;
    d = mul_inv(e, M);
    std::cout << " " << d << " ";
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie( 0 ), std::cout.tie( 0 );

    uint32_t x = 130409, y = 23159, e = 353883097;
    long d = 0, M = 0;
        std::cout << " *";
        std::cout << " " << x << " " << y << " ";
        gen_ed(y, x, e, d, M); // 1LL* e*d % M  = 1
        std::cout << "< ";
        std::cout << std::dec << 1LL* e*d % M;
        std::cout << " > " << std::endl;
//        std::cout << std::hex << e << " & " << d << std::endl;
//        std::cout << std::hex << x * y << " - " << x << ", " << y << std::endl;

    return 0;
}

//* 130409 23159 23158 , 130408 1509994232 353883097 1341512809 < 1 >
