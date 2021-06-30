// RSA_E_D.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <random> 
#include <string.h>
#include <bitset>
#include <fstream>


using namespace std;


class RSA {

private:
    uint32_t e = 0;

public:
    uint32_t N = 0;
    vector<int> small_primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                              47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101 };

    RSA(uint32_t N, uint32_t e) {
        this->N = N;
        this->e = e;
    }


    int find_inverse(long a, long b) {
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

    long find_Factorization(uint32_t N) {    // get p and q
        uint16_t start = 2;
        uint32_t end = sqrt(N);
        for (int i = start; i < end; i++) {
            if (is_prime(i)) {
                if (0 == (N % i)) {
                    // cout << "N and i: "<< N << " " << i << endl; // N = p . q, i = prime
                    return i;
                }
            }
        }
        return -1;
    }

    bool is_prime(int x) { // 根據費馬小定理判定質數 a^(x-1) = 1 mod x [if x is prime]
        if (x <= 101) {
            for (int p : this->small_primes)
                if (x == p) return true;
            return false;
        }
        for (int i = 0; i < 20; ++i) {
            int rnd = rand() % (x - 2) + 2; // [ 0, x - 3 ] -> [ 2, x - 1 ]
            if (int_pow_mod(rnd, x - 1, x) != 1)
                return false;
        }
        return true;
    }

    int int_pow_mod(int v, int p, int mod) {   // v^p % m 找餘數
        int res = 1;
        while (p) {
            if (p & 1) res = 1LL * res * v % mod;
            p >>= 1;
            v = 1LL * v * v % mod;
        }
        return res;
    }

	//void find_two_prime(uint32_t &x, uint32_t &y){

	//    uint32_t m = pow(std::numeric_limits<uint32_t>::max(), (double)10/21);
	//    uint32_t m = pow(std::numeric_limits<uint32_t>::max(), (double)12/21);
	////    std::cout << std::hex << " m " << m << std::endl;
	//    while(true){
	//        x = llrint( (double)g1() / g1.max() * ( m - m ) + m );
	//        if(is_prime(x)){
	//            break;
	//        }
	//    }
	//    uint32_t m2 = pow(std::numeric_limits<uint32_t>::max(), (double)9/21);
	//    uint32_t m2 = (std::numeric_limits<uint32_t>::max() / x);
	//    while(true){
	//        y = llrint( (double)g1() / g1.max() * ( m2 - m2 ) + m2 );
	//        if(is_prime(y)){
	//            break;
	//        }
	//    }
	//}

	THIS.
    /*int find_inverse(long a, long b) {
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
    }*/

    //BigInt find_Factorization(BigInt N) {    // get p and q
    //    BigInt start(2);
    //    BigInt end = sqrt(N);
    //    for (BigInt i = start; i < end; i += 1) {
    //        if (is_prime(i)) {
    //            if ((N % i) == 0) {
    //                // cout << "N and i: "<< N << " " << i << endl; // N = p . q, i = prime
    //                return i;
    //            }
    //        }
    //    }
    //    return BigInt(-1);
    //}
	THIS.

    long encryption_RSA(long p) {
        long plaintext = p;
        long cipher = (long long)(pow(plaintext, this->e)) % N;
        return cipher;
    }
     
    long decrption_RSA(long c) {
        // long cipher = c;
        // long plain = pow(c, )
        // return plain
        return 0;
    }

    void plaintext_input(char* text) {
        int length_text = (unsigned)strlen(text);
        // cout << length_text << endl;
        long* binary_text;
        binary_text = new long[500];
        for (int i = 0; i < length_text; i++) {
            // binary_text[i] = (long)bitset<32>(text[i]).to_string();
            // cout << bitset<8>(text[i]).to_string() << endl;  // 轉為binary
        }
        return;
    }

    char *convertDecimalToBinary(int n) {
        //char text[32] = { '\0' };
        char * text;
        text = new char[32] { '\0' };
        long long binaryNumber = 0;

        int remainder, i = 1, step = 1;

        while (n != 0)
        {
            remainder = n % 2;
            // cout << "Step " << step++ << ": " << n << "/2, Remainder = " << remainder << ", Quotient = " << n/2 << endl;
            text[step] = remainder;
            step++;
            n /= 2;
            binaryNumber += (remainder * i);
            i *= 10;
        }
        return text;
    }

    long long convertBinaryToDecimal(long long n) {
        int decimalNumber = 0, i = 0, remainder;
        while (n != 0) {
            remainder = n % 10;
            n /= 10;
            decimalNumber += remainder * pow(2, i);
            ++i;
        }
        return decimalNumber;
    }

    // a 32 bits of message cut with 24bit + 3's zero padding bits and random number 5 bits


    void binaryFileOperate() {
        int bias_a = 5, bias_r = 2; 

        FILE* pFile;
        fopen_s(&pFile, "bad.jpg", "rb");
        if (pFile != NULL) {
            int c;
            while (1) {
                c = fgetc(pFile);
                if (feof(pFile)) {
                    break;
                }
                cout  << hex << c << endl;
            }
        } else {
            cout << "Read file Fail !! \n";
        }
  
        fclose(pFile);

        //char buffer[500] = { '\0' };
        // bad.jpg = 217,862 bytes = 1,742,896 bits
        //ifstream readfile("D:/Dante/Learning/NCHU/Courses/計算機數學/VisualStudio/RSA/textFile.txt", ios::binary);    // textFile.txt || bad.jpg
        //if (readfile.is_open()){
        //    vector<unsigned char> buffer(std::istreambuf_iterator<char>(readfile), {});
        //    
        //    static int inputLength = buffer.size();
        //    cout << inputLength << endl;
        //    cout << (int*)&readfile << endl;
        //    char * tmp = (char*)&readfile;
        //    cout << sizeof(readfile);
        //    for (int i = 0; i < inputLength; i++) {
        //        for (int j = 0; j < 4; j++) {
        //            
        //        }
        //    }
        //}
        //readfile.close();
        //       
        //cout << sizeof(char32_t) << endl; // 4 bytes
        //ofstream writeFile("outputFile.txt", ios::out);

    }

};

int main() {

    // test
    /*int num = 14;
    int* pnum = (int*)num;
    cout << pnum << endl;*/



    RSA rsa_test = RSA(933192884, 645944641); // Test the 
    // cout << bitset<32>(933192884).to_string() << endl;
    uint32_t ans = rsa_test.find_Factorization((uint32_t)101723 * 28607);
    // cout << "ans: " << ans << endl;
    int charLength = 500;
    char* text;
    text = new char[charLength];
    // cin >> text;
    //rsa_test.plaintext_input(text);

    rsa_test.binaryFileOperate();
    cout << "this OK!\n";

    //cout << sizeof(text)/sizeof(text[0]) << endl;
    //char * a = rsa_test.convertDecimalToBinary(1195);
    //long long b = rsa_test.convertBinaryToDecimal(10010101011);
    //for (int i = 32 - 1; i >= 0; i--) { // 測試印出binary
    //    cout << (int)a[i] << " ";
    //}

    //system("Pause");



    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
