#include <chrono>

#include <gmpxx.h>
#include <vector>
#include <ctime>
#include <cmath>
#include <assert.h>
#include <random>
#include <iostream>
#include <fplll/fplll.h>


std::vector<mpz_class> xi(mpz_class& p,int gam,int eta,int rho,gmp_randstate_t state,int m){

    std::vector<mpz_class>ciphertexts(m,0);
    // std::vector<mpz_class>ciphertexts{mpz_class("56405845507494530020941008480572940286181689237258854"),mpz_class("39904821464460948494700284192336525523357407545067668"),mpz_class("56294991345433284900612805613249060787237279328022519")};

    mpz_class bound_q(mpz_class(1) << (gam-eta));
    mpz_class bound_r = mpz_class(1) << rho;
    mpz_class q;
    mpz_class r;
    for(int i=0;i<m;i++){
        mpz_urandomm(q.get_mpz_t(),state,bound_q.get_mpz_t());
        mpz_urandomm(r.get_mpz_t(),state,bound_r.get_mpz_t());
        ciphertexts[i] = p*q + r;
    }

    return ciphertexts;
}

fplll::ZZ_mat<mpz_t> create_lattice(const std::vector<mpz_class>& ciphertexts) {
    int n = ciphertexts.size();
    fplll::ZZ_mat<mpz_t> L(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                mpz_t one;
                mpz_init_set_ui(one, 1); 
                L[i][j] = one;
            } else {
                L[i][j] = mpz_t{0};
            }
        }
        mpz_t c;
        mpz_init(c);
        mpz_set(c,ciphertexts[i].get_mpz_t());
        L[i][n-1] = c;
    }
    return L;
}

fplll::ZZ_mat<mpz_t> create_lattice_prime(const fplll::MatrixRow<fplll::Z_NR<mpz_t>>& ai, const std::vector<mpz_class>& c, const mpz_class& C) {
    int n = ai.size();
    fplll::ZZ_mat<mpz_t> L(n+1, n+1);
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            if (i == j) {
                mpz_t one;
                mpz_init_set_ui(one, 1); 
                L[i][j] = one;
            } else {
                L[i][j] = mpz_t{0};

            }
        }
        if (i == n) {
            mpz_class sum{0};
            for (int j = 0; j < n; j++) {
                mpz_addmul(sum.get_mpz_t(),  c[j].get_mpz_t(), mpz_class(ai[j].get_si()).get_mpz_t());
            }
            sum *= C;
            mpz_t res;
            mpz_init(res);
            mpz_set(res,sum.get_mpz_t());
            L[i][n] = res;
        } else {
            mpz_class res = C * mpz_class(ai[i].get_si());
            mpz_t x;
            mpz_init(x);
            mpz_set(x,res.get_mpz_t());
            L[i][n] = x;
        }
    }
    return L;
}

mpz_class attack(const std::vector<mpz_class>& c, const int m,int iter){
    fplll::ZZ_mat<mpz_t> L = create_lattice(c);
    if(iter == 0){
        auto lll1 = std::chrono::steady_clock::now();

        fplll::lll_reduction(L);

        auto lll2 = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(lll2 - lll1).count();
        std::cout << "LLL1: " << total_time << " ms" << std::endl;
    }else{
        fplll::lll_reduction(L);

    }
   
    int rows = L.get_rows();
    int columns = L.get_cols();
    mpz_class gcd{0};

    for(int i =0; i < rows; i++){
        mpz_class sum{0};
        mpz_t am;
        for(int j = 0; j < m-1;j++){
            mpz_class ai;
            L[i][j].get_mpz(ai.get_mpz_t());
            mpz_addmul(sum.get_mpz_t(),  c[j].get_mpz_t(), ai.get_mpz_t());
        }
        mpz_class last;
        L[i][columns-1].get_mpz(last.get_mpz_t());
        mpz_class tmp{last - sum};
        mpz_init(am);
        mpz_fdiv_q(am,tmp.get_mpz_t(),c[m-1].get_mpz_t());

        L[i][columns-1] = am;
        mpz_class alpha = 0;
        for (int j = 0; j < m-1; j++) {
            mpz_class tmp;
            L[i][j].get_mpz(tmp.get_mpz_t());
            mpz_class log_val = mpz_class(mpz_sizeinbase(tmp.get_mpz_t(), 2));
            if (log_val > alpha) {
                alpha = log_val;
            }
        }
        mpz_class C {mpz_class(2) << 100};
        fplll::ZZ_mat<mpz_t> L_prime = create_lattice_prime(L[i],c,C);
        if(iter == 0){
            auto start_time = std::chrono::steady_clock::now();

            fplll::lll_reduction(L_prime);
            auto end_time = std::chrono::steady_clock::now();
            auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
            std::cout << "LLL2: " << total_time << " ms" << std::endl;
        }else{
            fplll::lll_reduction(L_prime);
        }
       
        int rows_p = L_prime.get_rows();
        int cols_p = L_prime.get_cols();
        for(int j=0;j<rows_p;j++){
            if(L_prime[j][cols_p-1].get_si() == 0){
                mpz_class n = c[0] - (abs(L_prime[j][0].get_si()));
                mpz_class m = c[1] - (abs(L_prime[j][1].get_si()));

                mpz_gcd(gcd.get_mpz_t(),n.get_mpz_t(),m.get_mpz_t());
                return gcd;
            }
        }

    }
    return gcd;
}

int main(){
    //(100,200) m=3 rho 12
    //(100,328) m=5 rho sqrt(eta)
    //(100,117) m=2 rho sqrt(eta)

    //(200,253) m=2 rho sqrt(eta)
    //(200,300) m=3 rho sqrt(eta)
    //(200,443) m=5 rho 12

    //(50,79) m=2 rho srqt(eta)
    //(50,210) m=5 rho sqrt(eta)
    //(50,123) m=3 rho sqrt(eta)
    std::random_device rd;  // obtain a random seed from the hardware
    std::mt19937 eng(rd());  // seed the generator
    std::uniform_int_distribution<> distr(100, 1000000000);  // define the range
    unsigned int gamma_j = 0;
    std::vector<std::tuple<int, int>> vec;

    std::vector<int> gamma {253,300,443};
    std::vector<int> m {2,3,5};

    for(int j =0;j<3;j++){
        int cnt = 0;
        auto start_time = std::chrono::steady_clock::now();

        for(int i =0;i<100;i++){

            gmp_randstate_t state;
            gmp_randinit_mt(state);
            gmp_randseed_ui(state, distr(eng));
            const unsigned int eta = 200;
            const unsigned int rho =(int) sqrt(eta);
            mpz_class prime;
            mpz_urandomb(prime.get_mpz_t(), state, eta);
            mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());

            std::vector<mpz_class> c = xi(prime,gamma[j],eta,rho,state,m[j]);
            mpz_class p = attack(c,m[j],i);

            if(mpz_cmp(prime.get_mpz_t(),p.get_mpz_t()) == 0){
                cnt +=1;
            }

        }
            auto end_time = std::chrono::steady_clock::now();
            auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
            std::cout << "Total time: " << total_time << " ms" << std::endl;
            std::cout << "Success: " << cnt << std::endl;
        }

    
    }
    