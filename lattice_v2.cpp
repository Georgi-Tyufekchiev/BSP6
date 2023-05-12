#include <iostream>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>
#include <fplll.h>
#include <random>

using namespace std;
using namespace fplll;
const mpz_class lambd{mpz_class(2) << 100};
vector<mpz_class> genXi(int rho, int eta, int gam,gmp_randstate_t state) {
    mpz_class x1 = 5;
    mpz_class x0 = 1;
    mpz_class prime;
    mpz_class gcd;
    mpz_class r;
    mpz_class q0;
    mpz_class q1;
    mpz_urandomb(prime.get_mpz_t(), state, eta);
    mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    mpz_gcd(gcd.get_mpz_t(),x1.get_mpz_t(),x0.get_mpz_t());

    while (x1 > x0 || gcd != 1) {
 
        mpz_gcd(gcd.get_mpz_t(),x1.get_mpz_t(),x0.get_mpz_t());
        // Generate r of rho bits
        mpz_urandomb(r.get_mpz_t(),state, rho);
        // Generate q0 of gamma - eta bits
        mpz_urandomb(q0.get_mpz_t(),state, gam - eta);
        // Generate q1 of gamma - eta  - 1bits
        mpz_urandomb(q1.get_mpz_t(),state, gam - eta - 1);
        // Compute x0 = prime * q0
        mpz_mul(x0.get_mpz_t(),prime.get_mpz_t(),q0.get_mpz_t());
        // Compute x1 = prime * q1 + r
        mpz_mul(x1.get_mpz_t(),prime.get_mpz_t(),q1.get_mpz_t());
        mpz_add(x1.get_mpz_t(),x1.get_mpz_t(),r.get_mpz_t());

    }
    vector<mpz_class> pk{x0,x1,prime};
    gmp_randclear(state);
    return pk;

}

vector<mpz_class> enc(vector<mpz_class>& pk, std::vector<int>& m, int n = 80) {
    int d = 2 * n;
    mpz_class two_to_d = mpz_class(1) << d;
    mpz_class two_to_d_minus_two = mpz_class(1) << (d - 2);
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(time(NULL));
    mpz_class store_mod;
    mpz_class store_mult;
    vector<mpz_class> c(m.size(),0);
    mpz_class N;
    for (int i = 0; i < m.size(); i++) {
        // Generate a new N for every m
        N = rand.get_z_range(two_to_d - two_to_d_minus_two) + two_to_d_minus_two;
        mpz_mul(store_mult.get_mpz_t(),N.get_mpz_t(),pk[1].get_mpz_t());
        store_mult += m[i];
        mpz_mod(store_mod.get_mpz_t(),store_mult.get_mpz_t(),pk[0].get_mpz_t());
        c[i] = store_mod;
    }
    return c;
}

ZZ_mat<mpz_t> create_lattice(int t, vector<mpz_class>& c) {
    ZZ_mat<mpz_t> M(t, t+1);
    mpz_t one;
    mpz_init_set_ui(one, 1); 
    // #pragma omp target teams distribute parallel for
    // #pragma omp parallel for simd
    for (int i=0; i<t; i++) {
        mpz_t ci;
        mpz_init(ci);
        mpz_set(ci,c[i].get_mpz_t());
        mpz_mul(ci,ci,lambd.get_mpz_t());
        M(i, 0) = ci;
    }

    // Populate the rest of the matrix as identity matrix
    // #pragma omp parallel for collapse(2)
    for (int i=0; i<t; i++) {
        for (int j=1; j<t+1; j++) {
            if (j-1 == i) {
                M(i, j) = one;
            } else {
                M(i, j) = mpz_t{0};
            }
        }
    }

    return M;
}

ZZ_mat<mpz_t> create_lattice_prime(const ZZ_mat<mpz_t>& A, int t) {
    int rows = t;
    int cols = t - 3 + t;
    ZZ_mat<mpz_t> M(rows, cols);
    mpz_t one;
    mpz_init_set_ui(one, 1); 
    // #pragma omp parralel for simd collapse(2)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < t - 3; j++) {
            mpz_class val;
            A[j][i+1].get_mpz(val.get_mpz_t());
            mpz_t tmp;
            mpz_init(tmp);
            mpz_set(tmp,val.get_mpz_t());
            mpz_mul(tmp,tmp,lambd.get_mpz_t());
            M[i][j] = tmp;

            
        }
    }

    // Set the remaining entries to the identity matrix
    for (int i = 0; i < rows; i++) {
        for (int j = t - 3; j < cols; j++) {
            if (j - t + 3 == i) {
            
                M[i][j] = one;
            }
            else {
                M[i][j] = mpz_t{0};
            }
        }
    }

    return M;
}

bool attack(gmp_randstate_t state,int t){
    const int rho = 128;
    const int eta = 6400;
    const int gamma = 512000;

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);
    std::vector<int> m(t);
    auto start_time = std::chrono::steady_clock::now();

    for (int i = 0; i < t; ++i) {
        m[i] = dis(gen);
    }
     // print vector
    // for (int i = 0; i < t; ++i) {
    //     cout << m[i] << " ";
    // }
    cout << endl;

    vector<mpz_class> pk = genXi(rho,eta,gamma,state);

    vector<mpz_class> c = enc(pk,m);

    ZZ_mat<mpz_t> L = create_lattice(t,c);
    mpz_class val;
    auto lll1 = std::chrono::steady_clock::now();

    lll_reduction(L);

    auto lll1e = std::chrono::steady_clock::now();
    auto lllt = std::chrono::duration_cast<std::chrono::milliseconds>(lll1e - lll1).count();
    std::cout << "Time LLL 1: " << lllt << " ms" << std::endl;

    ZZ_mat<mpz_t> L_prime = create_lattice_prime(L,t);
    L.clear();

    auto lll2 = std::chrono::steady_clock::now();

    lll_reduction(L_prime);

    auto lll2e = std::chrono::steady_clock::now();
    auto lll2t = std::chrono::duration_cast<std::chrono::milliseconds>(lll2e - lll2).count();
    std::cout << "Time LLL 2: " << lll2t << " ms" << std::endl;
   // print the first row of the matrix A
   int cols = L_prime.get_cols() - L_prime.get_rows();
    for (int i = cols; i < L_prime.get_cols(); i++) {
        L_prime[0][i].get_mpz(val.get_mpz_t());
        bool res = val.get_ui() == m[i - cols];
        if(!res){
            cout << "Result: " << boolalpha << res << endl;
            cout << endl;
            return false;
        }
        // gmp_printf("%Zd ", val.get_mpz_t());

    }
    auto end_time = std::chrono::steady_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Total time: " << total_time << " ms" << std::endl;
    cout << endl;
    return true;

}

int main(){
    vector<int> t {5};
    std::random_device rd;  // obtain a random seed from the hardware
    std::mt19937 eng(rd());  // seed the generator
    std::uniform_int_distribution<> distr(100, 1000000000);  // define the range
    for(int j=0;j<t.size();j++){
        cout << "t: " << t[j] << endl;
        for(int i=0;i<1;i++){
            gmp_randstate_t state;
            gmp_randinit_mt(state);
            gmp_randseed_ui(state, distr(eng));
            // cout << "i: " << i  << endl;
            attack(state,t[j]);

        }
        cout << endl;
    }
   

}



