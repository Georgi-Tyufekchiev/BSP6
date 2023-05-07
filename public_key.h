#ifndef PK_H
#define PK_H
#include <chrono>

#include <gmpxx.h>
#include <vector>
#include <ctime>
#include <cmath>
#include <assert.h>
#include <random>
#include <iostream>

class PKgenerator {
public:
    gmp_randclass rand;
    const time_t rand_seed {time(NULL)};
    gmp_randstate_t state;

private:
    const unsigned int noise;
    const unsigned int eta;
    const unsigned int gamma;
    const unsigned int tau;
    mpz_class prime;
    mpz_class q_zero;
    mpz_class x_zero;
    mpz_class pk_sum;
    std::vector<mpz_class> chi;
    std::vector<mpz_class> r_i;
    std::vector<mpz_class> delta;
   

public:
    PKgenerator(unsigned int bits,unsigned int tau,unsigned int rho,unsigned int gamma) :
        noise(rho),
        eta(bits),
        gamma(gamma),
        tau(tau),
        prime(0),
        q_zero(0),
        x_zero(0),
        pk_sum(0),
        chi(tau,0),
        r_i(tau,0),
        delta(tau,0),
        rand(gmp_randinit_default)
        
    {
        gmp_randinit_default(state);
        gmp_randseed_ui(state, time(NULL));

        rand.seed(rand_seed);

        // Generate random prime integer p of size eta bits
        mpz_urandomb(prime.get_mpz_t(), state, eta);
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());

        // Upper bound for q0 which is 2^y / prime
        mpz_class bound_y(mpz_class(1) << gamma);
        mpz_class bound_div;
        mpz_cdiv_q(bound_div.get_mpz_t(),bound_y.get_mpz_t(),prime.get_mpz_t());
        q_zero = rand.get_z_range(bound_div.get_ui());
        // mpz_urandomb(q_zero.get_mpz_t(), state,bound_div.get_ui());
        x_zero = (2*q_zero + 1) * prime;
        genXi();
        gmp_randclear(state);
    }
        ~PKgenerator(){
            chi.clear();
            r_i.clear();
            delta.clear();
        }

    void genXi(){
        mpz_class delta_mod;
        mpz_class tmp_delta;

        for (int i = 0; i < tau; i++) {
            // Generate integers r_i of size ( 2^2*noise) bits
            r_i[i] = rand.get_z_bits(noise);
        }

        for (int i = 0; i < tau; i++) {
            // Generate a set of integers Chi_i orf size [2^y] bits
            chi[i] = rand.get_z_bits(gamma);
        }
     

        for (int i = 0; i < tau; i++){
            // Compute the  deltas as delta = Chi mod p + xi * prime - r
            tmp_delta = chi[i] - 2*r_i[i];
            mpz_mod(delta_mod.get_mpz_t(), tmp_delta.get_mpz_t(), prime.get_mpz_t());
            delta[i] = delta_mod;
        }
    }

    std::vector<mpz_class> getXi(){
        std::vector<mpz_class> x(chi.size(),0);
        for(int i = 0; i < tau; i++){
            x[i] = chi[i] - delta[i];
        }

        return x;
    }


    void computePKsum(){
        mpz_class bound {2};
        mpz_class x;
        mpz_class epsilon;
        for(int i = 0; i < tau; i++){
            epsilon = rand.get_z_range(bound);
            x = chi[i] - delta[i];
            pk_sum += x * epsilon;
        }

        return;
    }


    mpz_class encrypt(mpz_class m){
        mpz_class r;
        mpz_class ciphertext_mod{0};
        mpz_class ciphertext;
        mpz_class test;

        r = rand.get_z_bits(noise);
        computePKsum();
        ciphertext = 2*pk_sum + 2 * r + m;
        mpz_mod(ciphertext_mod.get_mpz_t(), ciphertext.get_mpz_t(), x_zero.get_mpz_t());
        mpz_mod(test.get_mpz_t(), ciphertext_mod.get_mpz_t(), mpz_class(2).get_mpz_t());

        return ciphertext;

    }

    mpz_class decrypt(mpz_class c){
        mpz_class plaintext;
        mpz_class tmp;
        mpz_class modulus {2};
        mpz_mod(plaintext.get_mpz_t(), c.get_mpz_t(), prime.get_mpz_t());
        mpz_mod(tmp.get_mpz_t(), plaintext.get_mpz_t(), modulus.get_mpz_t());
        gmp_printf("bit: %Zd\n", tmp.get_mpz_t());
        return tmp;

    }

    mpz_class add(mpz_class c1,mpz_class c2){
        mpz_class reduced_c;
        mpz_class parity;
        mpz_class addition = c1 + c2;
        mpz_mod(reduced_c.get_mpz_t(),addition.get_mpz_t(),x_zero.get_mpz_t());


        return reduced_c;
    }

    mpz_class mutliply(mpz_class c1, mpz_class c2){
        mpz_class reduced_c;
        mpz_class multiplication = c1 * c2;
        mpz_mod(reduced_c.get_mpz_t(),multiplication.get_mpz_t(),x_zero.get_mpz_t());
        return reduced_c;
    }

    void checkOperations(){
        mpz_class result;
        for(int i=0;i<3;i++){
            mpz_class m1 = rand.get_z_range(mpz_class(2));
            mpz_class m2 = rand.get_z_range(mpz_class(2));

            mpz_class addition = decrypt(add(encrypt(m1),encrypt(m2)));
            mpz_class a = m1+m2;
            mpz_mod(result.get_mpz_t(),a.get_mpz_t(),mpz_class(2).get_mpz_t());
            
            assert(addition == result);
        }

        for(int i=0;i<3;i++){
            mpz_class m1 = rand.get_z_range(mpz_class(2));
            mpz_class m2 = rand.get_z_range(mpz_class(2));

            mpz_class mult = decrypt(mutliply(encrypt(m1),encrypt(m2)));
            mpz_class a = m1*m2;
            mpz_mod(result.get_mpz_t(),a.get_mpz_t(),mpz_class(2).get_mpz_t());
            
            assert(mult == result);
        }
        return;
    }




        
    mpz_class getPrime() const { return prime; }
    mpz_class getXZero() const { return x_zero; }
    mpz_class getQZero() const { return q_zero;}
    time_t getSeed() const {printf("seed: %s\n", ctime(&rand_seed));return rand_seed; }
    std::vector<mpz_class> getDelta() const { return delta; }

    void printPK(){
        printf("seed: %s\n", ctime(&rand_seed));
        char* str = new char[mpz_sizeinbase(x_zero.get_mpz_t(), 10) + 2];
        mpz_get_str(str, 10, x_zero.get_mpz_t());
        printf("x0: %s\n",str);
        for(mpz_class v: delta){
        char* str = new char[mpz_sizeinbase(v.get_mpz_t(), 10) + 2];
        mpz_get_str(str, 10, v.get_mpz_t());
        printf("delta: %s\n",str);
        }
        delete[] str;
        return;
    }
};


mpz_class prime_to_m(mpz_class n, mpz_class m){
    mpz_class g;
    while(true){
        mpz_gcd(g.get_mpz_t(), n.get_mpz_t(), m.get_mpz_t());
        if(g == 1){
            return n;
        }else{
            mpz_fdiv_q(n.get_mpz_t(),n.get_mpz_t(),g.get_mpz_t());
        }
        
    }
}
mpz_class xi(mpz_class p,int gam,int eta,int rho,gmp_randstate_t state){


    mpz_class bound_q(mpz_class(1) << (gam -eta));
    mpz_class q;
    mpz_class r;
    mpz_urandomm(q.get_mpz_t(),state,bound_q.get_mpz_t());
    mpz_class bound_r = mpz_class(1) << rho;
    mpz_urandomm(r.get_mpz_t(),state,bound_r.get_mpz_t());

    return (p * q) + r;
}


bool attackGACD(){
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, 1234);
    const unsigned int rho = 12;
    const unsigned int gamma = 1000;
    const unsigned int eta = 100;
    mpz_class prime{"1133866015397999278511396345017"};
    // mpz_urandomb(prime.get_mpz_t(), state, eta);
    // mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());

    unsigned int pBits = mpz_sizeinbase(prime.get_mpz_t(), 2);
    gmp_printf("prime: %Zd\n", prime.get_mpz_t());
    printf("p size:%d \n",pBits);
    int exp = (rho * (rho +1)) / (rho-1);
    mpz_class B;
    mpz_pow_ui(B.get_mpz_t(),mpz_class(2).get_mpz_t(),exp);
    B += 2200;

    mpz_class fa{1};
    mpz_fac_ui(fa.get_mpz_t(),B.get_ui());

    mpz_class g;
    for(int j=1;j<rho;j++){
        mpz_class z{1};
        for(int i=0;i< (1 << rho);i++){
            z = z * (xi(prime,gamma,eta,rho,state) - i);
        }

        if(j == 1){
             g = z;
            continue;
        }

        mpz_class gcd;
        mpz_gcd(gcd.get_mpz_t(), g.get_mpz_t(), z.get_mpz_t());
        g= prime_to_m(gcd,fa);
        gmp_printf("g: %Zd\n", g.get_mpz_t());

        // printf("j: %d \n",j);
        unsigned int gBits = mpz_sizeinbase(g.get_mpz_t(), 2);
        printf("gcd size:%d \n",gBits);

        if(gBits == pBits || gBits < pBits){

            break;
        }
        
    }

    if((mpz_cmp(g.get_mpz_t(),prime.get_mpz_t())) == 0){
        return true;
    }else{
        return false;
    }

}

#endif

  
