#ifndef PK_H
#define PK_H

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
    const unsigned int alpha;
    const unsigned int tau;
    mpz_class prime;
    mpz_class q_zero;
    mpz_class x_zero;
    mpz_class pk_sum;
    std::vector<mpz_class> chi;
    std::vector<mpz_class> r_i;
    std::vector<mpz_class> delta;
   

public:
    PKgenerator(unsigned int bits,unsigned int tau) :
        noise(54),
        eta(bits),
        gamma(150000),
        alpha(936),
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
        gmp_printf("parity: %Zd\n", test.get_mpz_t());

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
        mpz_mod(parity.get_mpz_t(),addition.get_mpz_t(),mpz_class(2).get_mpz_t());
        gmp_printf("parity sum: %Zd\n", parity.get_mpz_t());
        mpz_mod(parity.get_mpz_t(),reduced_c.get_mpz_t(),mpz_class(2).get_mpz_t());
        gmp_printf("parity reduced sum: %Zd\n", parity.get_mpz_t());


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

#endif

// big theta = 195
// small theta = 15
// mu = 56
