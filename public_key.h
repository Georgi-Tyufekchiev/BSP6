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
    unsigned int delta_size;

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
    std::vector<unsigned long int> r_i;
    std::vector<unsigned long int> delta;
   

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
        rand(gmp_randinit_default),
        delta_size{0}
        
    {
        gmp_randinit_default(state);
        rand.seed(rand_seed);

        // Generate random prime integer p of size eta bits
        auto start_time = std::chrono::steady_clock::now();

        mpz_urandomb(prime.get_mpz_t(), state, eta);
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());

        // Upper bound for q0 which is 2^y / prime
        mpz_class bound(mpz_class(2) << gamma - eta);
        mpz_urandomm(q_zero.get_mpz_t(), state, bound.get_mpz_t());
        x_zero = (2*q_zero + 1) * prime;
        genXi();
        chi.clear();
        r_i.clear();

        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        // std::cout << "KeyGen: " << total_time << " ms" << std::endl;
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
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> distr(1,noise); // distribution in range [1, 6]

        for (int i = 0; i < tau; i++) {
            // Generate integers r_i of size ( 2^2*noise) bits
            r_i[i] = distr(rng);
        }

        for (int i = 0; i < tau; i++) {
            // Generate a set of integers Chi_i orf size [2^y] bits
            chi[i] = rand.get_z_bits(gamma);
        }
     
        for (int i = 0; i < tau; i++){
            // Compute the  deltas as delta = Chi mod p + xi * prime - r
            tmp_delta = chi[i] - 2*r_i[i];
            mpz_mod(delta_mod.get_mpz_t(), tmp_delta.get_mpz_t(), prime.get_mpz_t());
            delta[i] = delta_mod.get_ui();
            int size = mpz_sizeinbase(delta_mod.get_mpz_t(), 2);
            delta_size += size;
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
        gmp_randinit_default(state);
        rand.seed(rand_seed);
        for(int i = 0; i < tau; i++){
            epsilon = rand.get_z_range(bound);
            x = rand.get_z_bits(gamma) - delta[i];
            pk_sum += x * epsilon;
        }

        return;
    }

    mpz_class encrypt(mpz_class m){
        mpz_class r;
        mpz_class ciphertext_mod{0};
        mpz_class ciphertext;
        auto start_time = std::chrono::steady_clock::now();

        r = rand.get_z_bits(noise);
        computePKsum();
        ciphertext = 2*pk_sum + 2 * r + m;
        mpz_mod(ciphertext_mod.get_mpz_t(), ciphertext.get_mpz_t(), x_zero.get_mpz_t());
        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Encrypt: " << total_time << " ms" << std::endl;
        return ciphertext;
    }

    mpz_class decrypt(mpz_class c){
        mpz_class plaintext;
        mpz_class tmp;
        mpz_t modulus;
        mpz_init_set_ui(modulus,2);
        auto start_time = std::chrono::steady_clock::now();

        mpz_mod(plaintext.get_mpz_t(), c.get_mpz_t(), prime.get_mpz_t());
        mpz_mod(tmp.get_mpz_t(), plaintext.get_mpz_t(), modulus);
        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "Decrypt: " << total_time << " ms" << std::endl;
        gmp_printf("bit: %Zd\n", tmp.get_mpz_t());
        return tmp;

    }

    mpz_class add(mpz_class c1,mpz_class c2){
        mpz_class reduced_c;
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

            mpz_class addition = decrypt(add((m1),encrypt(m2)));
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
    std::vector<unsigned long int> getDelta() const { return delta; }
    unsigned int pksize(){
        unsigned int size = mpz_sizeinbase(x_zero.get_mpz_t(), 2);
        unsigned int total_size = (size + delta_size) / 8000;
        return total_size;

    }

};


#endif

  
