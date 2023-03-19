#ifndef MODULO_H
#define MODULO_H

#include <vector>
#include <ctime>
#include <cmath>
#include <assert.h>
#include <random>
#include <iostream>
#include <stdexcept> 
#include <gmpxx.h>
#include "public_key.h"


class SwitchKey{
    private:
        const unsigned int levels;
        std::vector<unsigned int> eta_ladder;
        const unsigned int mu {56};
        std::vector<PKgenerator*> pk_list;
        mpz_class gamma{150000};
        unsigned int kappa;
        std::vector<mpz_class> y;
        const unsigned int theta;
        std::vector<bool> s;
        std::vector<bool> newS;
        std::vector<bool> s_prime;
        std::vector<mpz_class> small_sigma;

    public:
        SwitchKey(unsigned int levels,unsigned int theta):
            levels {levels},
            theta(theta),
            y(theta,0),
            eta_ladder(levels,0),
            s(theta,false),
            newS(theta,false),
            small_sigma(theta,0),
            s_prime(theta,false)         
    {
        // s.reserve(theta);
        // newS.reserve(theta);
        // s_prime.reserve(theta);
        // small_sigma.reserve(theta);
        decreasingModuli();
        generatePK();
    }

        ~SwitchKey() {
        for (auto pk : pk_list) {
            delete pk;
        }
        pk_list.clear();
        y.clear();
        eta_ladder.clear();
        s.clear();
        newS.clear();
        small_sigma.clear();
        s_prime.clear();
    }



 void generatePK(){
        // pk_list.reserve(levels);

        for(int i = 0; i < eta_ladder.size(); i++){
           PKgenerator* pk = new PKgenerator(eta_ladder[i],158);
           pk_list.push_back(pk);

        }
    }
    void decreasingModuli(){
        // create a list of decreasing eta values
        // eta_ladder.reserve(levels);

        eta_ladder[0] = (levels + 1) * mu;
        eta_ladder[eta_ladder.size() - 1] = 2 * mu;
        for (int i = levels-1; i > 1; i--) {
            unsigned int eta = (i + 1) * mu;
            eta_ladder[i-1] = eta;
        }
    }

    

    int switchKeyGen(unsigned int j){
        mpz_class eta = eta_ladder[j-1];
        mpz_class eta_prime = eta_ladder[j-2];
        kappa = 2 * gamma.get_ui() + eta.get_ui();
        genYvector(eta_prime);
        genSvector(eta,eta_prime,pk_list[0]->getPrime(),pk_list[1]->getPrime(),pk_list[1]->getXZero());
        return 0;
    }

    int genYvector(mpz_class eta_prime){
        gmp_randstate_t state;
        gmp_randinit_default(state);
        gmp_randseed_ui(state, time(NULL));
        const unsigned int base = 10;
        mpz_class modulus = mpz_class(1) << (eta_prime.get_ui() + 1);

        // generate random numbers modulo 2^(eta_prime + 1)
        // each number has kappa bits precision after the binary point
        mpf_set_default_prec(kappa);
        for (int i = 0; i < theta; i++) {
            mpz_class randomInt;
            mpz_urandomb(randomInt.get_mpz_t(), state, kappa);
            mpz_mod(randomInt.get_mpz_t(), randomInt.get_mpz_t(), modulus.get_mpz_t());
            y[i] = randomInt;
        }

        gmp_randclear(state);
        return 0;
    }

    int genSvector(mpz_class eta,mpz_class eta_prime,mpz_class p,mpz_class p_prime,mpz_class x0){
        gmp_randstate_t state;
        gmp_randinit_default(state);
        gmp_randseed_ui(state, time(NULL));
        gmp_randinit_mt(state);
        mp_exp_t* exponent = new mp_exp_t(1);
        
        mpz_class modulus = mpz_class(1) << (eta_prime.get_ui() + 1);


        // Generate a random vector s of theta bits
        for (int i = 0; i < theta; i++) {
            mp_bitcnt_t randomBit = gmp_urandomb_ui(state, 1);
            s[i] = (randomBit == 1);
        }
        

        // Calculate the dot product of s and y
        mpz_class dotProduct = 0;
        for (int i = 0; i < theta; i++) {
            if (s[i]) {
                dotProduct += y[i];
            }
        }

        // Calculate the error term e
        mpf_class floatError = mpf_class(0, kappa);
        mpf_urandomb(floatError.get_mpf_t(), state, kappa);
        char* strError = new char[kappa + 2];
        mpf_get_str(strError, exponent, 10, kappa, floatError.get_mpf_t());
        mpz_class error = mpz_class(strError) % (mpz_class(1) << eta.get_ui());
        error -= (mpz_class(1) << (eta.get_ui() - 1));
        delete[] strError;
        delete exponent;
        // printf("Error term ok\n");
        
        // Calculate the left-hand side of the equation
        mpz_class lhs = (mpz_class(1) << eta.get_ui()) / p;
        
        // Calculate the right-hand side of the equation
        mpz_class rhs = dotProduct + error;
        rhs %= modulus;

         // Calculate the vector s such that (2^eta)/p = <s, y> + e mod 2^(eta+1)
        for (int i = 0; i < theta; i++) {
            if (rhs >= (modulus / mpz_class(2))) {
                newS[i] = true;
                rhs -= y[i];
            } else {
                newS[i] = false;
            }
            rhs *= 2;
        }

        mpz_class q0 {x0 / p_prime};
        mpz_class rho{mpz_class(1) << 54};
        const unsigned int range = (eta_prime.get_ui() + 1) * theta;
        std::vector<mpz_class> q(range,0);
        std::vector<mpz_class> r(range,0);
        for(int i=0; i<range;i++){


            mpz_class q_i{0};

            mpz_urandomm(q_i.get_mpz_t(),state,q0.get_mpz_t());

            q[i] = q_i;
        }

        for(int i=0;i<range;i++){
            mpz_class r_i{0};
            
            mpz_urandomm(r_i.get_mpz_t(),state,rho.get_mpz_t());
            r[i] = r_i;
        }

        gmp_randclear(state);
        q.clear();
        r.clear();

        return 0;
    }

    void powersOf2(mpz_class eta){
        int k {1};
        for (int i = 0; i <= eta; i++) {
            s_prime.push_back(s[i] * (k << i));
        }
    }

};
#endif