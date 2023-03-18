#ifndef MODULO_H
#define MODULO_H

#include "public_key.h"

class SwitchKey{
    private:
        const unsigned int levels;
        std::vector<unsigned int> eta_ladder;
        const unsigned int mu {56};
        std::vector<PKgenerator> pk_list;
        mpz_class gamma;
        unsigned int kappa;
        std::vector<mpz_class> y;
        const unsigned int theta;
        std::vector<bool> s;
    
    public:
        SwitchKey(unsigned int L):
            levels {L},
            eta_ladder(levels),
            gamma {150000},
            theta {156},
            y(theta),
            s(theta)
    {
        pk_list.reserve(levels);
    }

    void generatePK(){
        for(int i = 0; i < eta_ladder.size(); i++){
           PKgenerator pk(eta_ladder[i]);
           pk_list.push_back(pk);

        }
    }

    void decreasingModuli(){
        eta_ladder[eta_ladder.size() - 1] = (levels + 1) * mu;
        mpz_class eta_ladder[0] = 2 * mu;
        for (int i = levels-1; i > 1; i--) {
            mpz_class eta = (i + 1) * mu;
            eta_ladder[i-1] = eta;
        }
    }

    void switchKeyGen(unsigned int j){
        mpz_class eta = eta_ladder[j];
        mpz_class eta_prime = eta_ladder[j-1];
        kappa = 2 * gamma.get_ui() + eta.get_ui();
    }

    void genYvector(mpz_class eta_prime){
        gmp_randstate_t state;
        gmp_randinit_default(state);
        gmp_randseed_ui(state, time(NULL));
        mp_exp_t* expo;
        size_t base(10);
        mpz_class modulus = mpz_class(1) << (eta_prime.get_ui() + 1);

        // generate random numbers modulo 2^(eta_prime + 1)
        // each number has kappa bits precision after the binary point
        for (int i = 0; i < theta; i++) {
            mpf_class floatRandom = mpf_class(0, kappa);
            mpf_urandomb(floatRandom.get_mpf_t(), state, kappa);
            std::string strRandom = mpf_get_str(NULL, expo, base,kappa, floatRandom.get_mpf_t());
            mpz_class random = mpz_class(strRandom.c_str()) % modulus;
            y[i] = random % modulus;
        }

        gmp_randclear(state);
    }

    void genSvector(mpz_class eta,mpz_class eta_prime,mpz_class p){
        gmp_randstate_t state;
        gmp_randinit_default(state);
        gmp_randseed_ui(state, time(NULL));
        gmp_randinit_mt(state);
        mp_exp_t exponent;
        mpz_class modulus = mpz_class(1) << (eta_prime.get_ui() + 1);


        // Generate a random vector s of theta bits
        for (int i = 0; i < theta; i++) {
            mp_bitcnt_t randomBit = gmp_urandomb_ui(state, 1);
            s[i] = (randomBit == 1);
        }
        
        gmp_randclear(state);

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
        mpf_get_str(strError, &exponent, 10, kappa, floatError.get_mpf_t());
        mpz_class error = mpz_class(strError) % (mpz_class(1) << eta.get_ui());
        error -= (mpz_class(1) << (eta.get_ui() - 1));
        delete[] strError;
        
        // Calculate the left-hand side of the equation
        mpz_class lhs = (mpz_class(1) << eta.get_ui()) / p;
        
        // Calculate the right-hand side of the equation
        mpz_class rhs = dotProduct + error;
        rhs %= modulus;

         // Calculate the vector s such that (2^eta)/p = <s, y> + e mod 2^(eta+1)
        std::vector<bool> newS(theta);
        for (int i = 0; i < theta; i++) {
            if (rhs >= (modulus / mpz_class(2))) {
                newS[i] = true;
                rhs -= y[i];
            } else {
                newS[i] = false;
            }
            rhs *= 2;
        }

    }

};


#endif