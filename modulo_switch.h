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
        std::vector<PKgenerator*> pk_list;
        mpz_class gamma{150000};
        unsigned int kappa;
        std::vector<mpf_class> y;
        const unsigned int theta;
        std::vector<bool> s;
        std::vector<mpz_class> q;
        mpz_class eta_prime;
    
    public:
        SwitchKey(unsigned int levels,unsigned int theta):
            levels {levels},
            theta(theta),
            y(theta,0),
            eta_ladder(levels,0),
            s(theta,false)
    {
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
        q.clear();
    }



    void generatePK(){
        for(int i = 0; i < eta_ladder.size(); i++){
           PKgenerator* pk = new PKgenerator(eta_ladder[i],158,54,150000);
           pk_list.push_back(pk);

        }
    }
    void decreasingModuli(){
        // create a list of decreasing eta values
        const unsigned int mu {56};

        eta_ladder[0] = (levels + 1) * mu;
        eta_ladder[eta_ladder.size() - 1] = 2 * mu;
        for (int i = levels-1; i > 1; i--) {
            unsigned int eta = (i + 1) * mu;
            eta_ladder[eta_ladder.size() - i] = eta;
        }
    }

    int switchKeyGen(unsigned int j){
        mpz_class eta = eta_ladder[1];
        eta_prime = eta_ladder[2];
        kappa = 2 * gamma.get_ui() + eta.get_ui();
        genYvector(eta_prime);
        genSvector(eta,eta_prime,pk_list[1]->getPrime(),pk_list[2]->getPrime(),pk_list[2]->getXZero(),pk_list[2]->getQZero());

        
        return 0;
    }

    mpz_class switchKey(mpz_class c){
        std::vector<mpz_class> expand_c(theta,0);
        std::vector<bool> bits;
        mpz_class modulus = mpz_class(1) << (eta_prime.get_ui() + 1);
        mpz_class c_mod;
        mpz_mod(c_mod.get_mpz_t(),c.get_mpz_t(),mpz_class(2).get_mpz_t());

        mpz_class mod_res;
        for(int i = 0; i < theta;i++){
            mpz_class product  = c * mpz_class(y[i].get_ui());
            // gmp_printf("product: %Zf\n", y[i].get_mpf_t());
            // mpf_out_str(stdout,10,10,y[i].get_mpf_t());

            mpz_mod(mod_res.get_mpz_t(),product.get_mpz_t(),modulus.get_mpz_t());
            expand_c[i] = mod_res;
        }

        bits.reserve(eta_prime.get_ui() + 1 * theta);
        for(auto elem : expand_c){
            // gmp_printf("elem: %Zd\n", elem.get_mpz_t());

            for (int i = eta_prime.get_ui(); i >= 0; --i) {

                bool bit = (elem.get_ui() >> i) & 1;
                bits.push_back(bit);
            }
        }

        // for(auto bit: bits){
        //     std::cout << "Bit:" << bit << std::endl;
        // }

        mpz_class c_prim = 0;
        for(int i = 0;i<q.size();i++){
            c_prim += q[i] * bits[i];
        }
        c_prim *= 2;
        // gmp_printf("cprim: %Zd\n", c_prim.get_mpz_t());  

        c_prim += c_mod;

        return c_prim;


    }

    int genYvector(mpz_class eta_prime){
        gmp_randstate_t state;
        gmp_randinit_default(state);
        gmp_randseed_ui(state, time(NULL));

        // Create the random number generator and seed it with the current time
        gmp_randclass rand(gmp_randinit_default);
        rand.seed(time(NULL));


        // Generate theta random float numbers with kappa bits of precision and store them in the y vector
        const unsigned int decimal_min = 1;
        const unsigned int decimal_max = (2 << (eta_prime.get_ui() + 1)) - 1;
        
        for (int i = 0; i < theta; i++) {
            // Generate a random decimal part in the specified range
            mpz_class decimal_part = rand.get_z_range(decimal_max - decimal_min + 1) + decimal_min;

            // Generate a random float with the specified precision and add the decimal part to it
            mpf_class f = rand.get_f(kappa) + mpf_class(decimal_part);
            y[i]=f;
        }
        gmp_randclear(state);
        return 0;
    }

    void genSvector(mpz_class eta,mpz_class eta_prime,mpz_class p,mpz_class p_prime,mpz_class x0,mpz_class q0){
        gmp_randstate_t state;
        gmp_randinit_default(state);
        gmp_randseed_ui(state, time(NULL));
        gmp_randinit_mt(state);
        mpf_set_default_prec(kappa);

        // Compute the divisor 2^eta / prime
        mpf_class divisor;
        mpf_class top(mpf_class(1) << eta_prime.get_ui());
        mpf_div(divisor.get_mpf_t(),top.get_mpf_t(),mpf_class(p).get_mpf_t());
     


        // Compute the threshold
        mpf_class eps;
        mpf_class power(mpf_class(1) << kappa);
        mpf_div(eps.get_mpf_t(),mpf_class(1).get_mpf_t(),power.get_mpf_t());

        mpf_class threshold = divisor - eps;
        mpz_class modu(mpz_class(1) << eta_prime.get_ui() + 1);
        // Generate the bit vector s
    
        mpz_class sum {0};
        for (int i = 0; i < y.size(); i++) {
            mpz_class decimal{y[i].get_ui()};
            mpz_class rem;
            mpz_mod(rem.get_mpz_t(),decimal.get_mpz_t(),modu.get_mpz_t());
            if (rem > threshold || mpf_cmp(mpf_class(sum).get_mpf_t(),threshold.get_mpf_t())>0) {
                s[i] = 0;
            } else {
                s[i] = 1;
                sum += y[i].get_ui();
                mpz_mod(sum.get_mpz_t(),sum.get_mpz_t(),modu.get_mpz_t());

            }
        }
        
        mpz_class rho{54};
        const unsigned int range = (eta_prime.get_ui() + 1) * theta;
        q.reserve(range);
        std::vector<mpz_class> r(range,0);
        for(int i=0; i<range;i++){
            mpz_class q_i{0};
            mpz_urandomm(q_i.get_mpz_t(),state,q0.get_mpz_t());
            q.push_back(q_i);
        }

        for(int i=0;i<range;i++){
            mpz_class r_i{0};
            mpz_urandomb(r_i.get_mpz_t(),state,rho.get_ui());
            r[i] = r_i;
        }

        for(int i=0; i< range;i++){
            q[i] *= p_prime;
            q[i] += r[i];
        }

        mpz_class diviso(mpz_class(1) << (eta_prime.get_ui() + 1));
        mpz_class scal {p_prime / diviso};
        std::vector<mpz_class> s_prime = powersOf2(eta_prime,range);
        for(int i=0;i<range;i++){
            s_prime[i] *= scal;
            q[i] += s_prime[i];
        }

        gmp_randclear(state);
        
        r.clear();

        return;
    }

    std::vector<mpz_class> powersOf2(mpz_class eta_p,unsigned int range){
        int k {1};
        std::vector<mpz_class> s_prime(range,0);
        for (int i = 1; i <=eta_p; i++) {
            unsigned int cnt = 0;
            for(int j=i*theta;j<(i+1)*theta;j++){
                s_prime[j] = s[cnt] * (k << i);
                cnt+=1;
            }
        }
        return s_prime;

    }

    mpz_class enc(mpz_class m){
        return pk_list[1]->encrypt(m);
    }
    
    mpz_class dec(mpz_class c,int i){
        return pk_list[2]->decrypt(c);
    }

    mpz_class add(mpz_class c1, mpz_class c2){
        mpz_class parity;
        mpz_class reduced;
        mpz_class c3 = c1 + c2;
    

        mpz_mod(reduced.get_mpz_t(),c3.get_mpz_t(),pk_list[1]->getXZero().get_mpz_t());
        // mpz_mod(parity.get_mpz_t(), reduced.get_mpz_t(), mpz_class(2).get_mpz_t());
        // gmp_printf("parity of reduced c3: %Zd\n", parity.get_mpz_t());

        mpz_class c_prim = switchKey(c3);
        mpz_class m = dec(c_prim,1);
        return m; 
    }

};
#endif