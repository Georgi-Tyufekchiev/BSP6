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
#include <thread> // for sleep_for


class SwitchKey{
    private:
        const unsigned int levels;
        std::vector<unsigned int> eta_ladder;
        std::vector<PKgenerator*> pk_list;
        const unsigned int gamma;
        unsigned int kappa;
        std::vector<mpf_class> y;
        const unsigned int theta;
        std::vector<bool> s;
        std::vector<unsigned long int> q;
        mpz_class eta_prime;
        int rho;
        unsigned int size;
        const unsigned int mu;

    
    public:
        SwitchKey(unsigned int levels,unsigned int theta, int rho, int gamma, int mu):
            levels {levels},
            theta(theta),
            y(theta,0),
            eta_ladder(levels,0),
            s(theta,false),
            rho(rho),
            gamma(gamma),
            size(0),
            mu(mu)
    {
        auto start_time = std::chrono::steady_clock::now();

        decreasingModuli();
        generatePK();

        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "KeyGen: " << total_time  << " ms" << std::endl;
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
            PKgenerator* pk = new PKgenerator(eta_ladder[i],158,rho,gamma);
            size += pk->pksize();
            pk_list.push_back(pk);
            // std::this_thread::sleep_for(std::chrono::seconds(1)); 

        }
        std::cout << "Pk Size: " << size << std::endl;
    }
    void decreasingModuli(){
        // create a list of decreasing eta values

        eta_ladder[0] = (levels + 1) * mu;
        eta_ladder[eta_ladder.size() - 1] = 2 * mu;
        for (int i = levels-1; i > 1; i--) {
            unsigned int eta = (i + 1) * mu;
            eta_ladder[eta_ladder.size() - i] = eta;
        }
    }

    int switchKeyGen(unsigned int j){
        auto start_time = std::chrono::steady_clock::now();
        mpz_class eta = eta_ladder[1];
        eta_prime = eta_ladder[2];
        kappa = 2 * gamma + eta.get_ui();
        genYvector(eta_prime);
        genSvector(eta,eta_prime,pk_list[1]->getPrime(),pk_list[2]->getPrime(),pk_list[2]->getXZero(),pk_list[2]->getQZero());
        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
        std::cout << "SwitchKeyGen: " << total_time * j << " s" << std::endl;
        return 0;
    }

    mpz_class switchKey(mpz_class c){
        auto start_time = std::chrono::steady_clock::now();

        std::vector<unsigned long int> expand_c(theta,0);
        std::vector<bool> bits;
        mpz_class modulus = mpz_class(1) << (eta_prime.get_ui() + 1);
        mpz_class c_mod;
        mpz_mod(c_mod.get_mpz_t(),c.get_mpz_t(),mpz_class(2).get_mpz_t());

        mpz_class mod_res;
        for(int i = 0; i < theta;i++){
            mpz_class product  = c * mpz_class(y[i].get_ui());
            mpz_mod(mod_res.get_mpz_t(),product.get_mpz_t(),modulus.get_mpz_t());
            expand_c[i] = mod_res.get_ui();
        }
        y.clear();

        bits.reserve(eta_prime.get_ui() + 1 * theta);
        for(auto elem : expand_c){
            for (int i = eta_prime.get_ui(); i >= 0; --i) {

                bool bit = (elem >> i) & 1;
                bits.push_back(bit);
            }
        }
        expand_c.clear();
        mpz_class c_prim = 0;
        for(int i = 0;i<q.size();i++){
            c_prim += q[i] * bits[i];
        }
        q.clear();
        bits.clear();
        c_prim *= 2;
        c_prim += c_mod;
        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        std::cout << "SwitchKey: " << total_time << " ms" << std::endl;
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
        
        const unsigned int range = (eta_prime.get_ui() + 1) * theta;
        q.reserve(range);
        std::vector<unsigned long int> r(range,0);
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist(1,q0.get_ui()); // distribution in range [1, 6]
        for(int i=0; i<range;i++){
            // mpz_class q_i{0};
            // mpz_urandomm(q_i.get_mpz_t(),state,q0.get_mpz_t());
            q.push_back(dist(rng));
        }

        std::uniform_int_distribution<std::mt19937::result_type> distr(1,rho); // distribution in range [1, 6]

        for(int i=0;i<range;i++){
            // mpz_class r_i{0};
            // mpz_urandomb(r_i.get_mpz_t(),state,rho);
            r[i] = distr(rng);
        }

        for(int i=0; i< range;i++){
            q[i] *= p_prime.get_ui();
            q[i] += r[i];
        }

        mpz_class diviso(mpz_class(1) << (eta_prime.get_ui() + 1));
        mpz_class scal {p_prime / diviso};
        std::vector<unsigned long int> s_prime = powersOf2(eta_prime,range);
        for(int i=0;i<range;i++){
            s_prime[i] *= scal.get_ui();
            q[i] += s_prime[i];
        }


        gmp_randclear(state);
        s_prime.clear();
        r.clear();

        return;
    }

    std::vector<unsigned long int> powersOf2(mpz_class eta_p,unsigned int range){
        int k {1};
        std::vector<unsigned long int> s_prime(range,0);
        for (int i = 1; i <=eta_p; i++) {
            unsigned int cnt = 0;
            for(int j=i*theta;j<(i+1)*theta;j++){
                s_prime[j] = s[cnt] * (k << i);
                cnt+=1;
            }
        }
        s.clear();
        return s_prime;
    }

    mpz_class enc(mpz_class m){return pk_list[1]->encrypt(m);}
    
    unsigned long dec(mpz_class c){return pk_list[2]->decrypt(c);}

    mpz_class XOR_gate(const mpz_class& ct1, const mpz_class& ct2){
        mpz_class c;
        mpz_add(c.get_mpz_t(),ct1.get_mpz_t(),ct2.get_mpz_t());
        mpz_class new_c = switchKey(c);
        return new_c;
    }

    mpz_class AND_gate(const mpz_class& ct1, const mpz_class& ct2){
        mpz_class c;
        mpz_mul(c.get_mpz_t(),ct1.get_mpz_t(),ct2.get_mpz_t());
        mpz_class new_c = switchKey(c);
        return new_c;
    }

    mpz_class OR_gate(const mpz_class& ct1, const mpz_class& ct2){
        mpz_class one = pk_list[1]->encrypt(mpz_class(1));
        mpz_class result = XOR_gate(ct1,one);
        mpz_class tmp = XOR_gate(ct2,one);
        result = AND_gate(result,tmp);   
        result = XOR_gate(result,one);
        return result;  
    }

    mpz_class NOT_gate(const mpz_class& ct1){
        mpz_class one = pk_list[1]->encrypt(mpz_class(1));
        mpz_class result = XOR_gate(ct1,one);
        return result;
    }

    std::vector<mpz_class> add(std::vector<mpz_class> cts){
        std::vector<mpz_class> new_cts;
        new_cts.reserve(cts.size());
        if(cts.size() == 2){
            new_cts.emplace_back(AND_gate(cts[0],cts[1]));
            new_cts.emplace_back(XOR_gate(cts[0],cts[1]));
        }

        if(cts.size() == 3){
            mpz_class tmp = OR_gate(cts[0],cts[1]);
            mpz_class tmp2 = OR_gate(cts[1],cts[2]);
            new_cts.emplace_back(AND_gate(tmp,tmp2));
            tmp = XOR_gate(cts[0],cts[1]);
            new_cts.emplace_back(XOR_gate(tmp,cts[2]));
        }
   

        return new_cts; 
    }

};
#endif