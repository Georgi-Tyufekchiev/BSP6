
// #include <gmpxx.h>
// #include <vector>
// #include <ctime>
// #include <cmath>
// #include <assert.h>
// #include <random>
// #include <iostream>
// class PKgenerator {
// private:
//     const unsigned int noise;
//     const unsigned int eta;
//     const unsigned int gamma;
//     const unsigned int alpha;
//     const unsigned int tau;
//     const unsigned int levels;
//     mpz_class prime;
//     mpz_class q_zero;
//     mpz_class x_zero;
//     mpz_class pk_sum;
//     std::vector<mpz_class> chi;
//     std::vector<mpz_class> r_i;
//     std::vector<mpz_class> delta;
//     gmp_randclass rand;
//     const time_t rand_seed {time(NULL)};
//     gmp_randstate_t state;
//     std::vector<mpz_class> ladder;
// public:
//     PKgenerator(unsigned int bits) :
//         noise(54),
//         eta(bits*3),
//         gamma(150000),
//         alpha(936),
//         tau(158),
//         prime(0),
//         q_zero(0),
//         x_zero(0),
//         pk_sum(0),
//         chi(tau),
//         r_i(tau),
//         delta(tau),
//         rand(gmp_randinit_default),
//         levels(5),
//         ladder(levels,0)
//     {
//         gmp_randinit_default(state);
//         rand.seed(rand_seed);

//         // Generate random prime integer p of size eta bits
//         mpz_urandomb(prime.get_mpz_t(), state, eta);
//         mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());

//         // Upper bound for q0 which is 2^y / prime
//         mpz_urandomb(q_zero.get_mpz_t(), state, gamma-1);
//         x_zero = (2*q_zero + 1) * prime;
//         genXi();
//         // gmp_randclear(state);
//     }

//     void genXi(){
//         mpz_class delta_mod;
//         mpz_class tmp_delta;

//         for (int i = 0; i < tau; i++) {
//             // Generate integers r_i of size ( 2^2*noise) bits
//             r_i[i] = rand.get_z_bits(noise);
//         }

//         for (int i = 0; i < tau; i++) {
//             // Generate a set of integers Chi_i orf size [2^y] bits
//             chi[i] = rand.get_z_bits(gamma);
//         }
     

//         for (int i = 0; i < tau; i++){
//             // Compute the  deltas as delta = Chi mod p + xi * prime - r
//             tmp_delta = chi[i] - 2*r_i[i];
//             mpz_mod(delta_mod.get_mpz_t(), tmp_delta.get_mpz_t(), prime.get_mpz_t());
//             delta[i] = delta_mod;
//         }
//     }

//     void computePKsum(){
//         mpz_class bound {2};
//         mpz_class x;
//         mpz_class epsilon;
//         for(int i = 0; i < tau; i++){
//             epsilon = rand.get_z_range(bound);
//             x = chi[i] - delta[i];
//             pk_sum += x * epsilon;
//         }

//         return;
//     }


//     mpz_class encrypt(mpz_class m){
//         mpz_class r;
//         mpz_class ciphertext_mod{0};
//         mpz_class ciphertext;

//         r = rand.get_z_bits(noise);
//         computePKsum();
//         ciphertext = pk_sum + 2 * r + m;
//         mpz_mod(ciphertext_mod.get_mpz_t(), ciphertext.get_mpz_t(), x_zero.get_mpz_t());
//         return ciphertext_mod;

//     }

//     mpz_class decrypt(mpz_class c){
//         mpz_class plaintext;
//         mpz_class tmp;
//         mpz_class modulus {2};
//         mpz_mod(plaintext.get_mpz_t(), c.get_mpz_t(), prime.get_mpz_t());
//         mpz_mod(tmp.get_mpz_t(), plaintext.get_mpz_t(), modulus.get_mpz_t());
//         // gmp_printf("num: %Zd\n", tmp.get_mpz_t());
//         return tmp;

//     }

//     mpz_class add(mpz_class c1,mpz_class c2){
//         mpz_class reduced_c;
//         mpz_class addition = c1 + c2;
//         mpz_mod(reduced_c.get_mpz_t(),addition.get_mpz_t(),x_zero.get_mpz_t());
//         return reduced_c;
//     }

//     mpz_class mutliply(mpz_class c1, mpz_class c2){
//         mpz_class reduced_c;
//         mpz_class multiplication = c1 * c2;
//         mpz_mod(reduced_c.get_mpz_t(),multiplication.get_mpz_t(),x_zero.get_mpz_t());
//         return reduced_c;
//     }

//     void checkOperations(){
//         mpz_class result;
//         for(int i=0;i<3;i++){
//             mpz_class m1 = rand.get_z_range(mpz_class(2));
//             mpz_class m2 = rand.get_z_range(mpz_class(2));

//             mpz_class addition = decrypt(add(encrypt(m1),encrypt(m2)));
//             mpz_class a = m1+m2;
//             mpz_mod(result.get_mpz_t(),a.get_mpz_t(),mpz_class(2).get_mpz_t());
            
//             assert(addition == result);
//         }

//         for(int i=0;i<3;i++){
//             mpz_class m1 = rand.get_z_range(mpz_class(2));
//             mpz_class m2 = rand.get_z_range(mpz_class(2));

//             mpz_class mult = decrypt(mutliply(encrypt(m1),encrypt(m2)));
//             mpz_class a = m1*m2;
//             mpz_mod(result.get_mpz_t(),a.get_mpz_t(),mpz_class(2).get_mpz_t());
            
//             assert(mult == result);
//         }
//         return;
//     }

//     void decreasingModuli(){
//         mpz_class mu = 56;
//         ladder[ladder.size() - 1] = (levels + 1) * mu;
//         mpz_class ladder[0] = 2 * mu;
//         for (int i = levels-1; i > 1; i--) {
//             mpz_class eta = (i + 1) * mu;
//             ladder[i-1] = eta;
//         }
//     }

        
//     mpz_class getPrime() const { return prime; }
//     mpz_class getXZero() const { return x_zero; }
//     time_t getSeed() const {printf("seed: %s\n", ctime(&rand_seed));return rand_seed; }
//     std::vector<mpz_class> getDelta() const { return delta; }

//     void printPK(){
//         printf("seed: %s\n", ctime(&rand_seed));
//         char* str = new char[mpz_sizeinbase(x_zero.get_mpz_t(), 10) + 2];
//         mpz_get_str(str, 10, x_zero.get_mpz_t());
//         printf("x0: %s\n",str);
//         for(mpz_class v: delta){
//         char* str = new char[mpz_sizeinbase(v.get_mpz_t(), 10) + 2];
//         mpz_get_str(str, 10, v.get_mpz_t());
//         printf("delta: %s\n",str);
//         }
//         delete[] str;
//         return;
//     }
// };

// // big theta = 195
// // small theta = 15
// // mu = 56
// int main(){
//     for(int i = 0 ; i < 10; i++){
//         PKgenerator test(40);
//         // mpz_class c = test.encrypt(mpz_class(0));
//         // test.decrypt(c);
//         // test.checkOperations();


//     }
    
//     printf("ok\n");
    
//     return 0; 
// }