#include <iostream>
#include "public_key.h"
#include "modulo_switch.h"
#include <cstdlib> // for srand()
#include <ctime> // for time()
#include <thread> // for sleep_for
#include <gperftools/profiler.h>

int main(){

    for(int i = 0;i<1 ;i++){
        int levels = 5;
        int theta = 735;
        int rho = 20;
        int gamma = 270000;
        int mu = 65;
        SwitchKey test(levels,theta,rho,gamma,mu);
        test.switchKeyGen(levels);
        mpz_class c1 = test.enc(mpz_class(1));
        mpz_class c2 = test.enc(mpz_class(0));

        mpz_class m = test.add(c1,c2);
        std::cout<< "\n";
        // std::this_thread::sleep_for(std::chrono::seconds(1)); 

    }


      // PKgenerator test(1026,158,27,150000);
      // mpz_class c1 = test.encrypt(mpz_class(1));
      // mpz_class d = test.decrypt(c1);
    //   for(int i = 0;i<1 ;i++){

    //     PKgenerator test(1026,158,27,150000);
    //     mpz_class c1 = test.encrypt(mpz_class(1));
    //     mpz_class c2 = test.encrypt(mpz_class(1));

    //     mpz_class m = test.add(c1,c2);
    //     mpz_class d = test.decrypt(m);

    //     std::cout<< "\n";
    //     std::this_thread::sleep_for(std::chrono::seconds(1)); // wait for 3 seconds

    // }


    return 0;

} 
