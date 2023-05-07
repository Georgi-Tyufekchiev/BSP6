#include <iostream>
#include "public_key.h"
#include "modulo_switch.h"
#include <cstdlib> // for srand()
#include <ctime> // for time()
#include <thread> // for sleep_for
#include <gperftools/profiler.h>

int main(){

    // for(int i = 0;i<10 ;i++){

    //     SwitchKey test(10,156);
    //     test.switchKeyGen(10);
    //     mpz_class c1 = test.enc(mpz_class(1));
    //     mpz_class c2 = test.enc(mpz_class(0));

    //     mpz_class m = test.add(c1,c2);
    //     std::cout<< "\n";
    //     // std::this_thread::sleep_for(std::chrono::seconds(1)); 

    // }


    //   for(int i = 0;i<10 ;i++){

    //     PKgenerator test(1000,156);
    //     mpz_class c1 = test.encrypt(mpz_class(1));
    //     mpz_class c2 = test.encrypt(mpz_class(1));

    //     mpz_class m = test.add(c1,c2);
    //     mpz_class d = test.decrypt(m);

    //     std::cout<< "\n";
    //     std::this_thread::sleep_for(std::chrono::seconds(1)); // wait for 3 seconds

    // }
    // ProfilerStart("cpu_profile.prof"); // Start CPU profiling

    std::srand(1234);
    int cnt = 0;
    for(int i = 0;i<1;i++){
        bool res = attackGACD();
        if(res){
            cnt+=1;
        }
        // printf("\n");
    }
     std::cout << "Success:" <<cnt << std::endl;

    

    return 0;

} 
