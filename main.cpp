#include <iostream>
#include "public_key.h"
#include "modulo_switch.h"
int main(){
    
    // PKgenerator test(50);
    // test.printPK();
    
    SwitchKey test(10,156);
    test.switchKeyGen(10);
    mpz_class c = test.enc(mpz_class(0));
    mpz_class m = test.dec(c);
    
    

    return 0;

} 



