#include <iostream>
#include "public_key.h"
#include "modulo_switch.h"
int main(){
    
    // PKgenerator test(50);
    // test.printPK();
    SwitchKey test(10,156);
    test.switchKeyGen(10);
    

    return 0;

} 



