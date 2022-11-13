#include <iostream>

using namespace std;


int *swap(int a, int b){
    int base[2];
    a = a - b;
    b = a + b;
    a = b - a;
    base[0]=a;
    base[1]=b;
    return base;
}

  
