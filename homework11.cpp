#include <iostream>
int Ctepen(int a, int b){
    int answer = 1;
    for(int k=0; k<b; k++){
        answer *= a;
        }
    return answer;
}
using namespace std;

float Electric_Capacitance(float x){
    const float c = 3*Ctepen(10,8);
    x *= Ctepen(10,9)/(c*c);
    return(x);
}

float Magnetic_Field_Strength(float x){
    float pi = 3.1415;
    x *= Ctepen(10,3)/(4*pi);
    return(x);
}
int main(){
int y;
cin >> y;
cout << Magnetic_Field_Strength(y);

    return(0);
}
