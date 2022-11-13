#include <iostream>
#include <cmath>


int main()
{
    float a; float b; float c;
    float x1; float x2; float x1_complex; float x2_complex;
    float D;
    std::cin >> a >> b >> c;

    D = - 4*a*c + b*b;
    if(D < 0){
        x1 = (-b)/(2*a);
        x1_complex = (sqrt(-D))/(2*a);
        x2_complex = (sqrt(-D))/(2*a);
        std::cout << "Dual complex solution:"<< '\n';
        if (x1 == 0) {
        std::cout << "-" << x1_complex<<"i" <<  '\n';
        std::cout << "+" << x1_complex<<"i" <<  '\n';
        }
        else{
        std::cout << x1 << "-" << x1_complex<<"i" <<  '\n';
        std::cout << x2 << "+" << x1_complex<<"i" <<  '\n';}
    }
    else if(D == 0){
        x1 = -b/(2*a);
        std::cout << "Single solution: " << x1;
    }
    else{
        x1 = (-b-sqrt(D))/(2*a);
        x2 = (-b+sqrt(D))/(2*a);
        std::cout << "Dual solution: " << x1 << ", " << x2;
    }
    return 0;
}