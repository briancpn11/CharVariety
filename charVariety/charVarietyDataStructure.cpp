// charVarietyDataStructure.cpp

#include "charVarietyDataStructure.h"

int main(){
    double eps = 1e-6;
    int param = 3;
    std::cout.precision(10);
    Interval x(param-eps, param+eps);
    std::cout << sin(x) << std::endl;
    return 0;
}