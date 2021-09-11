//charVariety/utils/helpers.cpp

#include "helpers.h"

#define MAXITER 30 // max number of iterations in a Newton scheme

// Actions at the tangent space level
Pair DTX(Pair pr){
    Pair prpr(pr.x, pr.z, (pr.x)*(pr.z) - (pr.y), 
                pr.v1, pr.v3, (pr.z)*(pr.v1) - pr.v2 + (pr.x)*(pr.v3));
    return prpr;
}

Pair DTY(Pair pr){
    Pair prpr(pr.z, pr.y, (pr.y)*(pr.z) - (pr.x), 
                pr.v3, pr.v2, -pr.v1 + (pr.z)*(pr.v2) + (pr.y)*(pr.v3));
    return prpr;
}

Pair DTX1(Pair pr){
    Pair prpr(pr.x, (pr.x)*(pr.y) - (pr.z), pr.y, 
                pr.v1, (pr.y)*(pr.v1) + (pr.x)*(pr.v2) - pr.v3, pr.v2);
    return prpr;
}

Pair DTY1(Pair pr){
    Pair prpr((pr.x)*(pr.y) - (pr.z), pr.y, pr.x, 
                (pr.y)*(pr.v1) + (pr.x)*(pr.v2) - pr.v3, pr.v2, pr.v1);
    return prpr;
}

Pair action(Pair aa, int* maps, int nIter){
    Pair a = aa;
    for (int j = nIter - 1; j >= 0; j--){
        switch (maps[j]){
            case 0: a = DTX(a); break;
            case 1: a = DTY(a); break;
            case 2: a = DTX1(a); break;
            case 3: a = DTY1(a); break;
        }
    }
    return a;
}

// Calculate the average log expansion of a list of matrices acting on a particular direction.
double avgLogExpansion(Matrix* A, int nMatrix, double angle){
    double ans = 0.0;
    for (int i = 0; i < nMatrix; i++){
        ans += 0.5 * log(A[i].expansionSq(angle));
    }
    return ans / double(nMatrix);
}

double minAvgExpansionGrid(Matrix* A, int nMatrix, double rr){
    double minAvgExpansion = 100000.0;
    for (double theta = 0; theta < PI; theta += rr){
        minAvgExpansion = std::min(minAvgExpansion, avgLogExpansion(A, nMatrix, theta));
    }
    return minAvgExpansion;
}

// Find the minimum of average expansion of a list of matrices
double minAvgExpansionNewton(Matrix* A, int nMatrix){

    double minAvgExpansion = 100000.0;
    double angleNew, angleOld;
    int co;

    double expSq, dExpSq, ddExpSq;
    double dTotalExpSq = 0.0, ddTotalExpSq = 0.0;

    for (int j = 0; j < nMatrix; j++){
        angleNew = A[j].contractDir();
        angleOld = angleNew - 1;
        co = 0;
        while (absVal(angleNew - angleOld) > EPS && co < MAXITER){
            angleOld = angleNew;
            dTotalExpSq = 0.0; 
            ddTotalExpSq = 0.0;
                        
            for (int i = 0; i < nMatrix; i++){
                expSq = A[i].coeff0() + A[i].coeff1() * cos(2 * angleOld) + A[i].coeff2() * sin(2 * angleOld);
                dExpSq = 2 * (-A[i].coeff1() * sin(2 * angleOld) + A[i].coeff2() * cos(2 * angleOld));
                ddExpSq = 4 * (-A[i].coeff1() * cos(2 * angleOld) - A[i].coeff2() * sin(2 * angleOld));
                
                dTotalExpSq += 0.5 * dExpSq / expSq;
                ddTotalExpSq += 0.5 * (ddExpSq / expSq - sqr(dExpSq / expSq));
            }
        
            angleNew = angleOld - dTotalExpSq / ddTotalExpSq;
            co++;
        }
        // cMax = (cMax < co ? co:cMax);
        // if (co == coMax) bad_coMax++;
        minAvgExpansion = std::min(minAvgExpansion, avgLogExpansion(A, nMatrix, angleNew));
    }
    return minAvgExpansion;
}




// int main(){
//     double eps = 1e-6;
//     int param = 3;
//     Interval x(param-eps, param+eps);
// 	Pair pp(1, 2, 3, 4, 5, 6);
// 	Pair qq(1, 2, 3, 4, 5, 6);
// 	Pair p = DTY1(pp);
// 	Matrix m(0, -1, 1, 0), m2(1, 2, 3, 4);
// 	RotationMatrix rho(PI/2);

// 	std::cout.precision(10);
//     // std::cout << "sin(x) = " << sin(x) << std::endl;
// 	std::cout << (pp == qq) << std::endl;
// 	pp.printPair();
// 	p.printPair();
// 	m.printMatrix();
// 	(m * m2).printMatrix();
// 	rho.printMatrix();
//     return 0;
// }