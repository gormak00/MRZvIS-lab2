#include "SIMDMatrixes.h"

using namespace std;

//deafult task c(ij)
double SIMDMatrixes::cij(int i, int j) {
    double c;
    callsOfMultiplying += 8;
    callsOfDifference += 3;
    callsOfSumm += 2;
    c = f_func(i, j) * (3. * G[i][j] - 2.) * G[i][j] + (d_func(i, j) + (4. * f_and_d(i, j) - 3. * d_func(i, j)) * G[i][j]) * (1. - G[i][j]);

    return c;
}

//deafult task f(ijk)
void SIMDMatrixes::Fijk(double** A, double** B, double* E, int i, int j, int k) {
    callsOfMultiplying += 7;
    callsOfDifference += 3;
    callsOfSumm += 2;
    F[i][j][k] = (a_to_b(A, B, i, k, j) * (2. * E[k] - 1.) * E[k] + b_to_a(A, B, i, k, j) * (1. + (4. * a_to_b(A, B, i, k, j) - 2.) * E[k]) * (1. - E[k]));
}

//deafult task d(ijk)
void SIMDMatrixes::Dijk(double** A, double** B, int i, int j, int k) {
    D[i][j][k] = a_and_b(A, B, i, k, j);
}






//F func
double SIMDMatrixes::f_func(int i, int j) {
    double result = 1;
    for (int k = 0; k < m; k++) {
        result *= F[i][j][k];
    }
    callsOfMultiplying += m - 1; 
    return result;
}


//D func
double SIMDMatrixes::d_func(int i, int j) {
    double result = 1;
    for (int k = 0; k < m; k++) {
        result *= 1 - D[i][j][k];
    }
    callsOfDifference += m + 1;
    callsOfMultiplying += m - 1;
    callsOfDifference++;
    return 1 - result;
}




//F and D
double SIMDMatrixes::f_and_d(int i, int j) {
    double res, half;
    half = (f_func(i, j) + d_func(i, j) - 1.);
    if (half > 0.) res = half;
    else res = 0.;
    callsOfDifference++;
    callsOfSumm++;
    callsOfCom++;
    return res;
}


//A to B
double SIMDMatrixes::a_to_b(double** A, double** B, int i, int k, int j) {
    double res;
    if ((1 - A[i][k]) > B[k][j]) res = (1 - A[i][k]);
    else res = B[k][j];
    callsOfDifference++;
    callsOfCom++;
    return res;

}


//B to A
double SIMDMatrixes::b_to_a(double** A, double** B, int i, int k, int j) {
    double res;
    if ((1 - B[k][j]) > A[i][k]) res = (1 - B[k][j]);
    else res = A[i][k];
    callsOfDifference++;
    callsOfCom++;
    return res;
}


//A and B
double SIMDMatrixes::a_and_b(double** A, double** B, int i, int k, int j) {
    double res;
    if (A[i][k] > B[k][j]) res = A[i][k];
    else res = B[k][j];
    callsOfDifference++;
    return res;
}