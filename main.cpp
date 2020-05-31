#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>

#include "SIMDMatrixes.h"

using namespace std;


int main() {
    while (true) {
        srand(time(nullptr));
        int i, j, k;

        int timeOfSumm = 1;
        int timeOfDifference = 1;
        int timeOfMultiplying = 1;
        int timeOfCom = 1;

        SIMDMatrixes lab2;

        cout << "Input m, p, q, n" << endl;
        cin >> lab2.m >> lab2.p >> lab2.q >> lab2.n;

        lab2.r = lab2.p * lab2.m * lab2.q;
        lab2.Tn = 0;
        lab2.callsOfSumm = 0;
        lab2.callsOfDifference = 0;
        lab2.callsOfMultiplying = 0;
        lab2.callsOfCom = 0;
        lab2.Lavg = 0;

        lab2.A = new double* [lab2.p];
        lab2.B = new double* [lab2.m];
        lab2.E = new double[lab2.m];
        lab2.G = new double* [lab2.p];
        lab2.C = new double* [lab2.p];


        for (i = 0; i < lab2.p; i++) {
            lab2.A[i] = new double[lab2.m];
            lab2.G[i] = new double[lab2.q];
            lab2.C[i] = new double[lab2.q];
        }

        for (i = 0; i < lab2.m; i++) {
            lab2.B[i] = new double[lab2.q];
        }

        for (i = 0; i < lab2.p; i++) {
            for (j = 0; j < lab2.m; j++) {
                lab2.A[i][j] = (double)(rand() % 20001) / 10000 - 1;
            }
        }
        for (i = 0; i < lab2.m; i++) {
            for (j = 0; j < lab2.q; j++) {
                lab2.B[i][j] = (double)(rand() % 20001) / 10000 - 1;
            }
        }
        for (i = 0; i < lab2.m; i++) {
            lab2.E[i] = (double)(rand() % 20001) / 10000 - 1;
        }
        for (i = 0; i < lab2.p; i++) {
            for (j = 0; j < lab2.q; j++) {
                lab2.G[i][j] = (double)(rand() % 20001) / 10000 - 1;
            }
        }


        lab2.F = new double** [lab2.p];
        for (i = 0; i < lab2.p; i++)
            lab2.F[i] = new double* [lab2.q];
        for (i = 0; i < lab2.p; i++)
            for (j = 0; j < lab2.q; j++)
                lab2.F[i][j] = new double[lab2.m];

        lab2.D = new double** [lab2.p];
        for (i = 0; i < lab2.p; i++)
            lab2.D[i] = new double* [lab2.q];
        for (i = 0; i < lab2.p; i++)
            for (j = 0; j < lab2.q; j++)
                lab2.D[i][j] = new double[lab2.m];

        cout << endl << "A:";
        for (i = 0; i < lab2.p; i++) {
            cout << endl;
            for (j = 0; j < lab2.m; j++) {
                cout << setw(7) << lab2.A[i][j] << " ";
            }
        }
        cout << endl << endl;


        cout << "B:";
        for (i = 0; i < lab2.m; i++) {
            cout << endl;
            for (j = 0; j < lab2.q; j++) {
                cout << setw(7) << lab2.B[i][j] << " ";
            }
        }
        cout << endl << endl;


        cout << "E:";
        cout << endl;
        for (i = 0; i < lab2.m; i++) {
            cout << setw(7) << lab2.E[i] << " ";
        }
        cout << endl << endl;


        cout << "G:";
        for (i = 0; i < lab2.p; i++) {
            cout << endl;
            for (j = 0; j < lab2.q; j++) {
                cout << setw(7) << lab2.G[i][j] << " ";
            }
        }

        for (i = 0; i < lab2.p; i++)
            for (j = 0; j < lab2.q; j++)
                for (k = 0; k < lab2.m; k++)
                    lab2.Fijk(lab2.A, lab2.B, lab2.E, i, j, k);
        int reductionTime = 3 * (timeOfMultiplying + timeOfDifference + timeOfSumm);
        int operationTime = 7 * timeOfMultiplying + 3 * timeOfDifference + 2 * timeOfSumm;
        lab2.Tn += (int)(reductionTime + operationTime) * (int)ceil((double)lab2.r / lab2.n);

        for (i = 0; i < lab2.p; i++)
            for (j = 0; j < lab2.q; j++)
                for (k = 0; k < lab2.m; k++)
                    lab2.Dijk(lab2.A, lab2.B, i, j, k);

        lab2.Tn += timeOfMultiplying * (int)ceil((double)lab2.r / lab2.n);


        for (i = 0; i < lab2.p; i++) {
            for (j = 0; j < lab2.q; j++) {
                lab2.C[i][j] = lab2.cij(i, j);
            }


        }


        int FTime = 2 * timeOfMultiplying * (lab2.m - 1);
        int DTime = 3 * (timeOfDifference * (lab2.m + 1) + timeOfMultiplying * (lab2.m - 1));
        operationTime = 8 * timeOfMultiplying + 4 * timeOfDifference + 2 * timeOfSumm;
        lab2.Tn += (FTime + DTime + operationTime) * (int)ceil((double)(lab2.p * lab2.q) / (double)lab2.n);


        cout << endl << endl;
        cout << "C:";
        for (i = 0; i < lab2.p; i++) {
            cout << endl;
            for (j = 0; j < lab2.q; j++) {
                cout << setw(12) << lab2.C[i][j] << " ";
            }
        }
        cout << endl << endl;

        cout << "----------------" << endl << "Parameters:" << endl << "----------------" << endl;
        lab2.T1 = timeOfSumm * lab2.callsOfSumm + timeOfDifference * lab2.callsOfDifference + timeOfMultiplying * lab2.callsOfMultiplying + timeOfCom * lab2.callsOfCom;
        if (lab2.Tn > lab2.T1) lab2.Tn = lab2.T1;
        lab2.Ky = (double)lab2.T1 / lab2.Tn;
        lab2.Eff = (double)lab2.Ky / lab2.n;
        //D
        lab2.Lavg = timeOfMultiplying * lab2.r;
        //F
        lab2.Lavg += (7 * timeOfMultiplying + 3 * timeOfDifference + 2 * timeOfSumm) * lab2.r;
        //C
        lab2.Lavg += (8 * timeOfMultiplying + 3 * timeOfDifference + 2 * timeOfSumm) * lab2.r;
        

        //F_func
        lab2.Lavg += (timeOfMultiplying * (lab2.m - 1)) * 2 * lab2.r;
        //D_func
        lab2.Lavg += (timeOfMultiplying * (lab2.m - 1) + timeOfDifference * (lab2.m + 1)) * 3 * lab2.r;
        
        //FandD
        lab2.Lavg += (timeOfCom + timeOfSumm + timeOfDifference) * (lab2.m - 1) * lab2.r;
        
        //AtoB
        lab2.Lavg += (timeOfCom + timeOfDifference) * (lab2.m - 1) * 2 * lab2.r;

        //BtoA
        lab2.Lavg += (timeOfCom + timeOfDifference) * (lab2.m - 1) * 2 * lab2.r;
       
        //AandB
        lab2.Lavg += (timeOfDifference) * (lab2.m - 1) * 2 * lab2.r;
         

        lab2.Lavg = ceil(lab2.Lavg / lab2.r);
        lab2.Diff = (double)lab2.Tn / lab2.Lavg;

        cout << "T1   = " << lab2.T1 << "\nTn   = " << lab2.Tn << "\nKy   = " << lab2.Ky << "\ne    = " <<
            lab2.Eff << "\nLsum = " << lab2.Tn << "\nLavg = " << lab2.Lavg << "\nD    = " << lab2.Diff << endl;
        cout << "----------------" << endl;

        cout << endl << endl << endl << endl << endl << endl;

    }
    return 0;
}

