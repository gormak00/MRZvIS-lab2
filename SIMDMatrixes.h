#pragma once

using namespace std;

class SIMDMatrixes {
    public:
        int m;
        int p;
        int q;
        int n;
        int r;

        double** A;
        double** B;
        double* E;
        double** G;
        double*** F;
        double*** D;
        double** C;


        int callsOfSumm;
        int callsOfDifference;
        int callsOfMultiplying;
        int callsOfCom;

        int T1;
        int Tn;
        double Ky;
        double Eff;
        double Diff;
        int Lavg;
    

        double cij(int i, int j);

        void Fijk(double** A, double** B, double* E, int i, int j, int k);

        void Dijk(double** A, double** B, int i, int j, int k);



    private:

        double a_and_b(double** A, double** B, int i, int k, int j);

        double a_to_b(double** A, double** B, int i, int k, int j);

        double b_to_a(double** A, double** B, int i, int k, int j);

        double f_func(int i, int j);

        double d_func(int i, int j);

        double f_and_d(int i, int j);

};