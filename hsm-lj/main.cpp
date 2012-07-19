    //
    //  main.cpp
    //  hsm
    //
    //  Created by Ionuț Georgescu on 6/9/12.
    //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
    //

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>

    //#include <vecLib/vecLib.h>

extern "C" {
    double dlamch_(char *);
    int dsyevr_(char *jobz, char *range, char *uplo, int *n, 
                double *a, int *lda, double *vl, double *vu, int *
                il, int *iu, double *abstol, int *m, double *w, 
                double *z__, int *ldz, int *isuppz, double *work, 
                int *lwork, int *iwork, int *liwork, int *info);
}

int Natom;
double  *x, *y, *z;
double U0, *H, *W;

char linebuf[256];

using namespace std;

static void hessian(double *x, double *y, double *z, double &U0, double *H)
{
    int i, j, k, l;
    double dx, dy, dz, A, B;
    double r2_, r22_, ir24_, ir28_, x_, y_, z_, H0[9];
    double *Hdiag = new double[9*Natom];
    
    U0 = 0.0;
    fill_n(Hdiag, 9*Natom, 0.0);
    for (i=0; i<Natom; i++) {
        for (j=i+1; j<Natom; j++) {
            dx = x[j] - x[i];
            dy = y[j] - y[i];
            dz = z[j] - z[i];
            
            r2_ = dx*dx + dy*dy + dz*dz;
            r22_ = r2_ * r2_;
            ir24_ = 1.0/(r22_ * r22_);
            ir28_ = ir24_*ir24_;
            
            U0 += (ir28_ - ir24_)*r2_;
            
            A = 672.0*ir28_ - 192.0*ir24_/r2_;
            B = 24.0*ir24_ - 48.0 * ir28_ * r2_;
            
            x_ = A * dx;
            y_ = A * dy;
            z_ = A * dz;
            
            H0[0] = x_*dx + B;
            H0[1] = x_ * dy;
            H0[2] = x_ * dz;
            
            H0[3] = y_*dx;
            H0[4] = y_ * dy + B;
            H0[5] = y_ * dz;
            
            H0[6] = z_*dx;
            H0[7] = z_ * dy;
            H0[8] = z_ * dz + B;
            
            for (k=0; k<9; k++) {
                Hdiag[9*i + k] += H0[k];
            }
            for (k=0; k<9; k++) {
                Hdiag[9*j + k] += H0[k];
            }
            
            for (k=0; k<3; k++) {
                for (l=0; l<3; l++) {
                    H[9*Natom*i + 3*j + 3*Natom*k+l] = -H0[3*k+l];
                }
            }
        }
        for (k=0; k<3; k++) {
            for (l=0; l<3; l++) {
                H[9*Natom*i + 3*i+3*Natom*k+l] = Hdiag[9*i + 3*k+l];
            }
        }
    }
    U0 *= 4.0;
    
    delete[] Hdiag;
}


int main(int argc, const char * argv[])
{
    int i;
    stringstream ss;
    string fname;
    
    if (argc==2) {
        ifstream fin(argv[1]);
        
        fin >> Natom;
        
        x = new double[Natom];
        y = new double[Natom];
        z = new double[Natom];
        H = new double[9*Natom*Natom];
        
        fin.getline(linebuf, sizeof(linebuf));
        fin.getline(linebuf, sizeof(linebuf));
        
        for (i=0; i<Natom; i++) {
            fin >> linebuf >> x[i] >> y[i] >> z[i];
        }
        fin.close();
        
        fname=argv[1];
    }
    else {
        ss << argv[1];
        ss >> Natom;
        
        x = new double[Natom];
        y = new double[Natom];
        z = new double[Natom];
        H = new double[9*Natom*Natom];
        
        ifstream fin(argv[2]);
        for (i=0; i<Natom; i++) {
            fin >> x[i] >> y[i] >> z[i];
        }
        fin.close();
        
        fname=argv[2];
    }
    
    hessian(x, y, z, U0, H);
    
    char jobz, range, uplo, cmach;
    int N, M, lwork, liwork, info, *iwork, *isuppz;
    double *work, *Z, abstol;
    double vl=numeric_limits<double>::min(),
    vu=numeric_limits<double>::max();
    int il=1, iu=3*Natom;
    
    jobz='N';
    range='A';
    uplo='L';
    N = 3*Natom;
    lwork = 26*N;
    liwork = 10*N;
    W = new double[N];
    Z = new double[N];
    work = new double[lwork];
    iwork = new int[liwork];
    isuppz = new int[2*N];
    cmach = 'S';
    abstol = dlamch_(&cmach);
    
    dsyevr_(&jobz, &range, &uplo, &N, H, &N, &vl, &vu, &il, &iu, &abstol, &M, W, 0, &N, isuppz, work, &lwork, iwork, &liwork, &info);
    
    string ofname = fname.substr(0, fname.length()-4) + "_spectr.dat";
    ofstream fout(ofname.c_str());
    for (i=0; i<3*Natom; i++) fout << W[i] << endl;
    fout.close();
}

