// SolverFFTW_N.cpp
// created by Claudia on 22-May-2018
// perform FFT using fftw

#include "SolverFFTW_N.h"

int SolverFFTW_N::sNmax = 0;
std::vector<PlanFFTW> SolverFFTW_N::sR2CPlans;
std::vector<PlanFFTW> SolverFFTW_N::sC2RPlans;
RMatXN SolverFFTW_N::sR2C_RMat;
CMatXN SolverFFTW_N::sR2C_CMat;
RMatXN SolverFFTW_N::sC2R_RMat;
CMatXN SolverFFTW_N::sC2R_CMat;

void SolverFFTW_N::initialize(int Nmax) {
    int ndim = 1;
    int xx = nPntElem * ndim;
    sNmax = Nmax;
    sR2CPlans.reserve(Nmax);
    sC2RPlans.reserve(Nmax);
    sR2C_RMat = RMatXN(Nmax, xx);
    sR2C_CMat = CMatXN(Nmax / 2 + 1, xx);
    sC2R_RMat = RMatXN(Nmax, xx);
    sC2R_CMat = CMatXN(Nmax / 2 + 1, xx);
    for (int NR = 1; NR <= Nmax; NR++) {
        int NC = NR / 2 + 1;
        int n[] = {NR};
        Real *r2c_r = &(sR2C_RMat(0, 0));
        Complex *r2c_c = &(sR2C_CMat(0, 0));
        sR2CPlans.push_back(planR2CFFTW(1, n, xx, r2c_r, n, 1, Nmax, complexFFTW(r2c_c), n, 1, Nmax / 2 + 1, SolverFFTW::mWisdomLearnOption));   
        Real *c2r_r = &(sC2R_RMat(0, 0));
        Complex *c2r_c = &(sC2R_CMat(0, 0));
        sC2RPlans.push_back(planC2RFFTW(1, n, xx, complexFFTW(c2r_c), n, 1, Nmax / 2 + 1, c2r_r, n, 1, Nmax, SolverFFTW::mWisdomLearnOption)); 
    }
}

void SolverFFTW_N::finalize() {
    for (int i = 0; i < sNmax; i++) {
        distroyFFTW(sR2CPlans[i]);
        distroyFFTW(sC2RPlans[i]);
    }
    sNmax = 0;
}

void SolverFFTW_N::computeR2C(int nr) {
    execFFTW(sR2CPlans[nr - 1]);
    Real inv_nr = one / (Real)nr;
    sR2C_CMat.topRows(nr / 2 + 1) *= inv_nr;
}

void SolverFFTW_N::computeC2R(int nr) {
    execFFTW(sC2RPlans[nr - 1]);
}
