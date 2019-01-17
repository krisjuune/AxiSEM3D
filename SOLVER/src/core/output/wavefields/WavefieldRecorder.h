#pragma once
#include "Element.h"
#include "Parameters.h"

class WavefieldRecorder {

public:
    WavefieldRecorder(double dt, double t_ini, std::vector<double> phi, int nPhi);

    static void buildInparam(WavefieldRecorder *&WFR, const Parameters &par);
    void initialize();
    void record(Real t, const std::vector<Element *> elems);

private:
    double mDt, mT_ini, mT_record;

    int mNphi;
    std::vector<double> mPhi;

    int mFrame;
    std::string mOutputPath;
};
