#include "WavefieldRecorder.h"
#include "Parameters.h"
#include "Element.h"
#include "XMPI.h"
#include "mpi.h"
#include <fstream>

WavefieldRecorder::WavefieldRecorder(double dt, double t_ini, std::vector<double> phi, int nPhi):
mDt(dt), mT_ini(t_ini), mPhi(phi), mNphi(nPhi) {}

void WavefieldRecorder::buildInparam(WavefieldRecorder *&WFR, const Parameters &par) {
    if (WFR) {
        delete WFR;
    }

    const std::string source = "WavefieldRecorder::buildInparam";

    int n = par.getValue<int>("OUT_WAVEFIELD_FRAMES");
    double t_max = par.getValue<double>("TIME_RECORD_LENGTH");
    double t_ini = par.getValue<double>("OUT_WAVEFIELD_TIME_START");
    double dt = (t_max - t_ini) / n;

    std::string mstr = par.getValue<std::string>("OUT_WAVEFIELD_PHI");
    std::vector<std::string> strs = Parameters::splitString(mstr, "$");
    int nPhi = strs.size();
    std::vector<double> phi;
    double p;
    for (int i = 0; i < nPhi; i++) {
        Parameters::castValue(p, strs[i], source);
        phi.push_back(p);
    }

    if (n >= 1) {
        WFR = new WavefieldRecorder(dt, t_ini, phi, nPhi);
    }
}

void WavefieldRecorder::initialize() {
    mOutputPath = Parameters::sOutputDirectory + "/wavefields/";

    std::fstream t_file;
    t_file.open(mOutputPath + "times.txt", std::fstream::out);
    t_file.close();

    mFrame = 0;
    mT_record = mT_ini;
}

void WavefieldRecorder::record(Real t, const std::vector<Element *> elements) {
    if (t >= mT_record) {
        for (int iphi = 0; iphi < mNphi; iphi++) {
            RDMatXX local_buffer = RDMatXX::Zero(elements.size(), 6);
            int i = 0;
            for (const auto &elem: elements) {
                local_buffer.block(i,0,1,3) = elem->getCenterCrds(degree * mPhi[iphi]).transpose();
                local_buffer.block(i,3,1,3) = elem->recordWF(degree * mPhi[iphi]);
                i++;
            }
            std::stringstream ss;
            ss << local_buffer << std::endl;
            std::vector<std::string> global_buffer;
            XMPI::gather(ss.str(), global_buffer, false);
            if (XMPI::root()) {
                std::fstream Wout;
                Wout.open(mOutputPath + "Phi_" + std::to_string(mPhi[iphi]) + "_Frame" + std::to_string(mFrame) + ".txt", std::fstream::out);
                for (const auto &local: global_buffer) {
                    Wout << local;
                }
                Wout.close();
            }
        }
        if (XMPI::root()) {
            std::fstream recTime;
            recTime.open(mOutputPath + "times.txt", std::fstream::app);
            recTime << t << std::endl;
            recTime.close();
        }
        mT_record += mDt;
        mFrame++;
    }
}
