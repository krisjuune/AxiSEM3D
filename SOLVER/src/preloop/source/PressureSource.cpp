// PressureSource.cpp
// created by Claudia on 6-Apr-2018
// axial pressure source

#include "PressureSource.h"
#include "Quad.h"
#include "SpectralConstants.h"
#include "PreloopFFTW.h"
#include <sstream>

#include "Relabelling.h"

PressureSource::PressureSource(double depth, double lat, double lon, double M0): Source(depth, lat, lon),
    mM0(M0) {}

void PressureSource::computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
    arPP_CMatX3 &fouriers) const {
        
    RDMatXN K = myQuad.getRho();
    int Nr = K.rows();
    
    RDMatXN p = RDMatXN::Zero(Nr,nPntElem);

    // particle relabelling
    RDColP J_PRT;
    if (myQuad.hasRelabelling()) {
        const RDMatXN &JJ = myQuad.getRelabelling().getStiffJacobian();
        J_PRT = JJ.block(0, 0, 1, nPntEdge).transpose();
    }
    // compute source pointwise
    int ipol_src = 0; // purely on-axis
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            // spatial delta function
            RDMatPP w = RDMatPP::Zero();
            w(ipol, jpol) = 1.;
            for (int jpol_src = 0; jpol_src <= nPol; jpol_src++) {
                double fact = interpFactZ(jpol_src);
                if (myQuad.hasRelabelling()) {
                    p.col(ipnt) += RDColX::Constant(Nr, 1, w(ipol_src, jpol_src) * fact * mM0 / (2. * pi) * J_PRT[jpol_src]);
                } else {
                    p.col(ipnt) += RDColX::Constant(Nr, 1, w(ipol_src, jpol_src) * fact * mM0 / (2. * pi));
                }
            }
        }
    }
    
    p = p.schur(K);
    int Nu = int(Nr/2) + 1;
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        fouriers[ipnt] = CMatX3::Zero(Nu, 3);
        RDColX &R2C_R = PreloopFFTW::getR2C_RMat(Nr);
        R2C_R = p.col(ipnt);
        PreloopFFTW::computeR2C(Nr);
        CDColX F = PreloopFFTW::getR2C_CMat(Nr);
        fouriers[ipnt].col(0) = F.topRows(Nu);
    }
    
}

std::string PressureSource::verbose() const {
    std::stringstream ss;
    ss << "\n========================== Source ==========================" << std::endl;
    ss << "  Type            =   " << "PressureSource" << std::endl;
    ss << "  Latitude        =   " << mLatitude << std::endl;
    ss << "  Longitude       =   " << mLongitude << std::endl;
    ss << "  Depth (km)      =   " << mDepth / 1e3 << std::endl;
    ss << "  Pressure (N/m2) =   " << (mM0 >= 0. ? " " : "") << mM0 << std::endl;
    ss << "========================== Source ==========================\n" << std::endl;
    return ss.str();
}
