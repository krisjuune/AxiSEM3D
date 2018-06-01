// PressureSource.cpp
// created by Claudia on 6-Apr-2018
// axial pressure source

#include "PressureSource.h"
#include "Quad.h"
#include "SpectralConstants.h"
#include <sstream>

#include "Relabelling.h"

PressureSource::PressureSource(double depth, double lat, double lon, double M0): Source(depth, lat, lon),
    mM0(M0) {}

void PressureSource::computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
    arPP_CMatX3 &fouriers) const {
    // set zero
    for (int ipnt = 0; ipnt < nPntElem; ipnt++) {
        fouriers[ipnt] = CMatX3::Zero(3, 3);
    }
    // Jacobian on axis
    std::array<RDMat22, nPntEdge> axJ;
    int ipol_src = 0;
    for (int jpol_src = 0; jpol_src <= nPol; jpol_src++) {
        const RDCol2 &xieta = SpectralConstants::getXiEta(ipol_src, jpol_src, true);
        axJ[jpol_src] = myQuad.jacobian(xieta);
        axJ[jpol_src] /= axJ[jpol_src].determinant();
    }
    // particle relabelling
    RDColP VX0, VX1, VX2, VX3;
    if (myQuad.hasRelabelling()) {
        const RDMatXN4 &X = myQuad.getRelabelling().getStiffX();
        VX0 = X.block(0, nPE * 0, 1, nPntEdge).transpose();
        VX3 = X.block(0, nPE * 3, 1, nPntEdge).transpose();
        // VX1 and VX2 are phi-independent. The following lines assume phi = 0
        VX1 = X.block(0, nPE * 1, 1, nPntEdge).transpose();
        VX2 = X.block(0, nPE * 2, 1, nPntEdge).transpose();
    }
    // compute source pointwise
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            // spatial delta function
            RDMatPP w = RDMatPP::Zero();
            w(ipol, jpol) = 1.;
            const RDMatPP &GU = SpectralConstants::getG_GLJ().transpose() * w;
            const RDMatPP &UG = w * SpectralConstants::getG_GLL();
            for (int jpol_src = 0; jpol_src <= nPol; jpol_src++) {
                double fact = interpFactZ(jpol_src);
                const RDMat22 &J = axJ[jpol_src];
                double dwds = J(1, 1) * GU(ipol_src, jpol_src) - J(1, 0) * UG(ipol_src, jpol_src);
                double dwdz = J(0, 0) * UG(ipol_src, jpol_src) - J(0, 1) * GU(ipol_src, jpol_src);
                if (myQuad.hasRelabelling() && !myQuad.isFluid()) {
                    double X0 = VX0(jpol_src);
                    double X1 = VX1(jpol_src);
                    double X2 = VX2(jpol_src);
                    double X3 = VX3(jpol_src);
                    // monopole
                    fouriers[ipnt](0, 0) += fact * dwds * 2 * mM0 * X0 / (2. * pi);
                    fouriers[ipnt](0, 2) += fact * dwdz * mM0 * X3 / (2. * pi);
                    // dipole
                    fouriers[ipnt](1, 0) += fact * dwdz * mM0 * (X1 - X2) / (4. * pi);
                    fouriers[ipnt](1, 1) += fact * dwdz * mM0 * (X1 - iid * X2) / (4. * pi) * iid;
                } else {
                    if (myQuad.hasRelabelling() && myQuad.isFluid()) {
                        std::cout << "Warning: relabelling of fluid source neglected." << std::endl;
                    }
                    // monopole
                    fouriers[ipnt](0, 0) += fact * dwds * 2 * mM0 / (2. * pi);
                    fouriers[ipnt](0, 2) += fact * dwdz * mM0 / (2. * pi);
                }
            }
        }
    }
}

std::string PressureSource::verbose() const {
    std::stringstream ss;
    ss << "\n========================== Source ==========================" << std::endl;
    ss << "  Type         =   " << "PressureSource" << std::endl;
    ss << "  Latitude     =   " << mLatitude << std::endl;
    ss << "  Longitude    =   " << mLongitude << std::endl;
    ss << "  Depth (km)   =   " << mDepth / 1e3 << std::endl;
    ss << "  Moment (N.m) =   " << (mM0 >= 0. ? " " : "") << mM0 << std::endl;
    ss << "========================== Source ==========================\n" << std::endl;
    return ss.str();
}
