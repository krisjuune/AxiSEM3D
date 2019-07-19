// Element.cpp
// created by Kuangdai on 27-Mar-2016 
// base class of AxiSEM3D spectral elements

#include "Element.h"
#include "Point.h"
#include "Acoustic.h"

#include "FieldFFT.h"
#include "SolverFFTW_N3.h"
#include "Gradient.h"
#include "PRT.h"
#include "XMath.h"

Element::Element(Gradient *grad, PRT *prt, const std::array<Point *, nPntElem> &points):
mGradient(grad), mPRT(prt), mHasPRT(prt != 0) {
    mMaxNr = mMaxNu = -1;
    for (int i = 0; i < nPntElem; i++) {
        mPoints[i] = points[i];
        mMaxNr = std::max(mMaxNr, mPoints[i]->getNr());
        mMaxNu = std::max(mMaxNu, mPoints[i]->getNu());
    }
    if (mHasPRT) {
        mPRT->checkCompatibility(mMaxNr);
    }
}

Element::~Element() {
    delete mGradient;
    if (mHasPRT) {
        delete mPRT;
    }
}

void Element::addSourceTerm(const arPP_CMatX3 &source) const {
    arPP_CMatX3 SF_source = source;
    if (fluid()) {
        const Acoustic *acoust = getAcoustic();
        RMatXN K = acoust->getRho();

        vec_ar3_CMatPP source_F = vec_ar3_CMatPP(mMaxNu + 1, zero_ar3_CMatPP);
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                int ipnt = ipol * nPntEdge + jpol;
                for (int alpha = 0; alpha < SF_source[ipnt].rows(); alpha++) {
                    for (int i = 0; i < 3; i++) {
                        source_F[alpha][i](ipol, jpol) = SF_source[ipnt](alpha, i);
                    }
                }
            }
        }

        RMatXN3 &source_P = SolverFFTW_N3::getR2C_RMat();
        FieldFFT::transformF2P(source_F, mMaxNr);
        source_P = SolverFFTW_N3::getC2R_RMat().topRows(mMaxNr);
        for (int alpha = 0; alpha <= mMaxNu; alpha++) {
            for (int i = 0; i < 3; i++) {
                for (int ipnt = 0; ipnt < nPE; ipnt++) {
                    source_P(alpha, i * nPE + ipnt) *= K(alpha, ipnt);
                }
            }
        }
        FieldFFT::transformP2F(source_F, mMaxNr);
        for (int i = 0; i < 3; i++) {
            for (int ipol = 0; ipol <= nPol; ipol++) {
                for (int jpol = 0; jpol <= nPol; jpol++) {
                    int ipnt = ipol * nPntEdge + jpol;
                    for (int alpha = 0; alpha < SF_source[ipnt].rows(); alpha++) {
                        SF_source[ipnt](alpha, i) = source_F[alpha][i](ipol, jpol);
                    }
                }
            }
        }
    }
    for (int i = 0; i < nPntElem; i++) {
        mPoints[i]->addToStiff(SF_source[i]);
    }
}

bool Element::axial() const {
    return mPoints[0]->axial();
}

std::string Element::costSignature() const {
    std::stringstream ss;
    ss << verbose() << "$DimAzimuth=" << mMaxNr << "$Axial=" << (axial() ? "T" : "F");
    return ss.str();
}

#include "Geodesy.h"
RDMatPP Element::formThetaMat() const {
    RDMatPP theta;
    int ipnt = 0;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &sz = mPoints[ipnt++]->getCoords();
            theta(ipol, jpol) = Geodesy::theta(sz);
        }
    }
    return theta;
}

RDMatXX Element::getCoordsOnSide(int side) const {
    int ipol0 = 0, ipol1 = 0, jpol0 = 0, jpol1 = 0;
    if (side == 0) {
        ipol0 = 0;
        ipol1 = nPol;
        jpol0 = jpol1 = 0;
    } else if (side == 1) {
        ipol0 = ipol1 = nPol;
        jpol0 = 0;
        jpol1 = nPol;
    } else if (side == 2) {
        ipol0 = 0;
        ipol1 = nPol;
        jpol0 = jpol1 = nPol;
    } else {
        ipol0 = ipol1 = 0;
        jpol0 = 0;
        jpol1 = nPol;
    }
    RDMatXX sz = RDMatXX::Zero(2, nPntEdge);
    int ipntedge = 0;
    for (int ipol = ipol0; ipol <= ipol1; ipol++) {
        for (int jpol = jpol0; jpol <= jpol1; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            sz.col(ipntedge) = mPoints[ipnt]->getCoords();
            ipntedge++;
        }
    }
    return sz;
}

RRow3 Element::recordWF(const double phi) const {
    RRow3 u_spz;
    computeGroundMotion(phi, mWFRweights, u_spz);
    return u_spz;
}

RDCol2 Element::getCenterCrds(const double phi) const {
    RDCol2 crds;
    RDRowN dz = XMath::computeFourierAtPhi(mDz, phi);
    crds << mCenterCrds(0), mCenterCrds(1) + dz(round((nPE-1)/2));
    return crds;
}
