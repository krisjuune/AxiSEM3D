// Element.h
// created by Kuangdai on 27-Mar-2016
// base class of AxiSEM3D spectral elements

#pragma once

class Point;
class Gradient;
class PRT;
class Acoustic;

#include "eigenc.h"
#include "eigenp.h"

class Element {
public:
    Element(Gradient *grad, PRT *prt, const std::array<Point *, nPntElem> &points);
    virtual ~Element();

    // compute stiffness term
    virtual void computeStiff() const = 0;

    // measure cost
    virtual double measure(int count) const = 0;

    // test stiffness
    virtual void test() const = 0;

    // compute Real displacement, used by receiver
    virtual void computeGroundMotion(Real phi, const RMatPP &weights, RRow3 &u_spz) const = 0;
    virtual void computeStrain(Real phi, const RMatPP &weights, RRow6 &strain) const = 0;
    virtual void computeCurl(Real phi, const RMatPP &weights, RRow3 &curl) const = 0;
    virtual void forceTIso() = 0;

    // side-wise
    virtual void feedDispOnSide(int side, CMatXX_RM &buffer, int row) const = 0;
    RDMatXX getCoordsOnSide(int side) const;

    // verbose
    virtual std::string verbose() const = 0;

    // get point ptr
    const Point *getPoint(int index) const {return mPoints[index];};

    // source
    void addSourceTerm(const arPP_CMatX3 &source) const;

    // get nr
    int getMaxNr() const {return mMaxNr;};
    int getMaxNu() const {return mMaxNu;};

    // signature for cost measurement
    std::string costSignature() const;

    // axial
    bool axial() const;

    virtual bool fluid() const {return false;};

    virtual Acoustic * getAcoustic() const {
        throw std::runtime_error("Element::getAcoustic || Error.");
    };

    // form theta for TIso
    RDMatPP formThetaMat() const;

    // reset
    virtual void resetZero() {};

    // for wavefield recorder
    void addWFRweights(RDMatPP weights, RDCol2 crds, RDMatXN dz) {
        mWFRweights = weights;
        mCenterCrds = crds;
        mDz = dz;
    }
    RRow3 recordWF(const double phi) const;
    RDCol2 getCenterCrds(const double phi) const;

protected:
    int mMaxNu;
    int mMaxNr;
    std::array<Point *, nPntElem> mPoints;
    Gradient *mGradient;
    PRT *mPRT;

    // flags
    bool mHasPRT;

    RDMatPP mWFRweights;
    RDCol2 mCenterCrds;
    RDMatXN mDz;

public:
    // domain tag, mainly for debug
    void setDomainTag(int tag) {mDomainTag = tag;};
    int getDomainTag() const {return mDomainTag;};

private:
    int mDomainTag;
};
