// Volumetric3D_block.h
// created by Claudia on 25-Jul-2019 
// a cuboid heterogeneity

#pragma once
#include "Volumetric3D.h"
#include "eigenp.h"

class Volumetric3D_block: public Volumetric3D {
public:
    
    void initialize(const std::vector<std::string> &params);
    
    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values, bool isFluid) const;

    std::string verbose() const;
    
    void setSourceLocation(double srcLat, double srcLon, double srcDep) {
        mSrcLat = srcLat;
        mSrcLon = srcLon;
        mSrcDep = srcDep;
    }
    
    bool makeFluid3D() const {return mFluid;};
    
private:
    // property name
    MaterialProperty mMaterialProp;
    
    // value inside the bubble
    double mValueInside; 

    // reference type
    MaterialRefType mReferenceType;
    
    // radius of the bubble
    double mRadius;

    // center of the bubble
    double mX1, mX2, mY1, mY2, mZ1, mZ2;
    
    // source-centered
    bool mCartesian = true;
    bool mSourceCentered = true;
    double mSrcLat = 0.;
    double mSrcLon = 0.;
    double mSrcDep = 0.;
    
    // 3d fluid
    bool mFluid = false;
    
    // halfwidth at half maximum of Gaussian 
    // how the perturbation fades outside the block
    double mHWHM = -1.;
    
    // temp variables for performance
    RDCol3 mXyzCorner1, mXyzCorner2;
};
