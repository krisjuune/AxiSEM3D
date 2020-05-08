// Volumetric3D_EMC.h
// created by Kuangdai on 16-May-2017 
// genetral Volumetric3D model with IRIS-EMC format

#pragma once

#include "Volumetric3D.h"
#include "eigenp.h"

class Geometric3D;

class Volumetric3D_EMC: public Volumetric3D {
public:

    void initialize();
    void initialize(const std::vector<std::string> &params);
    void setSourceLocation(double srcLat, double srcLon, double srcDep);
    void applyShift(double x, double y, double z) {
        if (mCartesian) {
            mGridDep = mGridDep.array() - z;
            mGridLat = mGridLat.array() + x;
            mGridLon = mGridLon.array() + y;
        } else {
            throw std::runtime_error("Volumetric3D_EMC::initialize || "
                "3D model shift not implemented for non-Cartesian models. Shift source position instead.");
        }
    };
    bool get3dPropertiesInternal(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties, 
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values, bool isFluid) const;
    bool makeFluid3D() const {return true;};  
    bool isEMC() const {return true;};
    bool isCartesian() const {return mCartesian;};
    std::string verbose() const;
    
    std::string getVarName() const {return mVarName;};
    void findDiscontinuity(RDColX &lat, RDColX &lon, RDMatXX &depth, bool &cartesian, double val, 
        double minDepth, double maxDepth, std::string &compType, bool from_bottom) const;
    void getSourceLocation(double &srcLat, double &srcLon, double &srcDep) const {
        srcLat = mSrcLat;
        srcLon = mSrcLon;
        srcDep = mSrcDep;
    };
    void setDiscontinuities(const std::vector<Geometric3D *> discontinuities) {
        mDiscontinuities = discontinuities;
        if (mDiscontinuities.size() > 0) mVerticalDiscontinuities = false;
    };
    
private:
    
    // file
    std::string mFileName;
    std::string mVarName;
    
    // property
    MaterialProperty mMaterialProp;
    MaterialRefType mReferenceType;
    
    // factor
    double mFactor = 1.0;
    
    // use geocentric or geographic
    bool mCartesian = false;
    bool mGeographic = false;
    
    // source coords
    double mSrcLat;
    double mSrcLon;
    double mSrcDep;
    
    // one-file-per-depth format
    bool mOneFilePerDepth = false;
    
    // consider vertical disc or not
    bool mVerticalDiscontinuities = true;
    
    // data
    std::vector<RDMatXX> mGridData;
    RDColX mGridDep;
    RDColX mGridLat;
    RDColX mGridLon;
    
    // special model flag
    // abs -- use absolute value of the perturbations
    // pow -- use power of the perturbations
    // to use these special model flags, reference type cannot be Absolute
    // for "abs" flag, the following factor is used to change the sign
    // for "pow" flag, the following factor specifies the power, e.g.,
    // 2.0 means "squared", which makes the model sharper  
    // 0.5 means "sqrt", which makes the model smoother
    // the pow flag keeps the sign and global absolute maximum of the perturbations
    std::string mModelFlag = "none";
    double mModelFlagFactor = 1.0;
    
    std::vector<Geometric3D *> mDiscontinuities;
};

