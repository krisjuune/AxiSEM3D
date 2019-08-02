// Volumetric3D_block.cpp
// created by Claudia on 25-Jul-2019 
// a cuboid heterogeneity

#include "Volumetric3D_block.h"
#include "Parameters.h"
#include "Geodesy.h"
#include <boost/algorithm/string.hpp>
#include <sstream>

void Volumetric3D_block::initialize(const std::vector<std::string> &params) {
    // need at least 9 parameters to make a block
    if (params.size() < 9) {
        throw std::runtime_error("Volumetric3D_block::initialize || "
            "Not enough parameters for a cuboid heterogeneity. Need 9 at least.");
    }
        
    const std::string source = "Volumetric3D_block::initialize";
    
    // property name
    bool found = false;
    for (int i = 0; i < Volumetric3D::MaterialPropertyString.size(); i++) {
        if (boost::iequals(params[0], Volumetric3D::MaterialPropertyString[i])) {
            mMaterialProp = Volumetric3D::MaterialProperty(i);
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_block::initialize || "
            "Unknown material property, name = " + params[0]);
    }
    
    // reference type
    found = false;
    for (int i = 0; i < Volumetric3D::MaterialRefTypeString.size(); i++) {
        if (boost::iequals(params[1], Volumetric3D::MaterialRefTypeString[i]) ||
            boost::iequals(params[1], Volumetric3D::MaterialRefTypeStringShort[i])) {
            mReferenceType = Volumetric3D::MaterialRefType(i);
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_block::initialize || "
            "Unknown material reference type, type = " + params[1]);
    }
    
    // value inside
    Parameters::castValue(mValueInside, params[2], source);
    
    // location
    Parameters::castValue(mX1, params[3], source);
    Parameters::castValue(mX2, params[4], source);
    Parameters::castValue(mY1, params[5], source);
    Parameters::castValue(mY2, params[6], source);
    Parameters::castValue(mZ1, params[7], source); mZ1 *= 1e3;
    Parameters::castValue(mZ2, params[8], source); mZ2 *= 1e3;

    // optional
    try {
        int ipar = 9;
        Parameters::castValue(mCartesian, params.at(ipar++), source);
        Parameters::castValue(mSourceCentered, params.at(ipar++), source);
        Parameters::castValue(mFluid, params.at(ipar++), source);
        Parameters::castValue(mHWHM, params.at(ipar++), source); mHWHM *= 1e3;
    } catch (std::out_of_range) {
        // nothing
    }    
    
    // save values for verbose output
    mXyzCorner1(0) = mX1;
    mXyzCorner1(1) = mY1;
    mXyzCorner1(2) = mZ1;
        
    mXyzCorner2(0) = mX2;
    mXyzCorner2(1) = mY2;
    mXyzCorner2(2) = mZ2;
    
    // compute xyz and put sides in right order
    RDCol3 xyzCorner1, xyzCorner2;
    if (!mCartesian) {
        RDCol3 rtpCorner1, rtpCorner2;
        if (mSourceCentered) {
            RDCol3 rtpCorner1Src, rtpCorner2Src;
            rtpCorner1Src(0) = Geodesy::getROuter() - mZ1;
            rtpCorner1Src(1) = mX1 * degree;
            rtpCorner1Src(2) = mY1 * degree;
            rtpCorner1 = Geodesy::rotateSrc2Glob(rtpCorner1Src, mSrcLat, mSrcLon, mSrcDep);
            
            rtpCorner2Src(0) = Geodesy::getROuter() - mZ2;
            rtpCorner2Src(1) = mX2 * degree;
            rtpCorner2Src(2) = mY2 * degree;
            rtpCorner2 = Geodesy::rotateSrc2Glob(rtpCorner2Src, mSrcLat, mSrcLon, mSrcDep);
        } else {
            rtpCorner1(0) = Geodesy::getROuter() - mZ1;
            rtpCorner1(1) = Geodesy::lat2Theta_d(mX1, mZ1);
            rtpCorner1(2) = Geodesy::lon2Phi(mY1);
            
            rtpCorner2(0) = Geodesy::getROuter() - mZ2;
            rtpCorner2(1) = Geodesy::lat2Theta_d(mX2, mZ2);
            rtpCorner2(2) = Geodesy::lon2Phi(mY2);
        }
        xyzCorner1 = Geodesy::toCartesian(rtpCorner1);
        xyzCorner2 = Geodesy::toCartesian(rtpCorner2);
    } else {
        xyzCorner1(0) = 1e3 * mX1;
        xyzCorner1(1) = 1e3 * mY1;
        xyzCorner1(2) = 1e3 * mZ1;
        
        xyzCorner2(0) = 1e3 * mX2;
        xyzCorner2(1) = 1e3 * mY2;
        xyzCorner2(2) = 1e3 * mZ2;
    }
    mX1 = std::min({xyzCorner1(0), xyzCorner2(0)});
    mX2 = std::max({xyzCorner1(0), xyzCorner2(0)});
    mY1 = std::min({xyzCorner1(1), xyzCorner2(1)});
    mY2 = std::max({xyzCorner1(1), xyzCorner2(1)});
    mZ1 = std::min({xyzCorner1(2), xyzCorner2(2)});
    mZ2 = std::max({xyzCorner1(2), xyzCorner2(2)});
    
    // use 20% of radius for HWHM if not specified
    if (mHWHM < 0.) {
        mHWHM = mRadius * .2;
    }
    
    // for Absolute models
    if (mReferenceType == Volumetric3D::MaterialRefType::Absolute) {
        // decay is not allowed
        mHWHM = 0.;
        // convert to SI
        mValueInside *= MaterialPropertyAbsSI[mMaterialProp];
    }
}

bool Volumetric3D_block::get3dProperties(double r, double theta, double phi, double rElemCenter,
    std::vector<MaterialProperty> &properties, 
    std::vector<MaterialRefType> &refTypes,
    std::vector<double> &values, bool isFluid) const {

    // header
    properties = std::vector<MaterialProperty>(1, mMaterialProp);
    refTypes = std::vector<MaterialRefType>(1, mReferenceType);
    values = std::vector<double>(1, 0.);
    
    // distance from point to axis
    RDCol3 rtpTarget;
    rtpTarget(0) = r;
    rtpTarget(1) = theta;
    rtpTarget(2) = phi;
    const RDCol3 &xyzTarget = Geodesy::toCartesian(rtpTarget);
    
    bool insideX = (mX1 <= xyzTarget(0) && xyzTarget(0) <= mX2);
    bool insideY = (mY1 <= xyzTarget(1) && xyzTarget(1) <= mY2);
    bool insideZ = (mZ1 <= xyzTarget(2) && xyzTarget(2) <= mZ2);
    
    double dx, dy, dz;
    if (insideX) {
        dx = 0;
    } else {
        dx = std::min({std::abs(xyzTarget(0) - mX1), std::abs(xyzTarget(0) - mX2)});
    }
    if (insideY) {
        dy = 0;
    } else {
        dy = std::min({std::abs(xyzTarget(1) - mY1), std::abs(xyzTarget(1) - mY2)});
    }
    if (insideZ) {
        dz = 0;
    } else {
        dz = std::min({std::abs(xyzTarget(2) - mZ1), std::abs(xyzTarget(2) - mZ2)});
    }
    
    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // outside range
    if (distance > 4. * mHWHM) {
        return false;
    }
    
    // compute Gaussian
    double stddev = mHWHM / sqrt(2. * log(2.));
    double gaussian = mValueInside * exp(-distance * distance / (stddev * stddev * 2.));
    
    // set perturbations    
    values[0] = gaussian;
    return true;    
}

std::string Volumetric3D_block::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name                   =   block" << std::endl;
    ss << "  Material Property            =   " << MaterialPropertyString[mMaterialProp] << std::endl;
    ss << "  Reference Type               =   " << MaterialRefTypeString[mReferenceType] << std::endl;
    if (mReferenceType == Volumetric3D::MaterialRefType::Absolute) {
        ss << "  Value Inside                 =   " << mValueInside / MaterialPropertyAbsSI[mMaterialProp] << std::endl;
    } else {
        ss << "  Value Inside                 =   " << mValueInside << std::endl;
    }
    ss << "  X (km) or Lat or Theta (deg) =   " << mXyzCorner1(0) << " - " << mXyzCorner2(0) << std::endl;
    ss << "  Y (km) or Lon or Phi (deg)   =   " << mXyzCorner1(1) << " - " << mXyzCorner2(1) << std::endl;
    ss << "  Depth (km)                   =   " << mXyzCorner1(2) << " - " << mXyzCorner2(2) << std::endl;
    ss << "  Cartesian                    =   " << (mCartesian ? "YES" : "NO") << std::endl;
    ss << "  Source-centered              =   " << (mSourceCentered ? "YES" : "NO") << std::endl;
    ss << "  HWHM / km                    =   " << mHWHM / 1e3 << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}
