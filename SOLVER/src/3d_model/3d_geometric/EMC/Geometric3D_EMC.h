// Geometric3D_EMC.h
// created by Kuangdai on 16-May-2017 
// topography on a boundary at any depth, with IRIS-EMC format

#pragma once

#include "Geometric3D.h"
#include "eigenp.h"

class Geometric3D_EMC: public Geometric3D {
public:

    void initialize();
    void initialize(const std::vector<std::string> &params);
    void initializeFromVolumetric(double rLayer, double rUpper, double rLower, RDColX lat, RDColX lon, RDMatXX depth, bool cartesian);
    void setSourceLocation(double srcLat, double srcLon, double srcDep);
    void applyShift(double x, double y) {
        if (mCartesian) {
            mGridLat = mGridLat.array() + x;
            mGridLon = mGridLon.array() + y;
        } else {
            throw std::runtime_error("Geometric3D_EMC::initialize || "
                "3D model shift not implemented for non-Cartesian models. Shift source position instead.");
        }
    };
    double getDeltaR(double r, double theta, double phi, double rElemCenter) const;
    bool isCartesian() const {return mCartesian;};
    std::string verbose() const;
    
    void handleDiscontinuities(double r, double theta, double phi, bool geographic, bool cartesian, double rElemCenter, double lat0, double lat1, 
        double lon0, double lon1, double dep0, double dep1, std::vector<bool> &use, bool v) const;
    
private:
    bool isAboveDiscontinuity(double lat, double lon, double r) const;
    
    // radius
    double mRLayer;
    double mRLower;
    double mRUpper;
    
    // file
    std::string mFileName;
    std::string mVarName;
    
    // factor
    double mFactor = 1.0;
    
    // use geocentric or geographic
    bool mGeographic = false;
    bool mCartesian = false;
    
    // source coords
    double mSrcLat;
    double mSrcLon;
    double mSrcDep;

    // data
    RDMatXX mGridData;
    RDColX mGridLat;
    RDColX mGridLon;

    std::string mVerboseType = "EMC";
};
