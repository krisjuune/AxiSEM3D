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
    void initializeOcean(double rLayer, double rUpper, double rLower, RDColX lat, RDColX lon, RDMatXX bathymetry);
    void setSourceLocation(double srcLat, double srcLon, double srcDep);
    double getDeltaR(double r, double theta, double phi, double rElemCenter) const;
    bool isCartesian() const {return mCartesian;};
    std::string verbose() const;
    
private:
    
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

    std::string mVerboseOcean = "";
};
