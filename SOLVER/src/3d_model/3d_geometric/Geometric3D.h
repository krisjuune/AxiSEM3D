// Geometric3D.h
// created by Kuangdai on 4-Jun-2016 
// base class of geometric 3D models, such as 
// topography, ellipticity and undulating Moho

#pragma once
#include <string>
#include <vector>
#include <eigenp.h>

class Parameters;
class AutoGeometricParams;

class Geometric3D {
public:

    virtual ~Geometric3D() {finalize();};
    
    // initialize internal variables if needed
    virtual void initialize() {};
    virtual void initialize(const std::vector<std::string> &params) {initialize();};
    virtual void initializeOcean(double rLayer, double rUpper, double rLower, RDColX lat, RDColX lon, RDMatXX bathymetry) {};
    virtual void setSourceLocation(double srcLat, double srcLon, double srcDep) {};

    // finalize internal variables if needed
    virtual void finalize() {};
    
    // get undulation (deltaR)
    // IMPORTANT NOTES: 
    // a) This function should be realized such that r/theta/phi are the geocentric 
    //    coordinates, without rotating the source to the north pole.  
    //    For models given in geographic coordinates, geocentric-to-geographic  
    //    conversions from (theta, phi) to (lat, lon) has to be performed internally.
    // b) All models should be defined independently with respect to the perfect sphere.
    virtual double getDeltaR(double r, double theta, double phi, double rElemCenter) const = 0;
    virtual bool isCartesian() const = 0;

    // verbose
    virtual std::string verbose() const = 0;
    
    // build from input parameters
    static void buildInparam(std::vector<Geometric3D *> &models,
        const Parameters &par, std::vector<AutoGeometricParams *> Vol2GeomModels, 
        double srcLat, double srcLon, double srcDep, int verbose);

};

