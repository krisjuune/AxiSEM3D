// Geometric3D_EMC.cpp
// created by Kuangdai on 16-May-2017 
// topography on a boundary at any depth, with IRIS-EMC format

#include "Geometric3D_EMC.h"
#include <sstream>
#include "XMPI.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "NetCDF_Reader.h"
#include "NetCDF_ReaderAscii.h"
#include "eigenp.h"

void Geometric3D_EMC::initialize() {
    // EMC is in float
    Eigen::Matrix<float, Eigen::Dynamic, 1> flat, flon;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> fdata;
    std::string latStr = mCartesian ? "x" : "latitude";
    std::string lonStr = mCartesian ? "y" : "longitude";
    
    // read file
    if (XMPI::root()) {
        std::string fname = Parameters::sInputDirectory + "/" + mFileName;
        if (NetCDF_Reader::checkNetCDF_isAscii(fname)) {
            NetCDF_ReaderAscii reader;
            reader.open(fname);
            reader.read1D(latStr, flat);
            reader.read1D(lonStr, flon);
            reader.read2D(mVarName, fdata);
            reader.close();
        } else {
            NetCDF_Reader reader;
            reader.open(fname);
            reader.read1D(latStr, flat);
            reader.read1D(lonStr, flon);
            reader.read2D(mVarName, fdata);
            reader.close();
        }
        if (fdata.rows() != flat.size() || fdata.cols() != flon.size()) {
            throw std::runtime_error("Geometric3D_EMC::initialize || "
                "Inconsistent data dimensions || File = " + fname);
        }
        if (!XMath::sortedAscending(flat) || !XMath::sortedAscending(flon)) {
            throw std::runtime_error("Geometric3D_EMC::initialize || "
                "Grid coordinates are not sorted ascendingly || File = " + fname);
        }
    }
    
    // broadcast
    XMPI::bcastEigen(flat);
    XMPI::bcastEigen(flon);
    XMPI::bcastEigen(fdata);
    
    // to double
    mGridLat = flat.cast<double>();
    mGridLon = flon.cast<double>();
    mGridData = fdata.cast<double>();

    if (!mCartesian) {
        // to SI
        mGridData *= 1e3;
    }

    // apply factor
    mGridData *= mFactor;

    ////////// plot computed data ////////////  
    // std::fstream fsdr;
    // fsdr.open("/Users/kuangdai/Desktop/crust1/cmb.txt", std::fstream::out);
    // double r = mRLayer; // double r = mRSurf;
    // int intGrid = 1;
    // for (int i = 0; i <= 180 * intGrid; i++) {
    //     double theta = i * degree / intGrid;
    //     for (int j = -180 * intGrid; j <= 180 * intGrid; j++) {
    //         double phi = j * degree / intGrid;
    //         if (phi < 0) phi += 2. * pi;
    //         fsdr << getDeltaR(r, theta, phi, r) << " ";
    //     }
    //     fsdr << std::endl;
    // }   
    // fsdr.close();
    // exit(0);
    ////////// plot computed data ////////////  
}

void Geometric3D_EMC::initialize(const std::vector<std::string> &params) {
    if (params.size() < 5) throw std::runtime_error("Geometric3D_EMC::initialize || "
        "Not enough parameters to initialize a Geometric3D_EMC object, at least 5 needed.");
    
    const std::string source = "Geometric3D_EMC::initialize";

    // initialize location
    Parameters::castValue(mRLayer, params[0], source);
    Parameters::castValue(mRLower, params[1], source);
    Parameters::castValue(mRUpper, params[2], source);

    // initialize data
    Parameters::castValue(mFileName, params[3], source);
    Parameters::castValue(mVarName, params[4], source);
    
    try {
        int ipar = 5;
        Parameters::castValue(mFactor, params.at(ipar++), source);
        Parameters::castValue(mGeographic, params.at(ipar++), source);
        Parameters::castValue(mCartesian, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }

    mRLayer *= 1e3;
    mRLower *= 1e3;
    mRUpper *= 1e3;

    if (mCartesian) {
        mRLayer = Geodesy::getROuter() - mRLayer;
        mRLower = Geodesy::getROuter() - mRLower;
        mRUpper = Geodesy::getROuter() - mRUpper;
    }

    initialize();
}

void Geometric3D_EMC::setSourceLocation(double srcLat, double srcLon, double srcDep) {
    mSrcLat = srcLat;
    mSrcLon = srcLon;
    mSrcDep = srcDep;
}

void Geometric3D_EMC::initializeFromVolumetric(double rLayer, double rUpper, double rLower, RDColX lat, RDColX lon, RDMatXX depth, bool cartesian) {
    mRLayer = Geodesy::getROuter() - rLayer;
    mRUpper = Geodesy::getROuter() - rUpper;
    mRLower = Geodesy::getROuter() - rLower;
    mGridLat = lat;
    mGridLon = lon;
    mGridData = depth;

    mCartesian = cartesian;

    mVerboseType = "Discontinuity from volumetric model (implementated as EMC)";
}

double Geometric3D_EMC::getDeltaR(double r, double theta, double phi, double rElemCenter) const {
    if (r > mRUpper || r < mRLower) { 
        return 0.;
    }

    double lat, lon;

    if (mCartesian) {
        RDCol3 rtpG;
        rtpG << r, theta, phi;
        RDCol3 xyz = Geodesy::Glob2Cartesian(rtpG, mSrcLat, mSrcLon, mSrcDep);
        lat = xyz(0);
        lon = xyz(1);
        r = xyz(2);
    } else {
        // to geocentric
        if (mGeographic) {
            // which radius to use?
            theta = pi / 2. - Geodesy::theta2Lat_d(theta, 0.) * degree;
        }

        // regularise
        lat = 90. - theta / degree;
        lon = phi / degree;
        XMath::checkLimits(lat, -90., 90.);
        if (mGridLon[0] < 0.) {
            // lon starts from -180.
            if (lon > 180.) {
                lon -= 360.;
            }
            XMath::checkLimits(lon, -180., 180.);
        } else {
            // lon starts from 0.
            XMath::checkLimits(lon, 0., 360.);
        }
    }
    
    // interpolation on sphere
    int llat0, llon0, llat1, llon1;
    double wlat0, wlon0, wlat1, wlon1;
    double tol = 0.01;
    XMath::interpLinearRobust(lat, mGridLat, llat0, wlat0, tol);
    XMath::interpLinearRobust(lon, mGridLon, llon0, wlon0, tol);    
    if (llat0 < 0 || llon0 < 0) {
        return 0.;
    }
    
    llat1 = llat0 + 1;
    llon1 = llon0 + 1;
    wlat1 = 1. - wlat0;
    wlon1 = 1. - wlon0;
    
    double dr = 0.;
    dr += mGridData(llat0, llon0) * wlat0 * wlon0;
    dr += mGridData(llat1, llon0) * wlat1 * wlon0;
    dr += mGridData(llat0, llon1) * wlat0 * wlon1;
    dr += mGridData(llat1, llon1) * wlat1 * wlon1;
    
    // interpolation along radius    
    if (rElemCenter < mRLayer) {
        return dr / (mRLayer - mRLower) * (r - mRLower);
    } else {
        return dr / (mRUpper - mRLayer) * (mRUpper - r);
    }
}

void Geometric3D_EMC::handleDiscontinuities(double r, double theta, double phi, bool geographic, bool cartesian, double rElemCenter, 
    double lat0, double lat1, double lon0, double lon1, double dep0, double dep1, std::vector<bool> &use, bool v) const {

    double r0 = Geodesy::getROuter() - dep1;
    double r1 = Geodesy::getROuter() - dep0;
    
    if (r0 > mRUpper || r1 < mRLower) {
        return;
    }

    double lat, lon;
    if (mCartesian) {
        RDCol3 rtpG;
        rtpG << r, theta, phi;
        RDCol3 xyz = Geodesy::Glob2Cartesian(rtpG, mSrcLat, mSrcLon, mSrcDep);
        lat = xyz(0);
        lon = xyz(1);
        r = xyz(2);
    } else {
        // to geocentric
        if (mGeographic) {
            // which radius to use?
            theta = pi / 2. - Geodesy::theta2Lat_d(theta, 0.) * degree;
        }

        // regularise
        lat = 90. - theta / degree;
        lon = phi / degree;
        XMath::checkLimits(lat, -90., 90.);
        if (mGridLon[0] < 0.) {
            // lon starts from -180.
            if (lon > 180.) {
                lon -= 360.;
            }
            XMath::checkLimits(lon, -180., 180.);
        } else {
            // lon starts from 0.
            XMath::checkLimits(lon, 0., 360.);
        }
    }
    
    std::vector<RDCol3> crds;
    if (mCartesian && !cartesian) {
        std::vector<RDCol3> rtp(8, RDCol3::Zero(3,1));
        
        lon0 = Geodesy::lon2Phi(lon0);
        lon1 = Geodesy::lon2Phi(lon1);
        
        double lat00 = Geodesy::lat2Theta_r(lat0, r0);
        double lat01 = Geodesy::lat2Theta_r(lat0, r1);
        double lat10 = Geodesy::lat2Theta_r(lat1, r0);
        double lat11 = Geodesy::lat2Theta_r(lat1, r1);
        
        rtp[0](0) = lat00; rtp[0](1) = lon0; rtp[0](2) = r0;
        rtp[1](0) = lat10; rtp[1](1) = lon0; rtp[1](2) = r0;
        rtp[2](0) = lat00; rtp[2](1) = lon1; rtp[2](2) = r0;
        rtp[3](0) = lat10; rtp[3](1) = lon1; rtp[3](2) = r0;
        rtp[4](0) = lat01; rtp[4](1) = lon0; rtp[4](2) = r1;
        rtp[5](0) = lat11; rtp[5](1) = lon0; rtp[5](2) = r1;
        rtp[6](0) = lat01; rtp[6](1) = lon1; rtp[6](2) = r1;
        rtp[7](0) = lat11; rtp[7](1) = lon1; rtp[7](2) = r1;
        
        for (auto &crd: rtp) {
            crds.push_back(Geodesy::Glob2Cartesian(crd, mSrcLat, mSrcLon, mSrcDep));
        }
    } else if (!mCartesian && cartesian) {
        std::vector<RDCol3> xyz(8, RDCol3::Zero(3,1));
        
        xyz[0](0) = lat0; xyz[0](1) = lon0; xyz[0](2) = r0;
        xyz[1](0) = lat1; xyz[1](1) = lon0; xyz[1](2) = r0;
        xyz[2](0) = lat0; xyz[2](1) = lon1; xyz[2](2) = r0;
        xyz[3](0) = lat1; xyz[3](1) = lon1; xyz[3](2) = r0;
        xyz[4](0) = lat0; xyz[4](1) = lon0; xyz[4](2) = r1;
        xyz[5](0) = lat1; xyz[5](1) = lon0; xyz[5](2) = r1;
        xyz[6](0) = lat0; xyz[6](1) = lon1; xyz[6](2) = r1;
        xyz[7](0) = lat1; xyz[7](1) = lon1; xyz[7](2) = r1;
        
        for (auto &crd: xyz) {
            RDCol3 rtp = Geodesy::Cartesian2Glob(crd, mSrcLat, mSrcLon, mSrcDep);
            rtp(0) = Geodesy::theta2Lat_r(rtp(0), rtp(2));
            rtp(1) = Geodesy::phi2Lon(rtp(1));
            crds.push_back(rtp);
        }
    } else {
        RDCol3 crd = RDCol3::Zero(3,1);
        
        crd(0) = lat0; crd(1) = lon0; crd(2) = r0;
        crds.push_back(crd);
        crd(0) = lat1; crd(1) = lon0; crd(2) = r0;
        crds.push_back(crd);
        crd(0) = lat0; crd(1) = lon1; crd(2) = r0;
        crds.push_back(crd);
        crd(0) = lat1; crd(1) = lon1; crd(2) = r0;
        crds.push_back(crd);
        crd(0) = lat0; crd(1) = lon0; crd(2) = r1;
        crds.push_back(crd);
        crd(0) = lat1; crd(1) = lon0; crd(2) = r1;
        crds.push_back(crd);
        crd(0) = lat0; crd(1) = lon1; crd(2) = r1;
        crds.push_back(crd);
        crd(0) = lat1; crd(1) = lon1; crd(2) = r1;
        crds.push_back(crd);
    }
    
    if (!mCartesian) {
        if  (mGridLon[0] < 0.) {
            // lon starts from -180.
            for (auto &crd: crds) {
                if (crd(1) > 180.) crd(1) -= 360;
                XMath::checkLimits(crd(1), -180., 180.);
            }
        } else {
            for (auto &crd: crds) {
                if (crd(1) < 0.) crd(1) += 360;
                XMath::checkLimits(crd(1), 0., 360.);
            }
        }
    }
    bool isQueryPointAbove = isAboveDiscontinuity(lat, lon, rElemCenter);
    for (int i = 0; i < 8; i++) {
        bool isCornerPointAbove = isAboveDiscontinuity(crds[i](0), crds[i](1), crds[i](2));
        if (isQueryPointAbove != isCornerPointAbove) use[i] = false;
    }
}

bool Geometric3D_EMC::isAboveDiscontinuity(double lat, double lon, double r) const {
    int llat0, llon0, llat1, llon1;
    double wlat0, wlon0, wlat1, wlon1;
    XMath::interpLinearRobust(lat, mGridLat, llat0, wlat0, 0.01);
    XMath::interpLinearRobust(lon, mGridLon, llon0, wlon0, 0.01);    
    if (llat0 < 0 || llon0 < 0) {
        throw std::runtime_error("Geometric3D_EMC::isAboveDiscontinuity || Discontinuity is smaller than mesh.");
    }
    
    llat1 = llat0 + 1;
    llon1 = llon0 + 1;
    wlat1 = 1. - wlat0;
    wlon1 = 1. - wlon0;
    
    double dr = 0.;
    dr += mGridData(llat0, llon0) * wlat0 * wlon0;
    dr += mGridData(llat1, llon0) * wlat1 * wlon0;
    dr += mGridData(llat0, llon1) * wlat0 * wlon1;
    dr += mGridData(llat1, llon1) * wlat1 * wlon1;
    
    return r > mRLayer + dr;
} 

std::string Geometric3D_EMC::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Geometric =======================" << std::endl;
    ss << "  Model Name          =   " << mVerboseType << std::endl;
    if (mCartesian) {
        ss << "  Layer Depth (m)     =   " << Geodesy::getROuter() - mRLayer << std::endl;
        ss << "  Depth Range (m)     =   [" << Geodesy::getROuter() - mRUpper << ", " << Geodesy::getROuter() - mRLower << "]" << std::endl;
        ss << "  X Coords Range (km) =   [" << mGridLat.minCoeff() / 1e3 << ", " << mGridLat.maxCoeff() / 1e3 << "]" << std::endl;
        ss << "  Y Coords Range (km) =   [" << mGridLon.minCoeff() / 1e3 << ", " << mGridLon.maxCoeff() / 1e3 << "]" << std::endl;
        ss << "  Data Range (m)      =   [" << mGridData.minCoeff() << ", " << mGridData.maxCoeff() << "]" << std::endl;
    } else {
        ss << "  Layer Radius (km)   =   " << mRLayer / 1e3 << std::endl;
        ss << "  Radius Range (km)   =   [" << mRLower / 1e3 << ", " << mRUpper / 1e3 << "]" << std::endl;
        ss << "  Latitude Range      =   [" << mGridLat.minCoeff() << ", " << mGridLat.maxCoeff() << "]" << std::endl;
        ss << "  Longitude Range     =   [" << mGridLon.minCoeff() << ", " << mGridLon.maxCoeff() << "]" << std::endl;
        ss << "  Data Range (km)     =   [" << mGridData.minCoeff() / 1e3 << ", " << mGridData.maxCoeff() / 1e3 << "]" << std::endl;
    }
    ss << "  Data File           =   " << mFileName << std::endl;
    ss << "  Variable Name       =   " << mVarName << std::endl;
    ss << "  Num. Dim1           =   " << mGridLat.size() << std::endl;
    ss << "  Num. Dim2           =   " << mGridLon.size() << std::endl;
    ss << "  Factor              =   " << mFactor << std::endl;
    ss << "  Use Geographic      =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "  Cartesian           =   " << (mCartesian ? "YES" : "NO") << std::endl;
    ss << "======================= 3D Geometric =======================\n" << std::endl;
    return ss.str();
}



