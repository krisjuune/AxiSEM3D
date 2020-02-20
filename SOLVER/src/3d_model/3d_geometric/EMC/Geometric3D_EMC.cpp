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
        mRLayer = 6371e3 - mRLayer;
        mRLower = 6371e3 - mRLower;
        mRUpper = 6371e3 - mRUpper;
    }

    initialize();
}

void Geometric3D_EMC::setSourceLocation(double srcLat, double srcLon, double srcDep) {
    mSrcLat = srcLat;
    mSrcLon = srcLon;
    mSrcDep = srcDep;
}

void Geometric3D_EMC::initializeOcean(double rLayer, double rUpper, double rLower, RDColX lat, RDColX lon, RDMatXX bathymetry) {
    mRLayer = rLayer;
    mRUpper = rUpper;
    mRLower = rLower;
    mGridLat = lat;
    mGridLon = lon;
    mGridData = bathymetry;

    mCartesian = true;

    mVerboseOcean = " (autogenerated seafloor)";
}

double Geometric3D_EMC::getDeltaR(double r, double theta, double phi, double rElemCenter) const {
    if (r > mRUpper || r < mRLower) { 
        return 0.;
    }

    double lat, lon;

    if (mCartesian) {
        RDCol3 rtpG;
        rtpG << r, theta, phi;
        RDCol3 rtpS = Geodesy::rotateGlob2Src(rtpG, mSrcLat, mSrcLon, mSrcDep);
        RDCol3 xyz = Geodesy::toCartesian(rtpS);
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
    XMath::interpLinear(lat, mGridLat, llat0, wlat0);
    XMath::interpLinear(lon, mGridLon, llon0, wlon0);    
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

std::string Geometric3D_EMC::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Geometric =======================" << std::endl;
    ss << "  Model Name          =   EMC" << mVerboseOcean << std::endl;
    if (mCartesian) {
        ss << "  Layer Depth (m)     =   " << 6371e3 - mRLayer << std::endl;
        ss << "  Depth Range (m)     =   [" << 6371e3 - mRUpper << ", " << 6371e3 - mRLower << "]" << std::endl;
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



