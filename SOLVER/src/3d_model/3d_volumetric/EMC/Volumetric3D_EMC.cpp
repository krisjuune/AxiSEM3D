// Volumetric3D_EMC.cpp
// created by Kuangdai on 16-May-2017 
// genetral Volumetric3D model with IRIS-EMC format

#include "Volumetric3D_EMC.h"
#include <sstream>
#include "XMPI.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "NetCDF_Reader.h"
#include "NetCDF_ReaderAscii.h"

#include "Geometric3D.h"
#include <dirent.h>
#include <sys/types.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <algorithm>

bool compareFunc(const std::pair<std::string, float> &a, 
    const std::pair<std::string, float> &b) {
    return a.second < b.second;
}

void Volumetric3D_EMC::initialize() {
    // meta data
    Eigen::Matrix<float, Eigen::Dynamic, 1> fdata, fdep, flat, flon;
    if (XMPI::root()) {
        std::vector<size_t> dims;
        std::string fname = Parameters::sInputDirectory + "/" + mFileName;
        std::string latStr = mCartesian ? "x" : "latitude";
        std::string lonStr = mCartesian ? "y" : "longitude";
 
        if (mOneFilePerDepth) {
            // get all files
            DIR *dir = opendir(fname.c_str());
            if (!dir) {
                throw std::runtime_error("Volumetric3D_EMC::initialize || " 
                    "Error opening directory of data files || Directory = " + fname);    
            }
            struct dirent *entry;
            std::vector<std::pair<std::string, float>> fileDepth;
            while ((entry = readdir(dir)) != NULL) {
                // file name
                std::string fn(entry->d_name);
                if (fn.find(".grd") == std::string::npos) {
                    continue;
                }
                // depth
                std::string depthStr(fn);
                boost::replace_first(depthStr, ".grd", "");
                std::size_t found = depthStr.find_last_of("_");
                if (found == std::string::npos) {
                    throw std::runtime_error("Volumetric3D_EMC::initialize || " 
                        "Error processing file name || File = " + fn);    
                }
                depthStr = depthStr.substr(found + 1);
                float depth = 0.;
                try {
                    depth = boost::lexical_cast<float>(depthStr);
                } catch (std::exception) {
                    throw std::runtime_error("Volumetric3D_EMC::initialize || " 
                        "Error processing file name || File = " + fn);
                }
                // add pair
                fileDepth.push_back(std::make_pair(fname + "/" + fn, depth));
            }
            closedir(dir);
            
            // sort depth
            std::sort(fileDepth.begin(), fileDepth.end(), compareFunc);
            
            // read lat and lon
            Eigen::Matrix<double, Eigen::Dynamic, 1> dlat, dlon;
            if (NetCDF_Reader::checkNetCDF_isAscii(fileDepth[0].first)) {
                NetCDF_ReaderAscii reader;
                reader.open(fileDepth[0].first);
                reader.read1D(latStr, dlat);
                reader.read1D(lonStr, dlon);
                reader.close();
            } else {
                NetCDF_Reader reader;
                reader.open(fileDepth[0].first);
                reader.read1D(latStr, dlat);
                reader.read1D(lonStr, dlon);
                reader.close();
            }
            flat = dlat.cast<float>();
            flon = dlon.cast<float>();
            
            // depths and data
            size_t depthLen = flat.size() * flon.size();   
            size_t totalLen = depthLen * fileDepth.size();
            fdep.resize(fileDepth.size());
            fdata.resize(totalLen);
            for (int i = 0; i < fileDepth.size(); i++) {
                fdep(i) = fileDepth[i].second;
                Eigen::Matrix<float, Eigen::Dynamic, 1> temp;
                if (NetCDF_Reader::checkNetCDF_isAscii(fileDepth[i].first)) {
                    NetCDF_ReaderAscii reader;
                    reader.open(fileDepth[i].first);
                    reader.readMetaData(mVarName, temp, dims);
                    reader.close();
                } else {
                    NetCDF_Reader reader;
                    reader.open(fileDepth[i].first);
                    reader.readMetaData(mVarName, temp, dims);
                    reader.close();
                }
                fdata.block(i * depthLen, 0, depthLen, 1) = temp;
            }
            dims.insert(dims.begin(), fdep.size());
            
            // perturbations are given in percentage
            fdata *= .01f;
        } else {
            if (NetCDF_Reader::checkNetCDF_isAscii(fname)) {
                NetCDF_ReaderAscii reader;
                reader.open(fname);
                reader.read1D("depth", fdep);
                reader.read1D(latStr, flat);
                reader.read1D(lonStr, flon);
                reader.readMetaData(mVarName, fdata, dims);
                reader.close();
            } else {
                NetCDF_Reader reader;
                reader.open(fname);
                reader.read1D("depth", fdep);
                reader.read1D(latStr, flat);
                reader.read1D(lonStr, flon);
                reader.readMetaData(mVarName, fdata, dims);
                reader.close();
            }
        }
        
        // check dimensions
        if (dims.size() != 3) {
            throw std::runtime_error("Volumetric3D_EMC::initialize || Inconsistent data dimensions || "
                "File/Directory = " + fname);
        }
        if (dims[0] != fdep.size() || dims[1] != flat.size() || dims[2] != flon.size()) {
            throw std::runtime_error("Volumetric3D_EMC::initialize || Inconsistent data dimensions || "
                "File/Directory = " + fname);
        }
        if (!XMath::sortedAscending(fdep) || !XMath::sortedAscending(flat) || !XMath::sortedAscending(flon)) {
            throw std::runtime_error("Volumetric3D_EMC::initialize || Grid coordinates are not sorted ascendingly || "
                "File/Directory = " + fname);
        }
    }
    XMPI::bcastEigen(fdep);
    XMPI::bcastEigen(flat);
    XMPI::bcastEigen(flon);
    XMPI::bcastEigen(fdata);
    
    mGridDep = fdep.cast<double>();
    mGridLat = flat.cast<double>();
    mGridLon = flon.cast<double>();
    RDColX data = fdata.cast<double>();
    
    // SI
    if (!mCartesian) mGridDep *= 1e3;
    if (mReferenceType == Volumetric3D::MaterialRefType::Absolute) {
        // convert to SI
        data *= MaterialPropertyAbsSI[mMaterialProp];
    }
    
    // apply factor
    data *= mFactor;
    
    // special flag
    if (!boost::iequals(mModelFlag, "none")) {
        if (boost::iequals(mModelFlag, "abs")) {
            data.array() = data.array().abs();
            data *= mModelFlagFactor;
        } else if (boost::iequals(mModelFlag, "pow")) {
            double absmax = data.array().abs().maxCoeff();
            double scalefact = std::abs(absmax / std::pow(absmax, mModelFlagFactor));
            data.array() = data.array().sign() * (data.array().abs().pow(mModelFlagFactor) * scalefact);
        } else {
            throw std::runtime_error("Volumetric3D_EMC::initialize || "
                "Unknown special model flag, flag = " + mModelFlag);
        }
    }
    
    // reshape data
    int pos = 0;
    for (int i = 0; i < mGridDep.size(); i++) {
        RDMatXX mat(mGridLat.size(), mGridLon.size());
        for (int j = 0; j < mGridLat.size(); j++) {
            for (int k = 0; k < mGridLon.size(); k++) {
                mat(j, k) = data(pos++);
            }
        }
        mGridData.push_back(mat);
    }
}

void Volumetric3D_EMC::initialize(const std::vector<std::string> &params) {
    if (params.size() < 4) throw std::runtime_error("Volumetric3D_EMC::initialize || "
        "Not enough parameters to initialize a Volumetric3D_EMC object, at least 4 needed.");
    
    const std::string source = "Volumetric3D_EMC::initialize";
        
    // initialize data
    Parameters::castValue(mFileName, params[0], source);
    Parameters::castValue(mVarName, params[1], source);
    
    // property name
    bool found = false;
    for (int i = 0; i < Volumetric3D::MaterialPropertyString.size(); i++) {
        if (boost::iequals(params[2], Volumetric3D::MaterialPropertyString[i])) {
            mMaterialProp = Volumetric3D::MaterialProperty(i);
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_EMC::initialize || "
            "Unknown material property, name = " + params[2]);
    }
    
    // reference type
    found = false;
    for (int i = 0; i < Volumetric3D::MaterialRefTypeString.size(); i++) {
        if (boost::iequals(params[3], Volumetric3D::MaterialRefTypeString[i]) ||
            boost::iequals(params[3], Volumetric3D::MaterialRefTypeStringShort[i])) {
            mReferenceType = Volumetric3D::MaterialRefType(i);
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_EMC::initialize || "
            "Unknown material reference type, type = " + params[3]);
    }
    
    try {
        int ipar = 4;
        Parameters::castValue(mFactor, params.at(ipar++), source);
        Parameters::castValue(mCartesian, params.at(ipar++), source);
        Parameters::castValue(mGeographic, params.at(ipar++), source);
        Parameters::castValue(mOneFilePerDepth, params.at(ipar++), source);
        Parameters::castValue(mVerticalDiscontinuities, params.at(ipar++), source);
        Parameters::castValue(mModelFlag, params.at(ipar++), source);
        Parameters::castValue(mModelFlagFactor, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    
    if (mCartesian && mGeographic) {
        throw std::runtime_error("Volumetric3D_EMC::initialize || "
            "Cartesian model has to be source-centered.");
    }
    
    if (!boost::iequals(mModelFlag, "none")) {
        if (mReferenceType == MaterialRefType::Absolute) {
            throw std::runtime_error("Volumetric3D_EMC::initialize || "
                "Imposing special model flag on an absolute model.");
        }
    }
    initialize();
}

void Volumetric3D_EMC::setSourceLocation(double srcLat, double srcLon, double srcDep) {
    mSrcLat = srcLat;
    mSrcLon = srcLon;
    mSrcDep = srcDep;
}

bool Volumetric3D_EMC::get3dPropertiesInternal(double r, double theta, double phi, double rElemCenter,
    std::vector<MaterialProperty> &properties, 
    std::vector<MaterialRefType> &refTypes,
    std::vector<double> &values, bool isFluid) const {

    // header
    properties = std::vector<MaterialProperty>(1, mMaterialProp);
    refTypes = std::vector<MaterialRefType>(1, mReferenceType);
    values = std::vector<double>(1, 0.);
    
    double dep, lat, lon;
    if (mCartesian) {
        RDCol3 rtpG;
        rtpG << r, theta, phi;
        RDCol3 xyz = Geodesy::Glob2Cartesian(rtpG, mSrcLat, mSrcLon, mSrcDep);
        lat = xyz(0);
        lon = xyz(1);
        dep = Geodesy::getROuter() - xyz(2);
        
        XMath::checkLimits(dep, 0., Geodesy::getROuter());
    } else {
        // to geocentric
        if (mGeographic) {
            // which radius to use?
            theta = pi / 2. - Geodesy::theta2Lat_d(theta, 0.) * degree;
        }
        
        // regularise
        dep = Geodesy::getROuter() - r;
        lat = 90. - theta / degree;
        lon = phi / degree;
        XMath::checkLimits(dep, 0., Geodesy::getROuter());
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
    
    // check center
    double dmin = mGridDep[0];
    double dmax = mGridDep[mGridDep.size() - 1];
    double dcenter = Geodesy::getROuter() - rElemCenter;
    if (dcenter < dmin || dcenter > dmax) {
        return false;
    }
    if (dep < dmin && dep > dmin * 0.999999) {
        dep = dmin;
    }
    if (dep > dmax && dep < dmax * 1.000001) {
        dep = dmax;
    }
    
    // interpolation
    double tol = 0.01;
    int ldep0, llat0, llon0, ldep1, llat1, llon1;
    double wdep0, wlat0, wlon0, wdep1, wlat1, wlon1;
    if (mVerticalDiscontinuities) {
        // use element center depth to locate layer
        XMath::interpLinear(dcenter, mGridDep, ldep0, wdep0);
        if (ldep0 < 0) {
            return false;
        }
        // use point depth to determine value
        wdep0 = 1. - 1. / (mGridDep(ldep0 + 1) - mGridDep(ldep0)) * (dep - mGridDep(ldep0));
    } else {
        XMath::interpLinearRobust(dep, mGridDep, ldep0, wdep0, tol);
    }
    XMath::interpLinearRobust(lat, mGridLat, llat0, wlat0, tol);
    XMath::interpLinearRobust(lon, mGridLon, llon0, wlon0, tol);    
    if (ldep0 < 0 || llat0 < 0 || llon0 < 0) {
        return false;
    }
    
    ldep1 = ldep0 + 1;
    llat1 = llat0 + 1;
    llon1 = llon0 + 1;

    std::vector<bool> use(8, true);
    for (const auto &model: mDiscontinuities) {
        model->handleDiscontinuities(r, theta, phi, mGeographic, mCartesian, rElemCenter, mGridLat(llat0), mGridLat(llat1), 
        mGridLon(llon0), mGridLon(llon1), mGridDep(ldep0), mGridDep(ldep1), use, 0);
    }
    
    if (std::none_of(use.begin(), use.end(), [] (const bool &b) {return b;} )) {
        throw std::runtime_error("Volumetric3D_EMC::get3dPropertiesInternal || "
            "Point on wrong side of discontinuity.");
    }
    
    RDColX weights = RDColX::Zero(8);
    if (std::any_of(use.begin(), use.end(), [] (const bool &b) {return !b;} )) {
        RDCol3 querycrd;
        querycrd << lat, lon, dep;
        if (!mCartesian) {
            RDCol3 rtpG;
            rtpG << r, theta, phi;
            querycrd = Geodesy::Glob2Cartesian(rtpG, mSrcLat, mSrcLon, mSrcDep);
            querycrd(2) = Geodesy::getROuter() - querycrd(2);
        }
        
        std::vector<RDCol3> crds;
        RDCol3 crd;
        crd << mGridLat(llat0) , mGridLon(llon0) , mGridDep(ldep0);
        crds.push_back(crd);
        crd << mGridLat(llat1) , mGridLon(llon0) , mGridDep(ldep0);
        crds.push_back(crd);
        crd << mGridLat(llat0) , mGridLon(llon1) , mGridDep(ldep0);
        crds.push_back(crd);
        crd << mGridLat(llat1) , mGridLon(llon1) , mGridDep(ldep0);
        crds.push_back(crd);
        crd << mGridLat(llat0) , mGridLon(llon0) , mGridDep(ldep1);
        crds.push_back(crd);
        crd << mGridLat(llat1) , mGridLon(llon0) , mGridDep(ldep1);
        crds.push_back(crd);
        crd << mGridLat(llat0) , mGridLon(llon1) , mGridDep(ldep1);
        crds.push_back(crd);
        crd << mGridLat(llat1) , mGridLon(llon1) , mGridDep(ldep1);
        crds.push_back(crd);
        
        for (int i = 0; i < 8; i++) {
            if (!use[i]) continue;
            RDCol3 crd = crds[i];
            if (!mCartesian) {
                crd(0) = Geodesy::lon2Phi(crd(0));
                crd(1) = Geodesy::lat2Theta_d(crd(1), crd(2));
                crd(2) = Geodesy::getROuter() - crd(2);
                crd = Geodesy::Glob2Cartesian(crd, mSrcLat, mSrcLon, mSrcDep);
                crd(2) = Geodesy::getROuter() - crd(2);
            }
            
            RDCol3 squaredDistance = querycrd - crd;
            squaredDistance = squaredDistance.array().pow(2);
            double dist = sqrt(squaredDistance.sum());
            if (dist < tinyDouble) {
                weights = RDColX::Zero(0);
                weights(i) = 1;
                break;
            } else {
                weights(i) = 1 / dist;
            }
        }
        double totalweight = weights.sum();
        weights = weights.array() / totalweight;
    } else {
        wdep1 = 1. - wdep0;
        wlat1 = 1. - wlat0;
        wlon1 = 1. - wlon0;
        
        weights(0) = wdep0 * wlat0 * wlon0;
        weights(1) = wdep0 * wlat1 * wlon0;
        weights(2) = wdep0 * wlat0 * wlon1;
        weights(3) = wdep0 * wlat1 * wlon1;
        weights(4) = wdep1 * wlat0 * wlon0;
        weights(5) = wdep1 * wlat1 * wlon0;
        weights(6) = wdep1 * wlat0 * wlon1;
        weights(7) = wdep1 * wlat1 * wlon1;
    }
    
    values[0] += mGridData[ldep0](llat0, llon0) * weights(0);
    values[0] += mGridData[ldep0](llat1, llon0) * weights(1);
    values[0] += mGridData[ldep0](llat0, llon1) * weights(2);
    values[0] += mGridData[ldep0](llat1, llon1) * weights(3);
    values[0] += mGridData[ldep1](llat0, llon0) * weights(4);
    values[0] += mGridData[ldep1](llat1, llon0) * weights(5);
    values[0] += mGridData[ldep1](llat0, llon1) * weights(6);
    values[0] += mGridData[ldep1](llat1, llon1) * weights(7);

    return true;
}

void Volumetric3D_EMC::findDiscontinuity(RDColX &lat, RDColX &lon, RDMatXX &depth, bool &cartesian, double val, 
    double minDepth, double maxDepth, std::string &compType, bool from_bottom) const {
    lat = mGridLat;
    lon = mGridLon;
    cartesian = mCartesian;

    int z0 = 0;
    int z1 = 0;
    
    if (mGridDep.maxCoeff() <= minDepth || mGridDep.minCoeff() >= maxDepth) {
        throw std::runtime_error("Volumetric3D_EMC::findDiscontinuity || "
            "Search range out of bounds of volumetric input model.");
    }
    
    while (mGridDep(++z0) < minDepth);
    z1 = z0;
    while (mGridDep(z1++) < maxDepth);
    
    int dir = (compType.compare("greater") == 0) ? 1 : -1;

    int dzi, zstart, zend;
    if (!from_bottom) {
        dzi = 1;
        zstart = z0;
        zend = z1;
    } else {
        dzi = -1;
        zstart = z1;
        zend = z0;
    }

    depth = RDMatXX::Zero(mGridLat.rows(), mGridLon.rows());
    for (int yi = 0; yi < mGridLon.rows(); yi++) {
        for (int xi = 0; xi < mGridLat.rows(); xi++) {
            bool found = false;
            for (int zi = z0; dzi * zi <= dzi * z1; zi += dzi) {
                if (dir * mGridData[zi](xi,yi) > dir * val) {
                    depth(xi,yi) = (mGridDep(zi - dzi) + mGridDep(zi)) / 2;
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error("Volumetric3D_EMC::findDiscontinuity || "
                    "Discontinuity no found at " + boost::lexical_cast<std::string>(mGridLat(xi)) + " , " + boost::lexical_cast<std::string>(mGridLon(yi)) + ".");
            }
        }
    }
}

std::string Volumetric3D_EMC::verbose() const {
    double dataMax = 0;
    double dataMin = 10e9;
    for (const auto &data: mGridData) {
        dataMax = std::max({dataMax, data.maxCoeff()});
        dataMin = std::min({dataMin, data.minCoeff()});
    }
    std::stringstream ss;
    std::string dim1str = mCartesian ? "XCoord" : "Latitude";
    std::string dim2str = mCartesian ? "YCoord" : "Longitude";
    std::string space1 = mCartesian ? "  " : "";
    std::string space2 = mCartesian ? "   " : "";
    ss << "\n======================= 3D Volumetric =======================" << std::endl;
    ss << "  Model Name           =   EMC" << std::endl;
    ss << "  Data File            =   " << mFileName << std::endl;
    ss << "  Variable Name        =   " << mVarName << std::endl;
    ss << "  Material Property    =   " << MaterialPropertyString[mMaterialProp] << std::endl;
    ss << "  Reference Type       =   " << MaterialRefTypeString[mReferenceType] << std::endl;
    ss << "  Values               =   [" << dataMin << ", " << dataMax << "]" << std::endl;
    ss << "  Num. Depths          =   " << mGridDep.size() << std::endl;
    ss << "  Num. " << dim1str << "s" << space1 << "       =   " << mGridLat.size() << std::endl;
    ss << "  Num. " << dim2str << "s" << space2 << "      =   " << mGridLon.size() << std::endl;
    ss << "  Depth Range          =   [" << mGridDep.minCoeff() << ", " << mGridDep.maxCoeff() << "]" << std::endl;
    ss << "  " << dim1str << " Range" << space1 << "       =   [" << mGridLat.minCoeff() << ", " << mGridLat.maxCoeff() << "]" << std::endl;
    ss << "  " << dim2str << " Range " << space2 << "     =   [" << mGridLon.minCoeff() << ", " << mGridLon.maxCoeff() << "]" << std::endl;
    ss << "  Factor               =   " << mFactor << std::endl;
    ss << "  Use Geographic       =   " << (mGeographic ? "YES" : "NO") << std::endl;
    ss << "  One File per Depth   =   " << (mOneFilePerDepth ? "YES" : "NO") << std::endl;
    if (!boost::iequals(mModelFlag, "none")) {
        ss << "  Special Model Flag   =   " << mModelFlag << std::endl;
        ss << "  Model Flag Factor    =   " << mModelFlagFactor << std::endl;
    }
    ss << "======================= 3D Volumetric =======================\n" << std::endl;
    return ss.str();
}