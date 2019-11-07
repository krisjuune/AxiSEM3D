// Volumetric3D_crust1.h
// created by Kuangdai on 16-May-2016
// crustal model CRUST 1.0
// http://igppweb.ucsd.edu/~gabi/crust1.html

#include "Volumetric3D_SEG.h"
#include <sstream>
#include <fstream>
#include "XMPI.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "AutoGeometricParams.h"
#include "ExodusModel.h"
#include <algorithm>
#include "eigenp.h"
#include <iostream>


std::vector<double> readDoubleArray(std::string fname) {
    std::fstream fs(fname, std::fstream::in);
    if (!fs) {
        throw std::runtime_error("Volumetric3D_SEG::readDoubleArray || "
                "Error opening parameter file: ||" + fname);
    }
    std::stringstream ss;
    ss << fs.rdbuf();
    fs.close();

    std::vector<double> array;
    double data;
    while (ss >> data) {
      array.push_back(data);
    }
    return array;
}

void readDoubleMat(std::string fname, int nz, RDMatXX &par) {
    std::fstream fs(fname, std::fstream::in);
    if (!fs) {
        throw std::runtime_error("Volumetric3D_SEG::readDoubleMat || "
                "Error opening parameter file: ||" + fname);
    }
    std::stringstream ss;
    ss << fs.rdbuf();
    fs.close();

    int i = 0;
    double data;
    while (ss >> data) {
        par(i % nz, int(i / nz)) = data;
        i++;
    }
}

void Volumetric3D_SEG::initialize(const std::vector<std::string> &params) {

    const std::string source = "Geometric3D_SEG::initialize";

    // initialize location
    if (params.size() > 0) {
        Parameters::castValue(mModelBathymetry, params[0], source);
    }
    if (params.size() > 1) {
        Parameters::castValue(mEqualizationDepth, params[1], source);
    }

    // read raw data
    std::string path = projectDirectory + "/src/3d_model/3d_volumetric/SEG/C3NA_data";
    if (XMPI::root()) {
        // std::fstream fsx, fsy, fsz;
        // fsx.open(path + "/SEG_C3NA.crdx", std::fstream::in);
        // fsy.open(path + "/SEG_C3NA.crdy", std::fstream::in);
        // fsz.open(path + "/SEG_C3NA.crdz", std::fstream::in);
        // if (!fsx || !fsy || !fsz) {
        //     throw std::runtime_error("Volumetric3D_SEG::initialize || "
        //         "Error opening SEG_C3NA.crd0 data files at directory: ||" + path);
        // }
        std::cout << "Volumetric3D_SEG::initialize || reading xyz.." << std::endl;
        std::vector<double> X = readDoubleArray(path + "/SEG_C3NA.crdx");
        std::vector<double> Y = readDoubleArray(path + "/SEG_C3NA.crdy");
        std::vector<double> Z = readDoubleArray(path + "/SEG_C3NA.crdz");

        mNx = X.size();
        mNy = Y.size();
        mNz = Z.size();

        // centering at source
        mX = RDColX::Zero(mNx);
        mY = RDColX::Zero(mNy);
        mZ = RDColX::Zero(mNz);

        mX(0) = X[0] - (X[0] + (X[mNx - 1] - X[0]) / 2);
        mY(0) = Y[0] - (Y[0] + (Y[mNy - 1] - Y[0]) / 2);
        mZ(0) = mRSurf - Z[0];

        for (int i = 1; i < mNx; i++) {
            mX(i) = X[i] - (X[0] + (X[mNx - 1] - X[0]) / 2);
            mModelHmin = std::min(mModelHmin, X[i] - X[i - 1]);
            mModelHmax = std::max(mModelHmax, X[i] - X[i - 1]);
        }
        for (int i = 1; i < mNy; i++) {mY(i) = Y[i] - (
            Y[0] + (Y[mNy - 1] - Y[0]) / 2);
            mModelHmin = std::min(mModelHmin, Y[i] - Y[i - 1]);
            mModelHmax = std::max(mModelHmax, Y[i] - Y[i - 1]);
        }
        for (int i = 1; i < mNz; i++) {
            mZ(i) = mRSurf - Z[i];
            mModelHmin = std::min(mModelHmin, Z[i] - Z[i - 1]);
            mModelHmax = std::max(mModelHmax, Z[i] - Z[i - 1]);
        }
    }

    XMPI::bcastEigen(mX);
    XMPI::bcastEigen(mY);
    XMPI::bcastEigen(mZ);
    mNx = mX.rows();
    mNy = mY.rows();
    mNz = mZ.rows();
    mVp = RDMatXX::Zero(mNz, mNx * mNy);
    mVs = RDMatXX::Zero(mNz, mNx * mNy);
    mRho = RDMatXX::Zero(mNz, mNx * mNy);

    int vp_mpirank = 0;
    int vs_mpirank = (XMPI::nproc() == 1) ? 0 : 1;
    int rho_mpirank = (XMPI::nproc() < 3) ? 0 : 2;

    if (XMPI::rank() == vp_mpirank) {
        std::cout << "Volumetric3D_SEG::initialize || reading vp.." << std::endl;
        readDoubleMat(path + "/SEG_C3NA.vp", mNz, mVp);
    }

    if (XMPI::rank() == vs_mpirank) {
        std::cout << "Volumetric3D_SEG::initialize || reading vs.." << std::endl;
        readDoubleMat(path + "/SEG_C3NA.vs", mNz, mVs);
    }

    if (XMPI::rank() == rho_mpirank) {
        std::cout << "Volumetric3D_SEG::initialize || reading rho.." << std::endl;
        readDoubleMat(path + "/SEG_C3NA.rho", mNz, mRho);
    }

    XMPI::cout << "Volumetric3D_SEG::initialize || broadcasting vp.." << XMPI::endl;
    XMPI::bcastEigen(mVp, vp_mpirank);
    XMPI::cout << "Volumetric3D_SEG::initialize || broadcasting vs.." << XMPI::endl;
    XMPI::bcastEigen(mVs, vs_mpirank);
    XMPI::cout << "Volumetric3D_SEG::initialize || broadcasting rho.." << XMPI::endl;
    XMPI::bcastEigen(mRho, rho_mpirank);

    XMPI::cout << "Volumetric3D_SEG::initialize || model complete." << XMPI::endl;
}

void Volumetric3D_SEG::modelBathymetry(std::vector<AutoGeometricParams *> &Vol2GeomModels) {
    if ((mVs.row(0).minCoeff() > 0) != (mMeshedOceanDepth == 0)) {
        throw std::runtime_error("Volumetric3D_SEG::modelBathymetry || Presence of ocean inconsistent between mesh and 3D model.");
    }
    if (mVs.row(0).minCoeff() > 0 || !mModelBathymetry) {
        mModelBathymetry = false;
        return;
    }

    AutoGeometricParams *v2g = new AutoGeometricParams(mMeshedOceanDepth,  mX, mY);

    RDCol2 Z_range;
    Z_range << mZ(0), std::max(mZ(mNz-1), mRSurf - mEqualizationDepth);
    v2g->setMeshLimits(Z_range);

    int zi_equal = 0;
    while (mZ(zi_equal) > Z_range(1)) {
        zi_equal++;
    }

    mZ_straightened = RDMatXX::Zero(mNz, mNx * mNy);
    double zmin = mRSurf;
    double zmax = 0;
    for (int yi = 0; yi < mNy; yi++) {
        for (int xi = 0; xi < mNx; xi++) {
            if (mVs(0, mNx * yi + xi) > 0) {
                throw std::runtime_error("Volumetric3D_SEG::modelBathymetry || Ocean is discontinuous.");
            }
            int zi_floor;
            bool found = false;
            for (int zi = 1; zi < mNz; zi++) {
                if (mVs(zi, mNx * yi + xi) > 0) {
                    found = true;
                    zi_floor = zi;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error("Volumetric3D_SEG::modelBathymetry || No ocean floor detected.");
            }

            v2g->setDepth(xi, yi, mZ(zi_floor));

            zmin = std::min(zmin, mZ(zi_floor));
            zmax = std::max(zmax, mZ(zi_floor));

            mZ_straightened.col(mNx * yi + xi) = mZ;
            double targetZ = mMeshedOceanDepth - (mZ(zi_floor - 1) - mZ(zi_floor)) / 2;
            double factor_above = (Z_range(0) - targetZ) / (Z_range(0) - mZ(zi_floor));
            double factor_below = (targetZ - Z_range(1)) / (mZ(zi_floor) - Z_range(1));
            for (int zi = 0; zi < zi_floor; zi++) {
                mZ_straightened(zi, mNx * yi + xi) = Z_range(0) + (mZ(zi) - Z_range(0)) * factor_above;
            }
            for (int zi = zi_floor; zi < zi_equal; zi++) {
                mZ_straightened(zi, mNx * yi + xi) = Z_range(1) + (mZ(zi) - Z_range(1)) * factor_below;
            }
        }
    }
    v2g->initialize();
    Vol2GeomModels.push_back(v2g);
}

void Volumetric3D_SEG::setupExodusModel(const ExodusModel *exModel) {
    mMeshedOceanDepth = exModel->getMeshedOcean();
    mHmin = exModel->getHmin();
    mHmax = exModel->getHmax();
}

bool Volumetric3D_SEG::get3dProperties(double r, double theta, double phi, double rElemCenter,
    std::vector<MaterialProperty> &properties,
    std::vector<MaterialRefType> &refTypes,
    std::vector<double> &values, bool isFluid) const {

    double vp_interp, vs_interp, rho_interp;

    RDCol3 rtp = RDCol3::Zero(3, 1);
    rtp(0) = r;
    rtp(1) = theta;
    rtp(2) = phi;
    RDCol3 xyz = Geodesy::toCartesian(rtp);

    // header
    properties.clear();
    properties.push_back(Volumetric3D::MaterialProperty::VP);
    properties.push_back(Volumetric3D::MaterialProperty::VS);
    properties.push_back(Volumetric3D::MaterialProperty::RHO);
    refTypes = std::vector<MaterialRefType>(3, Volumetric3D::MaterialRefType::Absolute);
    values = std::vector<double>(3, 0.);

    // allowing 1 cm error in grid points
    double tol = 0.01;

    int xi, yi, xj, yj;
    bool found = false;
    for (int i = 0; i < mNx; i++) {
        if (xyz(0) <= mX(i) + tol) {
            xj = i;
            xi = (xyz(0) >= mX(xj) - tol) ? xj : xj - 1;
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_SEG::get3dProperties || Mesh larger than model (x coordinate).");
    }

    found = false;
    for (int i = 0; i < mNy; i++) {
        if (xyz(1) <= mY(i) + tol) {
            yj = i;
            yi = (xyz(1) >= mY(yj) - tol) ? yj : yj - 1;
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("Volumetric3D_SEG::get3dProperties || Mesh larger than model (y coordinate).");
    }

    RDCol4 Z_cell_top = RDCol4::Zero(4,1);
    RDCol4 Z_cell_bot = RDCol4::Zero(4,1);
    RDCol4 vp_cell_top = RDCol4::Zero(4,1);
    RDCol4 vp_cell_bot = RDCol4::Zero(4,1);
    RDCol4 vs_cell_top = RDCol4::Zero(4,1);
    RDCol4 vs_cell_bot = RDCol4::Zero(4,1);
    RDCol4 rho_cell_top = RDCol4::Zero(4,1);
    RDCol4 rho_cell_bot = RDCol4::Zero(4,1);
    BCol4 exclude_top;
    BCol4 exclude_bot;
    exclude_top.fill(false);
    exclude_bot.fill(false);

    RDCol4 xyi;
    xyi << mNx * yi + xi, mNx * yi + xj, mNx * yj + xi, mNx * yj + xj;

    IRow4 zi, zj;

    if (!mModelBathymetry) {
        int z0, z1;

        found = false;
        for (int i = 0; i < mNz; i++) {
            if (xyz(2) >= mZ(i) - tol) {
               z1 = i;
               z0 = (xyz(2) <= mZ(xj)) ? i : i - 1;
               found = true;
               break;
            }
        }
        if (!found) {
           throw std::runtime_error("Volumetric3D_SEG::get3dProperties || Mesh larger than model (z coordinate).");
        }

        for (int corner = 0; corner < 4; corner++) {
            zi(corner) = z0;
            zj(corner) = z1;
            Z_cell_top(corner) = mZ(z0);
            Z_cell_bot(corner) = mZ(z1);
            vp_cell_top(corner) = mVp(z0, xyi(corner));
            vp_cell_bot(corner) = mVp(z1, xyi(corner));
            vs_cell_top(corner) = mVs(z0, xyi(corner));
            vs_cell_bot(corner) = mVs(z1, xyi(corner));
            rho_cell_top(corner) = mRho(z0, xyi(corner));
            rho_cell_bot(corner) = mRho(z1, xyi(corner));
            exclude_top(corner) = (isFluid != (mVs(z0, xyi(corner)) == 0));
            exclude_bot(corner) = (isFluid != (mVs(z1, xyi(corner)) == 0));
        }

    } else {
        int found = 0;
        for (int corner = 0; corner < 4; corner++) {
            for (int i = 0; i < mNz; i++) {
                if (xyz(2) >= mZ_straightened(i, xyi(corner)) - tol) {
                    zj(corner) = i;
                    Z_cell_bot(corner) = mZ_straightened(i, xyi(corner));
                    if (xyz(2) <= Z_cell_bot(corner)) {
                        Z_cell_top(corner) = Z_cell_bot(corner);
                        zi(corner) = i;
                    } else {
                        Z_cell_top(corner) = mZ_straightened(i - 1, xyi(corner));
                        zi(corner) = i - 1;
                    };

                    vp_cell_top(corner) = mVp(zi(corner), xyi(corner));
                    vs_cell_top(corner) = mVs(zi(corner), xyi(corner));
                    rho_cell_top(corner) = mRho(zi(corner), xyi(corner));
                    exclude_top(corner) = (isFluid != (mVs(zi(corner), xyi(corner)) == 0));

                    vp_cell_bot(corner) = mVp(zj(corner), xyi(corner));
                    vs_cell_bot(corner) = mVs(zj(corner), xyi(corner));
                    rho_cell_bot(corner) = mRho(zj(corner), xyi(corner));
                    exclude_bot(corner) = (isFluid != (mVs(zj(corner), xyi(corner)) == 0));

                    found += 1;
                    break;
                }
            }
        }
        if (found < 4) {
            throw std::runtime_error("Volumetric3D_SEG::get3dProperties || Mesh larger than model (z coordinate).");
        }
    }

    if (exclude_top.all() && exclude_bot.all()) {
        zi -= IRow4::Constant(4,1,1);
        zj -= IRow4::Constant(4,1,1);

        for (int corner = 0; corner < 4; corner++) {
            exclude_top(corner) = (isFluid != (mVs(zi(corner), xyi(corner)) == 0));
        }
        if (exclude_top.all()) {
            zi += IRow4::Constant(4,1,2);
            zj += IRow4::Constant(4,1,2);
            for (int corner = 0; corner < 4; corner++) {
                exclude_bot(corner) = (isFluid != (mVs(zj(corner), xyi(corner)) == 0));
            }
        }
        if (exclude_top.all() && exclude_bot.all()) {
            throw std::runtime_error("Volumetric3D_SEG::get3dProperties || Mesh solid-fluid boundary does not coincide with model.");
        }

        for (int corner = 0; corner < 4; corner++) {
            Z_cell_top(corner) = mZ_straightened(zi(corner), xyi(corner));
            vp_cell_top(corner) = mVp(zi(corner), xyi(corner));
            vs_cell_top(corner) = mVs(zi(corner), xyi(corner));
            rho_cell_top(corner) = mRho(zi(corner), xyi(corner));
            exclude_top(corner) = (isFluid != (mVs(zi(corner), xyi(corner)) == 0));

            Z_cell_bot(corner) = mZ_straightened(zj(corner), xyi(corner));
            vp_cell_bot(corner) = mVp(zj(corner), xyi(corner));
            vs_cell_bot(corner) = mVs(zj(corner), xyi(corner));
            rho_cell_bot(corner) = mRho(zj(corner), xyi(corner));
            exclude_bot(corner) = (isFluid != (mVs(zj(corner), xyi(corner)) == 0));
        }
    }

    RDCol4 d_top, d_bot;
    d_top << (xyz(0) - mX(xi)) * (xyz(0) - mX(xi)) + (xyz(1) - mY(xi)) * (xyz(1) - mY(xi)) + (xyz(2) - Z_cell_top(1)) * (xyz(2) - Z_cell_top(1)),
             (xyz(0) - mX(xj)) * (xyz(0) - mX(xj)) + (xyz(1) - mY(xi)) * (xyz(1) - mY(xi)) + (xyz(2) - Z_cell_top(2)) * (xyz(2) - Z_cell_top(2)),
             (xyz(0) - mX(xi)) * (xyz(0) - mX(xi)) + (xyz(1) - mY(xj)) * (xyz(1) - mY(xj)) + (xyz(2) - Z_cell_top(3)) * (xyz(2) - Z_cell_top(3)),
             (xyz(0) - mX(xj)) * (xyz(0) - mX(xj)) + (xyz(1) - mY(xj)) * (xyz(1) - mY(xj)) + (xyz(2) - Z_cell_top(4)) * (xyz(2) - Z_cell_top(4));
    d_bot << (xyz(0) - mX(xi)) * (xyz(0) - mX(xi)) + (xyz(1) - mY(xi)) * (xyz(1) - mY(xi)) + (xyz(2) - Z_cell_bot(1)) * (xyz(2) - Z_cell_bot(1)),
             (xyz(0) - mX(xj)) * (xyz(0) - mX(xj)) + (xyz(1) - mY(xi)) * (xyz(1) - mY(xi)) + (xyz(2) - Z_cell_bot(2)) * (xyz(2) - Z_cell_bot(2)),
             (xyz(0) - mX(xi)) * (xyz(0) - mX(xi)) + (xyz(1) - mY(xj)) * (xyz(1) - mY(xj)) + (xyz(2) - Z_cell_bot(3)) * (xyz(2) - Z_cell_bot(3)),
             (xyz(0) - mX(xj)) * (xyz(0) - mX(xj)) + (xyz(1) - mY(xj)) * (xyz(1) - mY(xj)) + (xyz(2) - Z_cell_bot(4)) * (xyz(2) - Z_cell_bot(4));

    vp_interp = multiv_interpolation(d_top, d_bot, vp_cell_top, vp_cell_bot, exclude_top, exclude_bot);
    vs_interp = multiv_interpolation(d_top, d_bot, vs_cell_top, vs_cell_bot, exclude_top, exclude_bot);
    rho_interp = multiv_interpolation(d_top, d_bot, rho_cell_top, rho_cell_bot, exclude_top, exclude_bot);

    if (vs_interp > vp_interp) {
            throw std::runtime_error("Volumetric3D_SEG::get3dProperties || Check 3D model: found vs larger than vp.");
    }

    values[0] = vp_interp;
    values[1] = vs_interp;
    values[2] = rho_interp;
    return true;
}

std::string Volumetric3D_SEG::verbose() const {
    std::stringstream ss;
    ss << "\n======================= 3D Volumetric ======================" << std::endl;
    ss << "  Model Name            =   SEG C3" << std::endl;
    ss << "  X Coords Range (km)   =   (" << round(mX(0)/1000) << ") - " << round(mX(mNx - 1)/1000) << std::endl;
    ss << "  Y Coords Range (km)   =   (" << round(mY(0)/1000) << ") - " << round(mX(mNy - 1)/1000) << std::endl;
    ss << "  Depth Range (km)      =   " << round((mRSurf - mZ(0))/1000) << " - " << round((mRSurf - mZ(mNz - 1))/1000) << std::endl;
    ss << "  Reference Type        =   Absolute" << std::endl;
    ss << "  Affected Properties   =   VP  (" << mVp.minCoeff() << " - " << mVp.maxCoeff() << ")" << std::endl;
    ss << "                            VS  (" << mVs.minCoeff() << " - " << mVs.maxCoeff() << ")" << std::endl;
    ss << "                            Rho (" << mRho.minCoeff() << " - " << mRho.maxCoeff() << ")" << std::endl;
    ss << "  Mesh vs Model Spacing =   min " << mHmin << " / " << mModelHmin << std::endl;
    ss << "                            max " << mHmax << " / " << mModelHmax << std::endl;
    std::string mb = mModelBathymetry ? "yes" : "no";
    ss << "  Ocean Bathymetry      =   " << mb << std::endl;
    ss << "======================= 3D Volumetric ======================\n" << std::endl;
    return ss.str();
}

double Volumetric3D_SEG::multiv_interpolation(const RDCol4 d_top, const RDCol4 d_bot,
    const RDCol4 val_top, const RDCol4 val_bot, const BCol4 exclude_top, const BCol4 exclude_bot) const {
    RDCol4 w_top = RDCol4::Zero(4, 1);
    RDCol4 w_bot = RDCol4::Zero(4, 1);

    for (int i = 0; i < 4; i++) {
       if (d_top(i) < tinyDouble && !exclude_top(i)) {
           return val_top(i);
       }
       if (d_bot(i) < tinyDouble && !exclude_bot(i)) {
           return val_bot(i);
       }
       w_top(i) = exclude_top(i) ? 0 : 1 / d_top(i);
       w_bot(i) = exclude_bot(i) ? 0 : 1 / d_bot(i);
    }

    double intval = 0;
    for (int i = 0; i < 4; i++) {
        intval += val_top(i) * w_top(i);
        intval += val_bot(i) * w_bot(i);
    }
    intval = intval / (w_top.sum() + w_bot.sum());
    return intval;
}
