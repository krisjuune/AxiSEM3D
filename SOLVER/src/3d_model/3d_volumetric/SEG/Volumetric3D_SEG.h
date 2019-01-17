// Volumetric3D_SEG.h
// created by Claudia on 14-Aug-2018
// reading in SEG models after conversion from sgy

#pragma once
#include "Volumetric3D.h"
#include "eigenp.h"

class AutoGeometricParams;
class ExodusModel;

class Volumetric3D_SEG: public Volumetric3D {

public:
    void initialize(const std::vector<std::string> &params);
    void setupExodusModel(const ExodusModel *exModel);
    void modelBathymetry(std::vector<AutoGeometricParams *> &Vol2GeomModels);

    bool get3dProperties(double r, double theta, double phi, double rElemCenter,
        std::vector<MaterialProperty> &properties,
        std::vector<MaterialRefType> &refTypes,
        std::vector<double> &values, bool isFluid) const;

    std::string verbose() const;
    //bool makeFluid3D() const {return true;};

private:
    double multiv_interpolation(const RDCol4 d_top, const RDCol4 d_bot,
        const RDCol4 val_top, const RDCol4 val_bot,
        const BCol4 exclude_bot, const BCol4 exclude_top) const;

private:
    // radus of reference sphere
    double mRSurf = 6371000.0;

    bool mModelBathymetry = true;
    double mMeshedOceanDepth = 0;
    double mEqualizationDepth = mRSurf;

    double mModelHmin = mRSurf;
    double mModelHmax = 0;
    double mHmin, mHmax;
    mutable double mVpmins = 8000;
    mutable double mVpmaxs = 0;
    mutable double mVsmins = 8000;
    mutable double mVsmaxs = 0;
    mutable double mVpminf = 8000;
    mutable double mVpmaxf = 0;
    mutable double mVsminf = 8000;
    mutable double mVsmaxf = 0;

    int mNx, mNy, mNz;
    RDColX mX, mY, mZ;
    RDMatXX mVp, mVs, mRho, mZ_straightened;
};
