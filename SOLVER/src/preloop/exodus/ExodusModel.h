// ExodusModel.h
// created by Kuangdai on 2-May-2016
// read Exodus mesh file

// Adapted from salvus project
// ExodusModel.h by Michael Afanasiev

#pragma once

#include <string>
#include <vector>
#include <map>
#include "eigenp.h"

class Parameters;
class AttParameters;

class ExodusModel {

public:
    ExodusModel(const std::string &fileName);
    void initialize();

    // general
    bool isIsotropic() const {return mElementalVariables.find("VP_0") != mElementalVariables.end();};
    bool hasAttenuation() const {return mGlobalVariables.find("nr_lin_solids") != mGlobalVariables.end();};
    int getNumQuads() const {return mConnectivity.rows();};
    int getNumQuadsInner() const {return mNumQuadsInner;};
    int getNumNodes() const {return mNodalS.rows();};
    double getROuter() const {
      try {
        return mGlobalVariables.at("radius");
      } catch (...) {return 6371000.0;}
    };
    bool isCartesian() const {return mGlobalRecords.at("crdsys") != "spherical";};
    bool hasABC() const {return (isCartesian() && mHasABC);}
    double getDistTolerance() const {return mDistTolerance;};
    double getHmax() const {return mHmax;};
    double getHmin() const {return mHmin;};
    RDCol2 getBoundaries() const {
      RDCol2 b(mNodalS.maxCoeff(),mNodalZ.minCoeff());
      return b;
    }
    double getMeshedOcean() const {return mMeshedOceanDepth;};

    // Node-wise
    double getNodalS(int nodeTag) const {return mNodalS(nodeTag);};
    double getNodalZ(int nodeTag) const {return mNodalZ(nodeTag);};
    double getAveGLLSpacing(int nodeTag) const {return mAveGLLSpacing(nodeTag);};

    // Quad-wise
    double getElementalVariables(const std::string &varName, int quadTag) const {
        return mElementalVariables.at(varName)(quadTag);
    };
    int getSideAxis(int quadTag) const {return mSideSets.at(mSSNameAxis)(quadTag);};
    int getSideSurface(int quadTag) const {return mSideSets.at(mSSNameSurface)(quadTag);};
    int getSideRightB(int quadTag) const {return mSideSets.at(mSSNameRightB)(quadTag);};
    int getSideLowerB(int quadTag) const {return mSideSets.at(mSSNameLowerB)(quadTag);};
    int getSideSolidFluid(int quadTag) const {
        if (mSideSets.find("solid_fluid_boundary") != mSideSets.end()) {
            return mSideSets.at("solid_fluid_boundary")(quadTag);
        }
        return -1;
    };
    const IMatX4 &getConnectivity() const {return mConnectivity;};
    IRow4 getConnectivity(int quadTag) const {return mConnectivity.row(quadTag);};
    IRow4 getVicinalAxis(int quadTag) const {return mVicinalAxis.row(quadTag);};

    int getNumAbsElements() const {return mN_ABC;};
    bool isABQuad(int quadTag) const {return (mABfield(quadTag, 1) > 0);};
    int getCopyTagAB(int quadTag) const {return mABfield(quadTag, 0);};
    int getABPosition(int quadTag) const {return mABfield(quadTag, 1);}; // 0 = edge of normal mesh; 1 = right ab boundary; 2 = lower ab boundary; 3 = ab corner
    double getABVmax() const {return mABC_Vmax;};

    std::string verbose() const;

    static void buildInparam(ExodusModel *&exModel, const Parameters &par,
        AttParameters *&attPar, int verbose);

    // double getR_CMB() const {return mR_CMB;};
    // double getR_ICB() const {return mR_ICB;};

private:

    void readRawData();
    void bcastRawData();
    void formStructured();
    void formAuxiliary();
    void AddAbsorbingBoundaryElements();

    // file name
    std::string mExodusFileName;

    // file properties
    std::string mExodusTitle;

    ///////////////////////////////////// raw data /////////////////////////////////////
    // global variables and records
    std::vector<std::string> mGlobalVariableNames;
    RDColX mGlobalVariableValues;
    std::vector<std::string> mGlobalRecordsRaw;

    // connectivity and coords
    IMatX4 mConnectivity;
    RDColX mNodalS, mNodalZ;

    // elemental variables
    std::vector<std::string> mElementalVariableNames;
    RDMatXX mElementalVariableValues;

    // side sets
    std::vector<std::string> mSideSetNames;
    IMatXX mSideSetValues;

    // ellipticity
    RDColX mEllipKnots;
    RDColX mEllipCoeffs;

    ///////////////////////////////////// structured /////////////////////////////////////
    // global variables and records
    std::map<std::string, double> mGlobalVariables;
    std::map<std::string, std::string> mGlobalRecords;

    // elemental variables
    std::map<std::string, RDColX> mElementalVariables;

    // side sets
    std::map<std::string, IColX> mSideSets;

    std::string mSSNameAxis = "t0";
    std::string mSSNameSurface = "r1";
    std::string mSSNameRightB = "t1";
    std::string mSSNameLowerB = "r0";

    ///////////////////////////////////// auxiliary /////////////////////////////////////
    // for Nr map
    RDColX mAveGLLSpacing;
    IMatX4 mVicinalAxis;
    double mDistTolerance;

    // for ABCs
    bool mHasABC;
    double mHmax, mHmin, mTSource;
    double mABC_Vmax = -1;
    int mN_ABC, mNumQuadsInner, mNumNodesInner;
    IMatX2 mABfield;

    // for automated oceantopography
    double mMeshedOceanDepth = 0;
};
