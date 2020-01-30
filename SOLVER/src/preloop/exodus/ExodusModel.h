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
    bool isIsotropic() const;
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
    bool hasStaceyABC() const {return (mHasStaceyABC);}
    bool hasSpongeBoundary() const {return (mHasSpongeBoundary);}
    bool hasMeshExtension() const {return (mHasMeshExtension);}
    bool hasParallelModelExtension() const {return (mHasModelExtension);}
    double getDistTolerance() const {return mDistTolerance;};
    double getHmax() const {return mHmax;};
    double getHmin() const {return mHmin;};
    RDCol2 getMeshBoundaryCoords() const {return mMeshEdgeCorner;};
    RDCol2 getSpongeStartCoords() const {return mInnerBoundaryCorner;};
    double getMeshedOcean() const {return mMeshedOceanDepth;};

    // Node-wise
    double getNodalS(int nodeTag) const {return mNodalS(nodeTag);};
    double getNodalZ(int nodeTag) const {return mNodalZ(nodeTag);};
    double getAveGLLSpacing(int nodeTag) const {return mAveGLLSpacing(nodeTag);};

    // Quad-wise
    double getElementalVariables(const std::string &varName, int quadTag) const;
    int getSideAxis(int quadTag) const {return mSideSets.at(mSSNameAxis)(quadTag);};
    int getSideSurface(int quadTag) const {return mSideSets.at(mSSNameSurface)(quadTag);};
    int getSideRightB(int quadTag) const {return mSideSets.at(mSSNameRightB)(quadTag);};
    int getSideLowerB(int quadTag) const {return mSideSets.at(mSSNameBottomB)(quadTag);};
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
    bool isSpongeQuad(int quadTag) const {return (mABfield(quadTag) > 0);}; // 0 = edge of normal mesh; 1 = right ab boundary; 2 = lower ab boundary; 3 = ab corner
    int getABPosition(int quadTag) const {return mABfield(quadTag);};
    
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
    void computeNodalRTheta();
    void addAbsorbingBoundaryElements();
    void setInnerBoundary();
    void setMeshEdge();
    void handleMeshRefinementForABC();
    void makeABField(RDCol2 inner);
    void extendRightSide(int nQuads, int nNodes, int N, double s);
    void extendBottomSide(int nQuads, int nNodes, int N, double s);
    
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
    RDColX mNodalS, mNodalZ, mNodalR, mNodalTheta;

    // elemental variables
    std::vector<std::string> mElementalVariableNames_all;
    // non-compact
    std::vector<std::string> mElementalVariableNames_elem;
    RDMatXX mElementalVariableValues_elem;
    // compact
    std::vector<std::string> mElementalVariableNames_axis;
    RDMatXX mElementalVariableValues_axis;
    RDColX mElementalVariableCoords_axis;
    
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
    // non-compact
    std::map<std::string, RDColX> mElementalVariables_elem;
    // compact
    std::map<std::string, RDColX> mElementalVariables_axis;
    
    // side sets
    std::map<std::string, IColX> mSideSets;

    std::string mSSNameAxis;
    std::string mSSNameSurface;
    std::string mSSNameRightB;
    std::string mSSNameBottomB;

    ///////////////////////////////////// auxiliary /////////////////////////////////////
    // for Nr map
    RDColX mAveGLLSpacing;
    IMatX4 mVicinalAxis;
    double mDistTolerance;

    // for ABCs
    bool mHasSpongeBoundary, mHasMeshExtension, mHasModelExtension, mHasStaceyABC;
    double mHmax, mHmin, mABC_Vmax, mTSource;
    int mNumQuadsInner, mNumNodesInner;
    RDCol2 mInnerBoundaryCorner, mMeshEdgeCorner;
    
    IColX mABfield;
    
    double mABCwidth = -1;
    int mN_maxWL_ABC = -1;
    int mN_ABC = -1;
    
    // for automated oceantopography
    double mMeshedOceanDepth = 0;
};
