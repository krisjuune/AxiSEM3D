// ExodusModel.cpp
// created by Kuangdai on 2-May-2016 
// read Exodus mesh file

// Adapted from salvus project 
// ExodusModel.cpp by Michael Afanasiev


#include "ExodusModel.h"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cfloat>

#include "XMPI.h"
#include "Mapping.h"
#include <boost/algorithm/string.hpp>

#include "AttBuilder.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "MultilevelTimer.h"

#include "NetCDF_Reader.h"

ExodusModel::ExodusModel(const std::string &fileName): mExodusFileName(fileName) {
    std::vector<std::string> substrs = Parameters::splitString(mExodusFileName, "/");
    mExodusTitle = substrs[substrs.size() - 1];
    boost::trim_if(mExodusTitle, boost::is_any_of("\t "));
}

void ExodusModel::initialize() {
    MultilevelTimer::begin("Read Exodus", 1);
    if (XMPI::root()) {
        readRawData();
    }
    MultilevelTimer::end("Read Exodus", 1);
    
    MultilevelTimer::begin("Bcast Exodus", 1);
    bcastRawData();
    MultilevelTimer::end("Bcast Exodus", 1);
    
    MultilevelTimer::begin("Process Exodus", 1);
    formStructured();
    formAuxiliary();
    MultilevelTimer::end("Process Exodus", 1);
}

void ExodusModel::readRawData() {
    // open file
    NetCDF_Reader reader;
    reader.open(mExodusFileName);
    RDMatXX dbuffer;
    
    // global
    reader.readString("name_glo_var", mGlobalVariableNames);
    reader.read1D("vals_glo_var", mGlobalVariableValues);
    reader.readString("info_records", mGlobalRecordsRaw);
    
    // connectivity and coords
    reader.read2D("connect1", mConnectivity);
    mConnectivity.array() -= 1;
    reader.read1D("coordx", mNodalS);
    reader.read1D("coordy", mNodalZ);
        
    // elemental variables
    reader.readString("name_elem_var", mElementalVariableNames);
    mElementalVariableValues = RDMatXX::Zero(getNumQuads(), mElementalVariableNames.size());
    for (int i = 0; i < mElementalVariableNames.size(); i++) {
        std::stringstream ss;
        ss << "vals_elem_var" << i + 1 << "eb1";
        reader.read2D(ss.str(), dbuffer);
        mElementalVariableValues.col(i) = dbuffer.transpose();
    }    
    
    // side sets
    reader.readString("ss_names", mSideSetNames);
    mSideSetValues = IMatXX::Zero(getNumQuads(), mSideSetNames.size());
    for (int i = 0; i < mSideSetNames.size(); i++) {
        std::stringstream sse, sss;
        IColX elems, sides;
        // elem
        sse << "elem_ss" << i + 1;
        reader.read1D(sse.str(), elems);
        // side
        sss << "side_ss" << i + 1;
        reader.read1D(sss.str(), sides);
        // processing
        IColX values = IColX::Constant(getNumQuads(), -1);
        for (int j = 0; j < sides.rows(); j++) {
            values(elems(j) - 1) = sides(j) - 1;
        }
        mSideSetValues.col(i) = values;
    }
    
    // ellipticity
    try {
      reader.read2D("ellipticity", dbuffer);
      mEllipKnots = dbuffer.row(0).transpose();
      mEllipCoeffs = dbuffer.row(1).transpose();
    } catch (...) {std::cout << "no ellipticity argument" << std::endl;}

    // close file
    reader.close();
}

void ExodusModel::bcastRawData() {
    XMPI::bcast(mGlobalVariableNames);
    XMPI::bcastEigen(mGlobalVariableValues);
    XMPI::bcast(mGlobalRecordsRaw);
    
    XMPI::bcastEigen(mConnectivity);
    XMPI::bcastEigen(mNodalS);
    XMPI::bcastEigen(mNodalZ);
    
    XMPI::bcast(mElementalVariableNames);
    XMPI::bcastEigen(mElementalVariableValues);
    
    XMPI::bcast(mSideSetNames);
    XMPI::bcastEigen(mSideSetValues);
    
    XMPI::bcastEigen(mEllipKnots);
    XMPI::bcastEigen(mEllipCoeffs);
}

void ExodusModel::formStructured() {
    // global variables
    for (int i = 0; i < mGlobalVariableNames.size(); i++) {
        std::string varName = mGlobalVariableNames[i];
        if (varName == "dt") {
            varName = "dt (nPol = 1)";
        }
        mGlobalVariables.insert(std::pair<std::string, double>(varName, mGlobalVariableValues(i)));
    }
    
    // global records
    std::vector<std::string> included = {"crdsys", "model"};
    for (int i = 0; i < mGlobalRecordsRaw.size(); i++) {
        std::vector<std::string> substrs = Parameters::splitString(mGlobalRecordsRaw[i], "=");
        if (std::find(included.begin(), included.end(), boost::trim_copy(substrs[0])) != included.end()) {
            mGlobalRecords.insert(std::pair<std::string, std::string>(
                boost::trim_copy(substrs[0]), boost::trim_copy(substrs[1])));
        }
    }
    
    // elemental variables
    for (int i = 0; i < mElementalVariableNames.size(); i++) {
        mElementalVariables.insert(std::pair<std::string, RDColX>(mElementalVariableNames[i], 
            mElementalVariableValues.col(i)));
    }
    
    // side sets
    for (int i = 0; i < mSideSetNames.size(); i++) {
        mSideSets.insert(std::pair<std::string, IColX>(mSideSetNames[i], 
            mSideSetValues.col(i)));
    }
    
    // name of axis and surface sets
    mSSNameAxis = isCartesian() ? "x0" : "t0";
    mSSNameSurface = isCartesian() ? "y1" : "r1";
    mSSNameRightB = isCartesian() ? "x1" : "t1";
    mSSNameLowerB = isCartesian() ? "y0" : "r0";

    // NOTE: we temporarily treat Cartesian meshes as special cases of spherical meshes
    //       by means of moving it to the "north pole". The introduced global 
    //       curvature should be ignorable, or the problem itself is ill-defined
    //       as a local problem.
    if (isCartesian()) {
        double R_EARTH = getROuter();
        double maxz = mNodalZ.maxCoeff();
        mNodalZ.array() += R_EARTH - maxz; 
    }
}

void ExodusModel::formAuxiliary() {
    MultilevelTimer::begin("Process Exodus DistTol", 2);
    // distance tolerance and maximum side length of elements
    mHmax = 0;
    double distTol = DBL_MAX;
    for (int i = 0; i < getNumQuads(); i++) {
        if (i % XMPI::nproc() != XMPI::rank()) {
            continue;
        }
        double s0 = mNodalS(mConnectivity(i, 0));
        double z0 = mNodalZ(mConnectivity(i, 0));
        double s1 = mNodalS(mConnectivity(i, 1));
        double z1 = mNodalZ(mConnectivity(i, 1));
        double s2 = mNodalS(mConnectivity(i, 2));
        double z2 = mNodalZ(mConnectivity(i, 2));
        double s3 = mNodalS(mConnectivity(i, 3));
        double z3 = mNodalZ(mConnectivity(i, 3));
        double dist0 = sqrt((s0 - s1) * (s0 - s1) + (z0 - z1) * (z0 - z1)) / 1000.;
        double dist1 = sqrt((s1 - s2) * (s1 - s2) + (z1 - z2) * (z1 - z2)) / 1000.;
        double dist2 = sqrt((s2 - s3) * (s2 - s3) + (z2 - z3) * (z2 - z3)) / 1000.;
        double dist3 = sqrt((s3 - s0) * (s3 - s0) + (z3 - z0) * (z3 - z0)) / 1000.;
        distTol = std::min({dist0, dist1, dist2, dist3, distTol});
        mHmax = std::max({dist0, dist1, dist2, dist3, mHmax});
    }
    mHmax = 1000 * mHmax;
    mDistTolerance = XMPI::min(distTol);
    MultilevelTimer::end("Process Exodus DistTol", 2);
    
    // rotate nodes of axial elements such that side 3 is on axis, side 1 is the right boundary and side 0 is the lower boundary
    MultilevelTimer::begin("Process Exodus Axis", 2);
    for (int axialQuad = 0; axialQuad < getNumQuads(); axialQuad++) {

        int shift=0;

        int axialSide = getSideAxis(axialQuad);
        int rightSide = getSideRightB(axialQuad);
        int bottomSide = getSideLowerB(axialQuad);

        if (hasABC() & rightSide != 1 & rightSide != -1) {shift = rightSide - 1;}
        if (hasABC() & bottomSide != 0 & bottomSide != -1) {shift = bottomSide;}
        if (axialSide != 3 & axialSide != -1) {shift = axialSide - 3;}

        if (shift==0) {
            continue;
        }

        // connectivity
        IRow4 con = mConnectivity.row(axialQuad);
        for (int j = 0; j < 4; j++) {
            mConnectivity(axialQuad, j) = con(Mapping::period0123(j + shift));
        }

        // elemental fields
        for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) {
            std::string vname = it->first;
            if (vname.substr(vname.length() - 2, 2) == std::string("_0")) {
            vname = vname.substr(0, vname.length() - 2);
                RDRow4 v_old;
                v_old(0) = mElementalVariables.at(vname + "_0")(axialQuad);
                v_old(1) = mElementalVariables.at(vname + "_1")(axialQuad);
                v_old(2) = mElementalVariables.at(vname + "_2")(axialQuad);
                v_old(3) = mElementalVariables.at(vname + "_3")(axialQuad);
                mElementalVariables.at(vname + "_0")(axialQuad) = v_old(Mapping::period0123(0 + shift));
                mElementalVariables.at(vname + "_1")(axialQuad) = v_old(Mapping::period0123(1 + shift));
                mElementalVariables.at(vname + "_2")(axialQuad) = v_old(Mapping::period0123(2 + shift));
                mElementalVariables.at(vname + "_3")(axialQuad) = v_old(Mapping::period0123(3 + shift));
            }
        }

        // side sets
        for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
            if (it->second(axialQuad) != -1) {
                it->second(axialQuad) = Mapping::period0123(it->second(axialQuad) - shift);
            }
        }
    // done
    // mSideSets.at(mSSNameAxis)(axialQuad) = 3;
    }
    MultilevelTimer::end("Process Exodus Axis", 2);

    MultilevelTimer::begin("Process Absorbing Boundaries", 2);
    // extending mesh to incorporate absorbing boundaries

    // Absorbing Boundaries
    for (int Quad = 0; Quad < getNumQuads(); Quad++) {
        mABfield.push_back({-1,-1}); // populate ABfield
    }
    if (hasABC()) {
        AddAbsorbingBoundaryElements();
    } else {
        mNumQuadsInner = getNumQuads();
        mNumNodesInner = getNumNodes();
    }
    MultilevelTimer::end("Process Absorbing Boundaries", 2);

    // average gll spacing
    MultilevelTimer::begin("Process Exodus GLL-Spacing", 2);
    std::vector<std::vector<int>> refElem(getNumNodes(), std::vector<int>());
    for (int i = 0; i < getNumQuads(); i++) {
        refElem[mConnectivity(i, 0)].push_back(i);
        refElem[mConnectivity(i, 1)].push_back(i);
        refElem[mConnectivity(i, 2)].push_back(i);
        refElem[mConnectivity(i, 3)].push_back(i);
    }
    mAveGLLSpacing = RDColX::Zero(getNumNodes());
    for (int i = 0; i < getNumNodes(); i++) {
        if (i % XMPI::nproc() != XMPI::rank()) {
            continue;
        }
        for (int j = 0; j < refElem[i].size(); j++) {
            int ielem = refElem[i][j];
            double s0 = mNodalS(mConnectivity(ielem, 0));
            double z0 = mNodalZ(mConnectivity(ielem, 0));
            double s1 = mNodalS(mConnectivity(ielem, 1));
            double z1 = mNodalZ(mConnectivity(ielem, 1));
            double s2 = mNodalS(mConnectivity(ielem, 2));
            double z2 = mNodalZ(mConnectivity(ielem, 2));
            double s3 = mNodalS(mConnectivity(ielem, 3));
            double z3 = mNodalZ(mConnectivity(ielem, 3));
            double dist0 = sqrt((s0 - s1) * (s0 - s1) + (z0 - z1) * (z0 - z1));
            double dist1 = sqrt((s1 - s2) * (s1 - s2) + (z1 - z2) * (z1 - z2));
            double dist2 = sqrt((s2 - s3) * (s2 - s3) + (z2 - z3) * (z2 - z3));
            double dist3 = sqrt((s3 - s0) * (s3 - s0) + (z3 - z0) * (z3 - z0));
            mAveGLLSpacing(i) += (dist0 + dist1 + dist2 + dist3) / 4. / nPol / refElem[i].size();
        }
    } 
    XMPI::sumEigenDouble(mAveGLLSpacing);
    MultilevelTimer::end("Process Exodus GLL-Spacing", 2);

    // locate Solid-fluid Boundaries
    MultilevelTimer::begin("Process Exodus Solid-Fluid Boundary", 2);
    // find nodes which are part of both solid and fluid quads
    std::vector<bool> SFNode;
    for (int i = 0; i < getNumNodes(); i++) {
        int SumFluid = 0;
        for (int j = 0; j < refElem[i].size(); j++) {
            SumFluid += getElementalVariables("fluid",refElem[i][j]);
        }
        SFNode.push_back(0 < SumFluid && SumFluid < refElem[i].size());
    }
    // identify side of quads which is on SF boundary
    IColX SideSets_SF(getNumQuads());
    for (int iQuad = 0; iQuad < getNumQuads(); iQuad++) {
        SideSets_SF(iQuad) = -1;
        for (int node0 = 0; node0 < 4; node0++) {
            int node1 = node0 + 1 < 4 ? node0 + 1 : 0;
            if (SFNode[mConnectivity(iQuad, node0)] && SFNode[mConnectivity(iQuad, node1)]) {
                SideSets_SF(iQuad) = node0;
            }
        }
    }
    mSideSets.insert(std::pair<std::string, IColX>("solid_fluid_boundary", SideSets_SF));
    MultilevelTimer::end("Process Exodus Solid-Fluid Boundary", 2);

    // find elements that are not axial but neighboring axial elements
    MultilevelTimer::begin("Process Exodus Vicinal", 2);
    mVicinalAxis = IMatX4::Constant(getNumQuads(), 4, -1);
    // first find near-axis nodes and axial quads
    std::vector<bool> nodeNearAxis(getNumNodes(), false);
    std::vector<bool> quadOnAxis(getNumQuads(), false);
    for (int axialQuad = 0; axialQuad < getNumQuads(); axialQuad++) {
        int axialSide = getSideAxis(axialQuad);
        if (axialSide == -1) {
            continue;
        }
        nodeNearAxis[mConnectivity(axialQuad, 0)] = true;
        nodeNearAxis[mConnectivity(axialQuad, 1)] = true;
        nodeNearAxis[mConnectivity(axialQuad, 2)] = true;
        nodeNearAxis[mConnectivity(axialQuad, 3)] = true;
        quadOnAxis[axialQuad] = true;
    }
    // loop over quads
    for (int iquad = 0; iquad < getNumQuads(); iquad++) {
        if (quadOnAxis[iquad]) continue;
        for (int j = 0; j < 4; j++) {
            int nTag = mConnectivity(iquad, j);
            if (nodeNearAxis[nTag]) mVicinalAxis(iquad, j) = j;
        }
    }
    MultilevelTimer::end("Process Exodus Vicinal", 2);
    
    // check if ocean presents in mesh
    MultilevelTimer::begin("Process Exodus Check Ocean", 2);
    std::string strVs = isIsotropic() ? "VS_0" : "VSV_0";
    for (int iQuad = 0; iQuad < getNumQuads(); iQuad++) {
        if (iQuad % XMPI::nproc() != XMPI::rank()) {
            continue;
        }
        if (getSideSurface(iQuad) == -1) {
            continue;
        }
        double vs = mElementalVariables.at(strVs)(iQuad);
        if (vs < tinyDouble) {
            throw std::runtime_error("ExodusModel::finishReading || "
                "Ocean is detected in mesh. By far, realistic ocean is not implemented in AxiSEM3D. ||"
                "Use non-ocean models in the Mesher and add ocean load in inparam.basic.");
        }
    }
    MultilevelTimer::end("Process Exodus Check Ocean", 2);
}

void ExodusModel::AddAbsorbingBoundaryElements() {

    int node0, node1, node2, node3;
    int Side, nNewNodes_right = 0, nNewNodes_bottom = 0;
    int iABCcopy = 0;
    bool set1, set2;
    RDCol2 node0coords, node1coords, node2coords, node3coords;
    std::vector<RDCol2> NewNodeCoords_right, NewNodeCoords_bottom;
    IRow4 newABCelem;
    std::vector<double> mapABCtoGlobal;

    int nQuads = getNumQuads();
    int nQuads_ini = nQuads;
    int nNodes = getNumNodes();

    // for verbose output to distinguish between inner mesh and boundaries
    mNumQuadsInner = nQuads_ini;
    mNumNodesInner = nNodes;

    // add N_ABC boundary elements to the right of each quad at x1
    for (int rightB_Quad = 0; rightB_Quad < nQuads_ini; rightB_Quad++) {

        if (getSideRightB(rightB_Quad) == 1) {
            mSideSets.at(mSSNameRightB)(rightB_Quad) = -1; // after adding new elements selected quad is no longer at the edge of the mesh

            mABfield[rightB_Quad] = {0,++iABCcopy};

            // find nodes at boundary of inner mesh
            node0 = mConnectivity(rightB_Quad,1);
            node3 = mConnectivity(rightB_Quad,2);

            // specified number of boundary elements added to the right side of the selected quad
            for (int n = 1; n <= mN_ABC; n++) {
                set1 = 0;
                set2 = 0;

                // coordinates for outer nodes of the new quad
                node1coords[0] = mNodalS(node0) + mHmax;
                node1coords[1] = mNodalZ(node0);
                node2coords[0] = mNodalS(node3) + mHmax;
                node2coords[1] = mNodalZ(node3);

                // check if outer nodes had been added previously
                for (int i = 0; i < nNewNodes_right; i++) {
                    if (NewNodeCoords_right[i] == node1coords) {
                        node1 = mapABCtoGlobal[i];
                        set1 = 1;
                    }
                    if (NewNodeCoords_right[i] == node2coords) {
                        node2 = mapABCtoGlobal[i];
                        set2 = 1;
                    }
                }

                // adding new nodes if required
                if (!set1) {
                    node1 = nNodes;
                    NewNodeCoords_right.push_back(node1coords);
                    mapABCtoGlobal.push_back(nNodes);
                    mNodalS.conservativeResize(nNodes + 1, 1);
                    mNodalS[nNodes] = node1coords[0];
                    mNodalZ.conservativeResize(nNodes + 1, 1);
                    mNodalZ[nNodes] = node1coords[1];
                    nNodes++;
                    nNewNodes_right++;
                }
                if (!set2) {
                    node2 = nNodes;
                    NewNodeCoords_right.push_back(node2coords);
                    mapABCtoGlobal.push_back(nNodes);
                    mNodalS.conservativeResize(nNodes + 1, 1);
                    mNodalS[nNodes] = node2coords[0];
                    mNodalZ.conservativeResize(nNodes + 1, 1);
                    mNodalZ[nNodes] = node2coords[1];
                    nNodes++;
                    nNewNodes_right++;
                }

                // add new quad to connectivity map
                newABCelem = {node0, node1, node2, node3};
                mConnectivity.conservativeResize(nQuads + 1, 4);
                mConnectivity.row(nQuads) = newABCelem;

                // update list of elemental variables and side sets (copies from nearest quad in inner mesh)
                for (int i = 0; i < mElementalVariableNames.size(); i++) {
                    std::string vname = mElementalVariableNames[i];
                    RDColX temp = mElementalVariables.at(vname);
                    temp.conservativeResize(nQuads + 1,1);
                    // handling gradients near the boundary -> only constant values from nodes at edge are passed on
                    if (vname.substr(vname.length() - 2, 2) == std::string("_0")) {
                        temp(nQuads) = mElementalVariables.at(vname.substr(0, vname.length() - 2) + "_1")(rightB_Quad);
                    } else if (vname.substr(vname.length() - 2, 2) == std::string("_3")) {
                        temp(nQuads) = mElementalVariables.at(vname.substr(0, vname.length() - 2) + "_2")(rightB_Quad);
                    } else {
                        temp(nQuads) = mElementalVariables.at(vname)(rightB_Quad);
                    }
                    mElementalVariables.at(vname) = temp;
                }
                for (int i = 0; i < mSideSetNames.size(); i++) {
                    IColX temp = mSideSets.at(mSideSetNames[i]);
                    temp.conservativeResize(nQuads + 1,1);
                    temp(nQuads) = mSideSets.at(mSideSetNames[i])(rightB_Quad);
                    mSideSets.at(mSideSetNames[i]) = temp;
                }

                // updating inner nodes for next iteration
                node0 = node1;
                node3 = node2;

                mABfield.push_back({1,iABCcopy});
                nQuads++;
            }

            mSideSets.at(mSSNameRightB)(nQuads - 1) = 1; // last boundary element is placed at edge of the mesh
        }
    }

    // update quad number to include right boundary zone
    nQuads_ini = nQuads;
    mapABCtoGlobal.clear();

    // add boundary at y0 (lower bound)
    // same principle - see comments of section above
    for (int lowerB_Quad = 0; lowerB_Quad < nQuads_ini; lowerB_Quad++) {

        if (getSideLowerB(lowerB_Quad) == 0) {
            mSideSets.at(mSSNameLowerB)(lowerB_Quad) = -1;

            if (mABfield[lowerB_Quad](0) == -1) {
                mABfield[lowerB_Quad] = {0,++iABCcopy};
            }

            node2 = mConnectivity(lowerB_Quad,1);
            node3 = mConnectivity(lowerB_Quad,0);

            for (int n = 1; n <= mN_ABC; n++) {
                set1 = 0;
                set2 = 0;

                node0coords[0] = mNodalS(node3);
                node0coords[1] = mNodalZ(node3) - mHmax;
                node1coords[0] = mNodalS(node2);
                node1coords[1] = mNodalZ(node2) - mHmax;

                for (int i = 0; i < nNewNodes_bottom; i++){
                    if (NewNodeCoords_bottom[i] == node0coords) {
                        node0 = mapABCtoGlobal[i];
                        set1 = 1;
                    }
                    if (NewNodeCoords_bottom[i] == node1coords) {
                        node1 = mapABCtoGlobal[i];
                        set2 = 1;
                    }
                }

                if (!set1) {
                    node0 = nNodes;
                    NewNodeCoords_bottom.push_back(node0coords);
                    mapABCtoGlobal.push_back(nNodes);
                    mNodalS.conservativeResize(nNodes + 1, 1);
                    mNodalS[nNodes] = node0coords[0];
                    mNodalZ.conservativeResize(nNodes + 1, 1);
                    mNodalZ[nNodes] = node0coords[1];
                    nNodes++;
                    nNewNodes_bottom++;
                }
                if (!set2) {
                    node1 = nNodes;
                    NewNodeCoords_bottom.push_back(node1coords);
                    mapABCtoGlobal.push_back(nNodes);
                    mNodalS.conservativeResize(nNodes + 1, 1);
                    mNodalS[nNodes] = node1coords[0];
                    mNodalZ.conservativeResize(nNodes + 1, 1);
                    mNodalZ[nNodes] = node1coords[1];
                    nNodes++;
                    nNewNodes_bottom++;
                }

                newABCelem = {node0, node1, node2, node3};
                mConnectivity.conservativeResize(nQuads + 1, 4);
                mConnectivity.row(nQuads) = newABCelem;

                for (int i = 0; i < mElementalVariableNames.size(); i++) {
                    std::string vname = mElementalVariableNames[i];
                    RDColX temp = mElementalVariables.at(vname);
                    temp.conservativeResize(nQuads + 1,1);
                    // handling gradients near the boundary -> only constant values from nodes at edge are passed on
                    if (vname.substr(vname.length() - 2, 2) == std::string("_2")) {
                        temp(nQuads) = mElementalVariables.at(vname.substr(0, vname.length() - 2) + "_1")(lowerB_Quad);
                    } else if (vname.substr(vname.length() - 2, 2) == std::string("_3")) {
                        temp(nQuads) = mElementalVariables.at(vname.substr(0, vname.length() - 2) + "_0")(lowerB_Quad);
                    } else {
                        temp(nQuads) = mElementalVariables.at(vname)(lowerB_Quad);
                    }
                    mElementalVariables.at(vname) = temp;
                }
                for (int i = 0; i < mSideSetNames.size(); i++) {
                    IColX temp = mSideSets.at(mSideSetNames[i]);
                    temp.conservativeResize(nQuads + 1,1);
                    temp(nQuads) = mSideSets.at(mSideSetNames[i])(lowerB_Quad);
                    mSideSets.at(mSideSetNames[i]) = temp;
                }
                mSideSets.at(mSSNameLowerB)(mSideSets.at(mSSNameLowerB).rows() - 1) = (n == mN_ABC) ? 0 : -1;

                node2 = node1;
                node3 = node0;

                mABfield.push_back({1,iABCcopy-1});
                nQuads++;
            }

            mSideSets.at(mSSNameLowerB)(nQuads - 1) = 0;
        }
    }
}


std::string ExodusModel::verbose() const {
    std::stringstream ss;
    ss << "\n======================= Exodus Model =======================" << std::endl;
    ss << "  Overview__________________________________________________" << std::endl;
    ss << "    Exodus Title      =   " << mExodusTitle << std::endl;
    ss << "    Mesh CS Type      =   " << (isCartesian() ? "Cartesian" : "Spherical") << std::endl;
    ss << "    Number of Nodes   =   " << getNumNodes() << std::endl;
    ss << "    Number of Quads   =   " << getNumQuads() << std::endl;
    ss << "  Global Variables__________________________________________" << std::endl;
    int widthname = -1;
    for (auto it = mGlobalVariables.begin(); it != mGlobalVariables.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mGlobalRecords.begin(); it != mGlobalRecords.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));    
    }
    for (auto it = mGlobalVariables.begin(); it != mGlobalVariables.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << "   =   " << it->second << std::endl;
    }
    for (auto it = mGlobalRecords.begin(); it != mGlobalRecords.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << "   =   " << it->second << std::endl;
    }
    ss << "  Connectivity______________________________________________" << std::endl;
    int width = (int)std::log10(std::max(getNumQuads(), getNumNodes())) + 1;
    ss << "    " << std::setw(width) << 0 << ": ";
    for (int j = 0; j < 4; j++) {
        ss << std::setw(width) << mConnectivity(0, j) << " ";
    }
    ss << std::endl << "    " << std::setw(width) << "..." << std::endl;
    ss << "    " << std::setw(width) << getNumQuads() - 1 << ": ";
    for (int j = 0; j < 4; j++) {
        ss << std::setw(width) << mConnectivity(getNumQuads() - 1, j) << " ";
    }
    ss << std::endl;
    ss << "  Coordinates_______________________________________________" << std::endl;
    ss << "    " << std::setw(width) << 0 << ": ";
    ss << std::setw(13) << mNodalS(0) << std::setw(13) << mNodalZ(0) << std::endl;
    ss << "    " << std::setw(width) << "..." << std::endl;
    ss << "    " << std::setw(width) << mNumNodesInner - 1 << ": ";
    ss << std::setw(13) << mNodalS(mNumNodesInner - 1) << std::setw(13) << mNodalZ(mNumNodesInner - 1) << std::endl;
    if (hasABC()) {
        ss << "    " << std::setw(width) << "..." << std::endl;
        ss << "    " << std::setw(width) << "absorbing boundary ends at" << std::endl;
        ss << "    " << std::setw(width) << getNumNodes() << ": ";
        ss << std::setw(13) << mNodalS(getNumNodes() - 1) << std::setw(13) << mNodalZ(getNumNodes() - 1) << std::endl;
    }
    ss << "  Elemental Variables_______________________________________" << std::endl;
    widthname = -1;
    for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mElementalVariables.begin(); it != mElementalVariables.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ": ";
        ss << std::setw(13) << it->second(0) << ", ..., ";
        ss << std::setw(13) << it->second(getNumQuads() - 1) << std::endl;        
    }
    ss << "  Side Sets_________________________________________________" << std::endl;
    widthname = -1;
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ":   ";
        int pair = 0;
        for (int q = 0; q < mNumQuadsInner; q++) {
            pair += (int)(it->second(q) >= 0);
        }
        ss << pair;
        if (hasABC()) {
            pair = 0;
            for (int q = mNumQuadsInner; q < getNumQuads(); q++) {
                pair += (int)(it->second(q) >= 0);
            }
            ss << " normal + " << pair << " absorbing";
        }
        ss << " edges" << std::endl;
    }
    ss << "  Miscellaneous_____________________________________________" << std::endl;
    ss << "    Distance Tolerance / m   =   " << mDistTolerance << std::endl;
    if (!isCartesian()) {
        ss << "  External__________________________________________________" << std::endl;
        ss << "    Num. Ellipticity Spline Knots   =   " << mEllipKnots.size() << std::endl;
    }
    ss << "======================= Exodus Model =======================\n" << std::endl;
    return ss.str();
}

void ExodusModel::buildInparam(ExodusModel *&exModel, const Parameters &par, 
    AttParameters *&attPar, int verbose) {
    if (exModel) {
        delete exModel;
    }

    std::string exfile = par.getValue<std::string>("MODEL_1D_EXODUS_MESH_FILE");
    exfile = Parameters::sInputDirectory + "/" + exfile;
    exModel = new ExodusModel(exfile);

    exModel->mN_ABC = par.getValue<int>("ABC_ELEMENTS");

    exModel->initialize();
    if (verbose) {
        XMPI::cout << exModel->verbose();
    }
    
    // form attenuation parameters
    if (exModel->hasAttenuation()) {
        int nr_lin_solids = (int)exModel->mGlobalVariables.at("nr_lin_solids");
        double f_min = exModel->mGlobalVariables.at("f_min");
        double f_max = exModel->mGlobalVariables.at("f_max");
        double f_ref = exModel->mGlobalVariables.at("f_ref");
        RDColX w(nr_lin_solids), y(nr_lin_solids);
        for (int i = 0; i < nr_lin_solids; i++) {
            std::stringstream sw, sy;
            sw << "w_" << i;
            sy << "y_" << i;
            w(i) = exModel->mGlobalVariables.at(sw.str());
            y(i) = exModel->mGlobalVariables.at(sy.str());
        }
        if (attPar) {
            delete attPar;
        }
        attPar = new AttParameters(nr_lin_solids, f_min, f_max, f_ref, w, y);
    } else {
        if (attPar) {
            delete attPar;
            attPar = 0;
        }
    }
    
    // ellipticity
    if (exModel->isCartesian()) {
        return;
    }
    std::string emode = par.getValue<std::string>("MODEL_3D_ELLIPTICITY_MODE");
    if (boost::iequals(emode, "off")) {
        // no ellipticity
        Geodesy::setup(exModel->getROuter(), 0., RDColX::Zero(0), RDColX::Zero(0));
    } else {
        double inv_f = par.getValue<double>("MODEL_3D_ELLIPTICITY_INVF");
        if (inv_f <= 0.) {
            throw std::runtime_error("ExodusModel::buildInparam || Invalid flattening.");
        } 
        Geodesy::setup(exModel->getROuter(), 1. / inv_f, exModel->mEllipKnots, exModel->mEllipCoeffs);
    }
}

