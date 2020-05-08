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

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include "XMath.h"

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
    
    // distance tolerance
    mDistTolerance = DBL_MAX;
    mHmax = 0;
    mHmin = 6371000;
    for (int i = 0; i < getNumQuads(); i++) {
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
        mDistTolerance = std::min({dist0, dist1, dist2, dist3, mDistTolerance});
        mHmax = std::max({dist0, dist1, dist2, dist3, mHmax});
        mHmin = std::min({dist0, dist1, dist2, dist3, mHmin});
    }
    mHmax = 1000 * mHmax;
    mHmin = 1000 * mHmin;
    
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
    
    // do not read ellipticity for cartesian
    bool cartesian = false;
    for (int i = 0; i < mGlobalRecordsRaw.size(); i++) {
        std::vector<std::string> substrs = Parameters::splitString(mGlobalRecordsRaw[i], "=");
        if (boost::iequals(boost::trim_copy(substrs[0]), "crdsys") &&
            boost::iequals(boost::trim_copy(substrs[1]), "cartesian")) {
            cartesian = true;
            break;
        }
    }
    if (cartesian) {
        mEllipKnots = RDColX::Zero(0);
        mEllipCoeffs = RDColX::Zero(0);
        // set radius to PREM
        mGlobalVariables.insert(std::pair<std::string, double>("radius", 6371e3));
    } else {
        // ellipticity
        reader.read2D("ellipticity", dbuffer);
        mEllipKnots = dbuffer.row(0).transpose();
        mEllipCoeffs = dbuffer.row(1).transpose();
    }
    
    // element var names
    reader.readString("name_elem_var", mElementalVariableNames_all);    
    
    // axis is needed
    if (cartesian) {
        mSSNameAxis = "x0";
    } else {
        mSSNameAxis = "t0";
        // delete this "if" after Martin fixes the issue
        if (std::find(mSideSetNames.begin(), mSideSetNames.end(), "t1") != mSideSetNames.end()) {
            mSSNameAxis = "t1";
        }
    }
    int axisIdx = std::find(mSideSetNames.begin(), mSideSetNames.end(), mSSNameAxis) - mSideSetNames.begin();
    IColX axial = mSideSetValues.col(axisIdx);
    
    // depth-dependent variable values
    std::vector<double> coords;
    std::vector<std::pair<int, int>> quad_nodes;
    for (int iquad = 0; iquad < getNumQuads(); iquad++) {
        int axisSide = axial(iquad);
        if (axisSide > -1) {
            // coords
            int otherNode = (axisSide == 3) ? 0 : axisSide + 1;
            double z1 = mNodalZ(mConnectivity(iquad, axisSide));
            double z2 = mNodalZ(mConnectivity(iquad, otherNode));
            // negative values not needed
            if (std::min(z1, z2) < -mDistTolerance) {
                continue;
            }
            // move both ends inward by tolerance
            z1 += (z2 - z1) / std::abs(z2 - z1) * mDistTolerance;
            z2 -= (z2 - z1) / std::abs(z2 - z1) * mDistTolerance;
            auto ir1 = coords.insert(std::upper_bound(coords.begin(), coords.end(), z1), z1);
            quad_nodes.insert(quad_nodes.begin() + (ir1 - coords.begin()),
                              std::pair<int, int>(iquad, axisSide));
            auto ir2 = coords.insert(std::upper_bound(coords.begin(), coords.end(), z2), z2);
            quad_nodes.insert(quad_nodes.begin() + (ir2 - coords.begin()),
                              std::pair<int, int>(iquad, otherNode));
        }
    }
    mElementalVariableCoords_axis = RDColX(coords.size());
    for (int j = 0; j < coords.size(); j++) {
        mElementalVariableCoords_axis(j) = coords[j];
    }
    
    // names of axial variables (remove index (0-3 for nodes) at the end of var name)
    for (int iname = 0; iname < mElementalVariableNames_all.size(); iname++) {
        std::string varName = mElementalVariableNames_all[iname];
        if (varName.substr(varName.length() - 2, 1) == std::string("_")) {
            std::string vname = varName.substr(0, varName.length() - 2);
            std::string inode_str = varName.substr(varName.length() - 1, 1);
            if (inode_str == "0") {
                mElementalVariableNames_axis.push_back(vname);
            }
        } else {
            mElementalVariableNames_axis.push_back(varName);
        }
    }
    
    // elemental variables only on axis
    mElementalVariableValues_axis = RDMatXX::Zero(mElementalVariableCoords_axis.size(), 
        mElementalVariableNames_axis.size());
    mABC_Vmax = -1;
    for (int iname = 0; iname < mElementalVariableNames_axis.size(); iname++) {
        std::string varName = mElementalVariableNames_axis[iname];
        if (std::find(mElementalVariableNames_all.begin(), mElementalVariableNames_all.end(), varName + "_0") 
        != mElementalVariableNames_all.end()) {
            // nodal dependent
            std::array<RDColX, 4> buffer;
            for (int inode = 0; inode < 4; inode++) {
                std::stringstream ss;
                ss << varName << "_" << inode;
                std::string varNameToRead = ss.str();
                int index = std::find(mElementalVariableNames_all.begin(), 
                    mElementalVariableNames_all.end(), varNameToRead) - 
                    mElementalVariableNames_all.begin() + 1;
                ss.str("");
                ss << "vals_elem_var" << index << "eb1";
                reader.read1D(ss.str(), buffer[inode]);
            }
            for (int idep = 0; idep < coords.size(); idep++) {
                mElementalVariableValues_axis(idep, iname) = 
                    buffer[quad_nodes[idep].second](quad_nodes[idep].first);
            }
        } else {
            // nodal independent
            RDColX buffer;
            std::stringstream ss;
            int index = std::find(mElementalVariableNames_all.begin(), 
                mElementalVariableNames_all.end(), varName) - 
                mElementalVariableNames_all.begin() + 1;
            ss << "vals_elem_var" << index << "eb1";
            reader.read1D(ss.str(), buffer);
            for (int idep = 0; idep < coords.size(); idep++) {
                mElementalVariableValues_axis(idep, iname) = 
                    buffer(quad_nodes[idep].first);
            }
        }
        if (varName.substr(0, 2) == "VP") {
            mABC_Vmax = std::max(mABC_Vmax, mElementalVariableValues_axis.col(iname).maxCoeff());
        }
    }            
    
    // element-wise
    mElementalVariableNames_elem.push_back("element_type");
    mElementalVariableNames_elem.push_back("dt");
    mElementalVariableValues_elem = RDMatXX::Zero(getNumQuads(), 
        mElementalVariableNames_elem.size());
    for (int i = 0; i < mElementalVariableNames_elem.size(); i++) {
        RDColX buffer;
        std::stringstream ss;
        int index = std::find(mElementalVariableNames_all.begin(), 
            mElementalVariableNames_all.end(), mElementalVariableNames_elem[i]) - 
            mElementalVariableNames_all.begin() + 1;
        ss << "vals_elem_var" << index << "eb1";
        reader.read1D(ss.str(), buffer);
        mElementalVariableValues_elem.col(i) = buffer;
    }  
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
    XMPI::bcast(mDistTolerance);
    XMPI::bcast(mHmax);
    XMPI::bcast(mABC_Vmax);
    
    XMPI::bcast(mElementalVariableNames_all);
    XMPI::bcast(mElementalVariableNames_elem);
    XMPI::bcast(mElementalVariableNames_axis);
    XMPI::bcastEigen(mElementalVariableValues_elem);
    XMPI::bcastEigen(mElementalVariableValues_axis);
    XMPI::bcastEigen(mElementalVariableCoords_axis);
    
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
    for (int i = 0; i < mElementalVariableNames_axis.size(); i++) {
        mElementalVariables_axis.insert(std::pair<std::string, RDColX>(mElementalVariableNames_axis[i], 
            mElementalVariableValues_axis.col(i)));
    }
    for (int i = 0; i < mElementalVariableNames_elem.size(); i++) {
        mElementalVariables_elem.insert(std::pair<std::string, RDColX>(mElementalVariableNames_elem[i], 
            mElementalVariableValues_elem.col(i)));
    }

    // side sets
    for (int i = 0; i < mSideSetNames.size(); i++) {
        mSideSets.insert(std::pair<std::string, IColX>(mSideSetNames[i],
            mSideSetValues.col(i)));
    }

    // name of axis and surface sets
    std::string sphere_axis_name = "t0";
    if (!isCartesian()) {
        if (mSideSets.find("t1") != mSideSets.end()) {
            sphere_axis_name = "t1";
        }
    }
    
    mSSNameAxis = isCartesian() ? "x0" : sphere_axis_name;
    mSSNameSurface = isCartesian() ? "y1" : "r1";
    
    if (!isCartesian()) {
        if (mSideSets.find("t1") != mSideSets.end()) {
            mSSNameRightB = "t0";
        } else {
            mSSNameRightB = "N/A";
        }
        if (mSideSets.find("r0") != mSideSets.end()) {
            mSSNameBottomB = "r0";
        } else {
            mSSNameBottomB = "N/A";
        }
    } else {
        mSSNameRightB = "x1";
        mSSNameBottomB = "y0";
    }

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

    // rotate nodes of axial elements such that side 3 is on axis, side 1 is the right boundary and side 0 is the lower boundary
    MultilevelTimer::begin("Process Exodus Axis", 2);
    for (int axialQuad = 0; axialQuad < getNumQuads(); axialQuad++) {

        int shift=0;
        int axialSide = getSideAxis(axialQuad);
        if (axialSide != 3 & axialSide != -1) {shift = axialSide - 3;}
        if (shift==0) {
            continue;
        }
        
        // side sets
        for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
            if (it->second(axialQuad) != -1) {
                it->second(axialQuad) = Mapping::period0123(it->second(axialQuad) - shift);
            }
        }
    }
    MultilevelTimer::end("Process Exodus Axis", 2);

    MultilevelTimer::begin("Process Absorbing Boundaries", 2);
    // extending mesh to incorporate absorbing boundaries
    mABfield = IColX::Constant(getNumQuads(), -1);
    mNumQuadsInner = getNumQuads();
    mNumNodesInner = getNumNodes();

    if (!isCartesian() && (mHasSpongeBoundary || mHasStaceyABC)) {
        computeNodalRTheta();
    }

    if (mHasSpongeBoundary) {
        addAbsorbingBoundaryElements();
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
}

void ExodusModel::addAbsorbingBoundaryElements() {
    // calculate boundary width in case of wavelengths input
    if (mN_maxWL_ABC > 0) {
        mABCwidth = mN_maxWL_ABC * mTSource * mABC_Vmax / 1000;
    }
    
    if (mHasMeshExtension) {
        if (!isCartesian()) {
            throw std::runtime_error("ExodusModel::AddAbsorbingBoundaryElements || "
                "Mesh extension not implemented for non-Cartesian meshes.");
        }
        if (mElementalVariables_axis.at("edge_aspect_ratio").maxCoeff() > 1.5) {
            throw std::runtime_error("ExodusModel::AddAbsorbingBoundaryElements || "
                "Mesh extension not implemented for meshes with refinements.");
        }
        
        //number of extended elements
        if (mN_ABC < 0) {
            mN_ABC = ceil(round(10000 * mABCwidth / mHmax) / 10);
        }
        setMeshEdge();
        mInnerBoundaryCorner = mMeshEdgeCorner;
        XMPI::cout << "Building right boundary extension..." << XMPI::endl;
        extendRightSide(getNumQuads(), getNumNodes(), mN_ABC, mHmax);
        XMPI::cout << "Building bottom boundary extension..." << XMPI::endl << XMPI::endl;
        extendBottomSide(getNumQuads(), getNumNodes(), mN_ABC, mHmax);
        setMeshEdge();
    } else {
        setMeshEdge();
        setInnerBoundary();
        makeABField(mInnerBoundaryCorner);
    }
}

void ExodusModel::computeNodalRTheta() {
    mNodalR = mNodalS.schur(mNodalS) + mNodalZ.schur(mNodalZ);
    mNodalR = mNodalR.cwiseSqrt();
    
    mNodalTheta = RDColX::Zero(getNumNodes());
    for (int i = 0; i < getNumNodes(); i++) {
        mNodalTheta(i) = mNodalR(i) > tinyDouble ? mNodalZ(i) / mNodalR(i) : 0;
        mNodalTheta(i) = acos(mNodalTheta(i));
    }
}

void ExodusModel::setMeshEdge() {
    if (isCartesian()) {
        mMeshEdgeCorner(0) = mNodalS.maxCoeff();
        mMeshEdgeCorner(1) = mNodalZ.minCoeff();    
    } else {
        mMeshEdgeCorner(0) = mNodalTheta.minCoeff();;
        mMeshEdgeCorner(1) = mNodalR.maxCoeff();
    }
}

void ExodusModel::setInnerBoundary() {
    RDCol2 innerBoundary;
    RDCol2 maxdiff(0, 0);
    
    RDCol2 crds0, crds1, crds2, crds3;
    for (int i = 0; i < getNumQuads(); i++) {
        if (isCartesian()) {
            crds0 << mNodalS(mConnectivity(i, 0)), mNodalZ(mConnectivity(i, 0));
            crds1 << mNodalS(mConnectivity(i, 1)), mNodalZ(mConnectivity(i, 1));
            crds2 << mNodalS(mConnectivity(i, 2)), mNodalZ(mConnectivity(i, 2));
            crds3 << mNodalS(mConnectivity(i, 3)), mNodalZ(mConnectivity(i, 3));
        } else {
            crds0 << mNodalTheta(mConnectivity(i, 0)), mNodalR(mConnectivity(i, 0));
            crds1 << mNodalTheta(mConnectivity(i, 1)), mNodalR(mConnectivity(i, 1));
            crds2 << mNodalTheta(mConnectivity(i, 2)), mNodalR(mConnectivity(i, 2));
            crds3 << mNodalTheta(mConnectivity(i, 3)), mNodalR(mConnectivity(i, 3));
        }
        RDCol2 diff0 = crds0 - crds1;
        RDCol2 diff1 = crds1 - crds2;
        RDCol2 diff2 = crds2 - crds3;
        RDCol2 diff3 = crds3 - crds0;
        
        diff0 = diff0.cwiseAbs();
        diff1 = diff1.cwiseAbs();
        diff2 = diff2.cwiseAbs();
        diff3 = diff3.cwiseAbs();
        
        maxdiff(0) = std::max({diff0(0), diff1(0), diff2(0), diff3(0), maxdiff(0)}); // s or theta
        maxdiff(1) = std::max({diff0(1), diff1(1), diff2(1), diff3(1), maxdiff(1)}); // z or r
    }
    
    RDCol2 dir_fac(-1, 1);
    if (mABCwidth > 0) { // wavelength or distance input
        RDCol2 N;
        N(0) = XMath::RobustRoundUp(mABCwidth / maxdiff(0), tinyDouble);
        N(1) = XMath::RobustRoundUp(mABCwidth / maxdiff(1), tinyDouble);
        if (!isCartesian()) {
            N(1) = XMath::RobustRoundUp((mABCwidth / mMeshEdgeCorner(1)) / maxdiff(1), tinyDouble); // theta
        }
        N(0) = std::max({N(0), 1.});
        N(1) = std::max({N(1), 1.});
        
        mN_ABC = N.minCoeff();
        
        innerBoundary(0) = mMeshEdgeCorner(0) + dir_fac(0) * N(0) * maxdiff(0); // s or theta
        innerBoundary(1) = mMeshEdgeCorner(1) + dir_fac(1) * N(1) * maxdiff(1); // z or r
    } else { // elements input
        innerBoundary(0) = mMeshEdgeCorner(0) + dir_fac(0) * mN_ABC * maxdiff(0);
        innerBoundary(1) = mMeshEdgeCorner(1) + dir_fac(1) * mN_ABC * maxdiff(1);
    }
    
    mInnerBoundaryCorner = innerBoundary;
}

void ExodusModel::makeABField(RDCol2 inner) {
    if (!isCartesian()) {
        inner = inner.colwise().reverse();
    }
    
    RDCol2 crds0, crds1, crds2, crds3;
    for (int i = 0; i < getNumQuads(); i++) {
        if (isCartesian()) {
            crds0 << mNodalS(mConnectivity(i, 0)), mNodalZ(mConnectivity(i, 0));
            crds1 << mNodalS(mConnectivity(i, 1)), mNodalZ(mConnectivity(i, 1));
            crds2 << mNodalS(mConnectivity(i, 2)), mNodalZ(mConnectivity(i, 2));
            crds3 << mNodalS(mConnectivity(i, 3)), mNodalZ(mConnectivity(i, 3));
        } else {
            crds0 << mNodalTheta(mConnectivity(i, 0)), mNodalR(mConnectivity(i, 0));
            crds1 << mNodalTheta(mConnectivity(i, 1)), mNodalR(mConnectivity(i, 1));
            crds2 << mNodalTheta(mConnectivity(i, 2)), mNodalR(mConnectivity(i, 2));
            crds3 << mNodalTheta(mConnectivity(i, 3)), mNodalR(mConnectivity(i, 3));
        }
        RDCol2 outerCrds(std::max({crds0(0), crds1(0), crds2(0), crds3(0)}), std::min({crds0(1), crds1(1), crds2(1), crds3(1)}));
        
        if ((outerCrds(0) > inner(0)) && (outerCrds(1) < inner(1))) {
            mABfield(i) = 3; // corner
        } else if (outerCrds(0) > inner(0)) { // right bry
            mABfield(i) = 1;
        } else if (outerCrds(1) < inner(1)) {
            mABfield(i) = 2; // bottom bry
        } else {
            continue;
        }
        mNumQuadsInner -= 1;
    }
}

void ExodusModel::extendRightSide(int nQuads, int nNodes, int N, double h) {
    // this function adds new elements onto:
    // - nNodalS, nNodalZ
    // - mConnectivity
    // - mSideSets
    // - mElementalVariables_elem
    // - mABfield
    
    std::vector<int> rightBQuadTags;
    for (int tag = 0; tag < nQuads; tag++) {
        if (getSideRightB(tag) >= 0) {
            rightBQuadTags.push_back(tag);
        }
    }
    int ExtNr = rightBQuadTags.size() * N;
    mNodalS.conservativeResize(nNodes + ExtNr + N, 1);
    mNodalZ.conservativeResize(nNodes + ExtNr + N, 1);
    mConnectivity.conservativeResize(nQuads + ExtNr, 4);
    mABfield.conservativeResize(nQuads + ExtNr, 1);
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        it->second.conservativeResize(nQuads + ExtNr, 1);
        it->second.segment(nQuads, ExtNr) = IColX::Constant(ExtNr, 1, -1);
    }
    for (auto it = mElementalVariables_elem.begin(); it != mElementalVariables_elem.end(); it++) {
        it->second.conservativeResize(nQuads + ExtNr, 1);
    }

    int BottomRightNode = mConnectivity(rightBQuadTags[0], getSideRightB(rightBQuadTags[0]));
    mNodalS.segment(nNodes, N) = RDColX::LinSpaced(N, mNodalS(BottomRightNode) + h, mNodalS(BottomRightNode) + N * h);
    mNodalZ.segment(nNodes, N) = RDColX::Constant(N, 1, mNodalZ(BottomRightNode));

    for (int i = 0; i < rightBQuadTags.size(); i++) {
        int myQuad = rightBQuadTags[i];
        int side = getSideRightB(myQuad);
        int node3_ini = mConnectivity(myQuad, side + 1);
        int node1_ini = nNodes + i * N;
        mConnectivity.row(nQuads + N * i) << mConnectivity(myQuad, side), node1_ini, node1_ini + N, node3_ini;
        for (int j = 0; j < N - 1; j++) {
            int node0 = node1_ini + j;
            mConnectivity.row(nQuads + i * N + j + 1) << node0, node0 + 1, node0 + N + 1, node0 + N;
        }
        mNodalS.segment(node1_ini + N, N) = RDColX::LinSpaced(N, mNodalS(node3_ini) + h, mNodalS(node3_ini) + N * h);
        mNodalZ.segment(node1_ini + N, N) = RDColX::Constant(N, 1, mNodalZ(node3_ini));
        mSideSets.at(mSSNameRightB)(myQuad) = -1;
        for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
            it->second.segment(nQuads + i * N, N) = IColX::Constant(N, 1, it->second(myQuad));
        }
        mSideSets.at(mSSNameRightB)(nQuads + (i + 1) * N - 1) = 1;
        mABfield(myQuad) = 0;
        mABfield.block(nQuads + i * N, 0, N, 1) = IColX::Constant(N, 1, 1);
        for (auto it = mElementalVariables_elem.begin(); it != mElementalVariables_elem.end(); it++) {
            it->second.segment(nQuads + i * N, N) = RDColX::Constant(N, 1, it->second(myQuad));
        }
    }
}

void ExodusModel::extendBottomSide(int nQuads, int nNodes, int N, double h) {
    // this function adds new elements onto:
    // - nNodalS, nNodalZ
    // - mConnectivity
    // - mSideSets
    // - mElementalVariables_elem
    // - mABfield
    
    // handling of corner elements is implemented in this function 
    
    std::vector<int> lowerBQuadTags;
    for (int tag = 0; tag < nQuads; tag++) {
        if (getSideLowerB(tag) >= 0) {
            lowerBQuadTags.push_back(tag);
        }
    }
    int ExtNr = lowerBQuadTags.size() * N;
    mNodalS.conservativeResize(nNodes + ExtNr + N, 1);
    mNodalZ.conservativeResize(nNodes + ExtNr + N, 1);
    mConnectivity.conservativeResize(nQuads + ExtNr, 4);
    mABfield.conservativeResize(nQuads + ExtNr, 1);
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        it->second.conservativeResize(nQuads + ExtNr, 1);
        it->second.segment(nQuads, ExtNr) = IColX::Constant(ExtNr, 1, -1);
    }
    for (auto it = mElementalVariables_elem.begin(); it != mElementalVariables_elem.end(); it++) {
        it->second.conservativeResize(nQuads + ExtNr, 1);
    }

    int BottomLeftNode = mConnectivity(lowerBQuadTags[0], 0);
    mNodalS.segment(nNodes, N) = RDColX::Constant(N, 1, mNodalS(BottomLeftNode));
    mNodalZ.segment(nNodes, N) = RDColX::Constant(N, 1, mNodalZ(BottomLeftNode)) - RDColX::LinSpaced(N, h, N * h);

    int CornerQuadTag;
    for (int i = 0; i < lowerBQuadTags.size(); i++) {
        int myQuad = lowerBQuadTags[i];
        int side = getSideLowerB(myQuad);
        int node2_ini = mConnectivity(myQuad, side + 1);
        int node0_ini = nNodes + i * N;
        mConnectivity.row(nQuads + N * i) << node0_ini, node0_ini + N, node2_ini, mConnectivity(myQuad, side);
        for (int j = 0; j < N - 1; j++) {
            int node0 = node0_ini + j + 1;
            mConnectivity.row(nQuads + i * N + j + 1) << node0, node0 + N, node0 + N - 1, node0 - 1;
        }
        mNodalS.segment(node0_ini + N, N) = RDColX::Constant(N, 1, mNodalS(node2_ini));
        mNodalZ.segment(node0_ini + N, N) = RDColX::Constant(N, 1, mNodalZ(node2_ini)) - RDColX::LinSpaced(N, h, N * h);
        mSideSets.at(mSSNameBottomB)(myQuad) = -1;
        for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
            it->second.segment(nQuads + i * N, N) = IColX::Constant(N, 1, it->second(myQuad));
        }
        mSideSets.at(mSSNameBottomB)(nQuads + (i + 1) * N - 1) = 0;

        if (mABfield(myQuad, 1) == -1) {
            mABfield(myQuad, 1) = 0;
            mABfield.block(nQuads + i * N, 0, N, 1) = IColX::Constant(N, 1, 2);
        } else if (mABfield(myQuad, 1) == 0) {
            mABfield.block(nQuads + i * N, 0, N, 1) = IColX::Constant(N, 1, 2);
        } else if (mABfield(myQuad, 1) == 1) {
            mABfield.block(nQuads + i * N, 0, N, 1) = IColX::Constant(N, 1, 3);
        }
        for (auto it = mElementalVariables_elem.begin(); it != mElementalVariables_elem.end(); it++) {
            it->second.segment(nQuads + i * N, N) = RDColX::Constant(N, 1, it->second(myQuad));
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
    ss << "    " << std::setw(width) << mNumNodesInner << ": ";
    ss << std::setw(13) << mNodalS(mNumNodesInner - 1) << std::setw(13) << mNodalZ(mNumNodesInner - 1) << std::endl;
    if (mHasSpongeBoundary && mHasMeshExtension) {
        ss << "    " << std::setw(width) << "..." << std::endl;
        ss << "    " << std::setw(width) << "extended boundary ends at" << std::endl;
        ss << "    " << std::setw(width) << getNumNodes() << ": ";
        ss << std::setw(13) << mNodalS(getNumNodes() - 1) << std::setw(13) << mNodalZ(getNumNodes() - 1) << std::endl;
    } else if (mHasSpongeBoundary) {
        ss << "    " << std::setw(width) << "..." << std::endl;
        ss << "    " << std::setw(width) << "intruding boundary ends at" << std::endl;
        double s,z;
        if (isCartesian()) {
            s = mInnerBoundaryCorner(0);
            z = mInnerBoundaryCorner(1);
        } else {
            s = mInnerBoundaryCorner(1) * sin(mInnerBoundaryCorner(0));
            z = mInnerBoundaryCorner(1) * cos(mInnerBoundaryCorner(0));
        }
        ss << std::setw(19) << s << std::setw(13) << z << std::endl;
    }
    ss << "  Elemental Variables_______________________________________" << std::endl;
    widthname = -1;
    for (auto it = mElementalVariables_elem.begin(); it != mElementalVariables_elem.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mElementalVariables_elem.begin(); it != mElementalVariables_elem.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ": ";
        ss << std::setw(13) << it->second(0) << ", ..., ";
        ss << std::setw(13) << it->second(getNumQuads() - 1) << std::endl;
    }
    
    ss << "  Depth-dependent Elemental Variables_______________________" << std::endl;
    widthname = -1;
    for (auto it = mElementalVariables_axis.begin(); it != mElementalVariables_axis.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mElementalVariables_axis.begin(); it != mElementalVariables_axis.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ": ";
        ss << std::setw(13) << it->second(0) << ", ..., ";
        ss << std::setw(13) << it->second(mElementalVariableCoords_axis.size() - 1) << std::endl;        
    }
    
    ss << "  Side Sets_________________________________________________" << std::endl;
    widthname = -1;
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        widthname = std::max(widthname, (int)(it->first.length()));
    }
    for (auto it = mSideSets.begin(); it != mSideSets.end(); it++) {
        ss << "    " << std::setw(widthname) << it->first << ":   ";
        int pair = 0;
        for (int q = 0; q < getNumQuads(); q++) {
            pair += (int)(it->second(q) >= 0);
        }
        ss << pair << " edges" << std::endl;
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

    exModel->mHasSpongeBoundary = par.getValue<bool>("ABC_SPONGE_BOUNDARIES");
    exModel->mHasStaceyABC = par.getValue<bool>("ABC_STACEY_BOUNDARIES");

    if (exModel->mHasSpongeBoundary) {
        exModel->mHasMeshExtension = par.getValue<bool>("ABC_SPONGE_BOUNDARIES_EXTEND_MESH");
        exModel->mHasModelExtension = par.getValue<bool>("ABC_SPONGE_BOUNDARIES_EXTEND_MODEL");
        exModel->mTSource = 2 * par.getValue<double>("SOURCE_STF_HALF_DURATION");
        
        std::string mstr = par.getValue<std::string>("ABC_SPONGE_BOUNDARIES_WIDTH");
        std::vector<std::string> strs = Parameters::splitString(mstr, "$");
        std::string format(strs[0]);
        if (boost::iequals(format, "wavelengths")) {
            exModel->mN_maxWL_ABC = boost::lexical_cast<int>(strs[1]);
        } else if (boost::iequals(format, "distance")) {
            exModel->mABCwidth = boost::lexical_cast<double>(strs[1]);
        } else if (boost::iequals(format, "elements")) {
            exModel->mN_ABC = boost::lexical_cast<int>(strs[1]);
        } else {
            throw std::runtime_error("ExodusModel::buildInparam || "
                "Unknown ABC extension format " + format + ".");
        }
    }

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
        Geodesy::setup(exModel->getROuter(), 0., RDColX::Zero(0), RDColX::Zero(0));
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

bool ExodusModel::isIsotropic() const {
    return std::find(mElementalVariableNames_all.begin(), 
        mElementalVariableNames_all.end(), "VP_0") != mElementalVariableNames_all.end();
}

double ExodusModel::getElementalVariables(const std::string &varName, int quadTag) const {
    if (varName == "element_type" || varName == "dt") {
        return mElementalVariables_elem.at(varName)(quadTag);
    }
    
    if (varName.substr(varName.length() - 2, 1) == std::string("_")) {
        // nodal dependent
        std::string vname = varName.substr(0, varName.length() - 2);
        int inode = boost::lexical_cast<int>(varName.substr(varName.length() - 1, 1));
        double coord = 0.;
        double coord_cen = 0.;
        if (isCartesian()) {
            coord = mNodalZ(mConnectivity(quadTag, inode));
            coord_cen += mNodalZ(mConnectivity(quadTag, 0));
            coord_cen += mNodalZ(mConnectivity(quadTag, 1));
            coord_cen += mNodalZ(mConnectivity(quadTag, 2));
            coord_cen += mNodalZ(mConnectivity(quadTag, 3));
        } else {
            coord = sqrt(pow(mNodalS(mConnectivity(quadTag, inode)), 2) + pow(mNodalZ(mConnectivity(quadTag, inode)), 2));
            coord_cen += sqrt(pow(mNodalS(mConnectivity(quadTag, 0)), 2) + pow(mNodalZ(mConnectivity(quadTag, 0)), 2));
            coord_cen += sqrt(pow(mNodalS(mConnectivity(quadTag, 1)), 2) + pow(mNodalZ(mConnectivity(quadTag, 1)), 2));
            coord_cen += sqrt(pow(mNodalS(mConnectivity(quadTag, 2)), 2) + pow(mNodalZ(mConnectivity(quadTag, 2)), 2));
            coord_cen += sqrt(pow(mNodalS(mConnectivity(quadTag, 3)), 2) + pow(mNodalZ(mConnectivity(quadTag, 3)), 2));
        }
        coord_cen /= 4.;
        // move coord slightly to center
        if (coord < coord_cen) {
            coord += mDistTolerance;
        } else {
            coord -= mDistTolerance;
        }
        RDColX diff = (mElementalVariableCoords_axis.array() - coord).array().abs();
        int minIndex = 0;
        double min = diff.minCoeff(&minIndex);
        return mElementalVariables_axis.at(vname)(minIndex);
    } else {
        // nodal independent
        double coord_cen = 0.;
        if (isCartesian()) {
            coord_cen += mNodalZ(mConnectivity(quadTag, 0));
            coord_cen += mNodalZ(mConnectivity(quadTag, 1));
            coord_cen += mNodalZ(mConnectivity(quadTag, 2));
            coord_cen += mNodalZ(mConnectivity(quadTag, 3));
        } else {
            coord_cen += sqrt(pow(mNodalS(mConnectivity(quadTag, 0)), 2) + pow(mNodalZ(mConnectivity(quadTag, 0)), 2));
            coord_cen += sqrt(pow(mNodalS(mConnectivity(quadTag, 1)), 2) + pow(mNodalZ(mConnectivity(quadTag, 1)), 2));
            coord_cen += sqrt(pow(mNodalS(mConnectivity(quadTag, 2)), 2) + pow(mNodalZ(mConnectivity(quadTag, 2)), 2));
            coord_cen += sqrt(pow(mNodalS(mConnectivity(quadTag, 3)), 2) + pow(mNodalZ(mConnectivity(quadTag, 3)), 2));
        }
        coord_cen /= 4.;
        RDColX diff = (mElementalVariableCoords_axis.array() - coord_cen).array().abs();
        int minIndex = 0;
        double min = diff.minCoeff(&minIndex);
        return mElementalVariables_axis.at(varName)(minIndex);
    }
}

bool ExodusModel::findDiscontinuity(double &layer, const std::string &varName, double val, double upper, 
    double lower, const std::string &compType, const bool from_bottom) const {

    std::vector<double> z(mElementalVariableCoords_axis.data(), 
                          mElementalVariableCoords_axis.data() + mElementalVariableCoords_axis.rows());
    
    auto it0 = std::upper_bound(z.begin(), z.end(), getROuter() - lower);
    auto it1 = std::upper_bound(z.begin(), z.end(), getROuter() - upper);
    
    if (it0 == z.end() || it1 == z.begin()) {
        throw std::runtime_error("ExodusModel::findDiscontinuity || Search range " 
        + boost::lexical_cast<std::string>(getROuter() - lower) + " - " + boost::lexical_cast<std::string>(getROuter() - upper)
        + " out of bounds of Exodus Model.");
    }
    
    int i0 = it0 - z.begin();
    int i1 = it1 - z.begin();
    
    RDColX vars = mElementalVariables_axis.at(varName).segment(i0, i1-i0);
    if (compType.compare("less") == 0) {
        vars *= -1;
        val *= -1;
    } else if (compType.compare("greater") != 0) {
        throw std::runtime_error("ExodusModel::findDiscontinuity || Invalid comparison type. Permitted values are 'greater' and 'less'.");
    }
    
    int idisc;
    bool found;
    if (from_bottom) {
        for (int i = 0; i < vars.size(); i++) {
            if (vars(i) > val) {
                idisc = i;
                found = true;
                break;
            }
        }
    } else {
        for (int i = vars.size() - 1; i >= 0; i--) {
            if (vars(i) > val) {
                idisc = i;
                found = true;
                break;
            }
        }
    }
    
    double z0 = z[i0 + idisc];
    double z1 = z[i0 + idisc + 1];
    
    if (z1 - z0 < 2 * mDistTolerance + tinyDouble) {
        layer = z0 + mDistTolerance;
    } else {
        layer = z0 - mDistTolerance;
    }
    
    layer = getROuter() - layer;
    return found;    
}

void ExodusModel::findClosestMeshLine(double &layer) const {
    double zlayer = getROuter() - layer;
    std::vector<double> z(mElementalVariableCoords_axis.data(), 
                          mElementalVariableCoords_axis.data() + mElementalVariableCoords_axis.rows());
                          
    auto it = std::upper_bound(z.begin(), z.end(), zlayer);
    
    if (abs(*(it) - zlayer) < abs(*(it - 1) - zlayer)) {
        layer = getROuter() - *(it);
    } else {
        layer = getROuter() - *(it - 1);
    }
};




