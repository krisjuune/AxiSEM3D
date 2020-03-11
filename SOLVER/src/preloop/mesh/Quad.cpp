// Quad.cpp
// created by Kuangdai on 5-May-2016 
// general Quad element 
// created from Exodus model, elements of AxiSEM3D mesh

#include "Quad.h"
#include "ExodusModel.h"

#include "ABCParameters.h"

#include "SphericalMapping.h"
#include "LinearMapping.h"
#include "SemiSphericalMapping.h"
#include "SpectralConstants.h"

#include "Material.h"
#include "Relabelling.h"

#include "Elastic.h"
#include "Acoustic.h"

#include "Gradient.h"
#include "PreloopGradient.h"
#include "SolidElement.h"
#include "FluidElement.h"
#include "Domain.h"

#include "GLLPoint.h" 
#include "NrField.h"

#include "XMath.h"
#include "Geodesy.h"

#include "OceanLoad3D.h"

#include "PreloopFFTW.h"
#include <cfloat>

Quad::Quad(const ExodusModel &exModel, int quadTag, const NrField &nrf): 
mQuadTag(quadTag) {
    
    // connectivity and coords
    const IMatX4 &connect = exModel.getConnectivity();
    for (int i = 0; i < 4; i++) {
        int nodeTag = connect(mQuadTag, i);
        mNodalCoords(0, i) = exModel.getNodalS(nodeTag);
        mNodalCoords(1, i) = exModel.getNodalZ(nodeTag);
        mGlobalNodeTags(i) = nodeTag;
    }
    
    // geometric mapping
    double distTol = exModel.getDistTolerance();
    double etype = exModel.getElementalVariables("element_type", mQuadTag);
    if (etype < .5) {
        // etype = 0.0, spherical
        mMapping = new SphericalMapping();
        double r0 = mNodalCoords.col(0).norm();
        double r1 = mNodalCoords.col(1).norm();
        double r2 = mNodalCoords.col(2).norm();
        double r3 = mNodalCoords.col(3).norm();
        if (std::abs(r0 - r1) < distTol && std::abs(r2 - r3) < distTol) {
            mCurvedOuter = r2 > r0 ? 2 : 0;
        } else if (std::abs(r1 - r2) < distTol && std::abs(r3 - r0) < distTol) {
            mCurvedOuter = r3 > r1 ? 3 : 1;       
        } else {
            throw std::runtime_error("Quad::Quad || Invalid spherical element shape.");
        }
    } else if (etype < 1.5) {
        // etype = 1.0, linear
        mMapping = new LinearMapping();
        mCurvedOuter = -1;
    } else {
        // etype = 2.0, semi-spherical
        mMapping = new SemiSphericalMapping();
        double r0 = mNodalCoords.col(0).norm();
        double r1 = mNodalCoords.col(1).norm();
        double r2 = mNodalCoords.col(2).norm();
        double r3 = mNodalCoords.col(3).norm();
        if (std::abs(r0 - r1) < distTol && r0 > r2) {
            mCurvedOuter = 0;
        } else if (std::abs(r1 - r2) < distTol && r1 > r3) {
            mCurvedOuter = 1;
        } else if (std::abs(r2 - r3) < distTol && r2 > r0) {
            mCurvedOuter = 2;
        } else if (std::abs(r3 - r0) < distTol && r3 > r1) {
            mCurvedOuter = 3;
        } else {
            throw std::runtime_error("Quad::Quad || Invalid semi-spherical element shape.");
        }
    }
    
    // solid fluid
    mIsFluid = exModel.getElementalVariables("fluid", mQuadTag) > .5;
    
    // axial boundary
    mIsAxial = false;
    mAxialSide = exModel.getSideAxis(mQuadTag);
    if (mAxialSide >= 0) {
        mIsAxial = true;
        if (mMapping->getType() != Mapping::MappingTypes::Linear && (mAxialSide + mCurvedOuter) % 2 == 0) {
            throw std::runtime_error("Quad::Quad || Conflict in axial setting.");
        }
        if (mAxialSide != 3) {
            throw std::runtime_error("Quad::Quad || Axial side must be 3.");
        }
    }

    mOnSFBoundary = false;
    mSFSide = exModel.getSideSolidFluid(mQuadTag);
    if (mSFSide >= 0) {
        mOnSFBoundary = true;
        if (mMapping->getType() == Mapping::MappingTypes::Linear && !exModel.isCartesian()) {
            throw std::runtime_error("Quad::Quad || Conflict in solid-fluid boundary.");
        }
        if (mMapping->getType() == Mapping::MappingTypes::Spherical && (mSFSide + mCurvedOuter) % 2 != 0) {
            throw std::runtime_error("Quad::Quad || Conflict in solid-fluid boundary.");
        }
        if (mMapping->getType() == Mapping::MappingTypes::SemiSpherical && mSFSide != mCurvedOuter) {
            throw std::runtime_error("Quad::Quad || Conflict in solid-fluid boundary.");
        }
    }

    // surface boundary
    mOnSurface = false;
    mSurfaceSide = exModel.getSideSurface(mQuadTag);
    if (mSurfaceSide >= 0) {
        mOnSurface = true;
        if (mMapping->getType() == Mapping::MappingTypes::Linear && !exModel.isCartesian()) {
            throw std::runtime_error("Quad::Quad || Conflict in surface setting.");
        }
        if (mMapping->getType() == Mapping::MappingTypes::Spherical && mSurfaceSide != mCurvedOuter) {
            throw std::runtime_error("Quad::Quad || Conflict in surface setting.");
        }
        if (mMapping->getType() == Mapping::MappingTypes::SemiSpherical && mSurfaceSide != mCurvedOuter) {
            throw std::runtime_error("Quad::Quad || Conflict in surface setting.");
        }
    }
    
    // nr field 
    for (int j = 0; j < 4; j++) {
        mNearAxisNodes(j) = exModel.getVicinalAxis(quadTag)[j];
        mNodalAveGLLSpacing(j) = exModel.getAveGLLSpacing(mGlobalNodeTags(j));
    } 
    formNrField(nrf, distTol);
    
    // integral factor
    formIntegralFactor();    
    
    // 1D material
    mMaterial = new Material(this, exModel);
    
    // relabelling
    mRelabelling = 0;
    
    // dt 
    mDeltaTRef = exModel.getElementalVariables("dt", mQuadTag);
    
    // init ocean depth
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            mOceanDepth[ipnt] = RDColX::Zero(mPointNr(ipol, jpol));    
        }
    }

    // absorbing boundaries
    mIsSpongeQuad = exModel.isSpongeQuad(mQuadTag);
    
    // set Copycoords for 3D model 
    if (mIsSpongeQuad & exModel.hasParallelModelExtension()) {
        RDCol2 SpongeStartCoords = exModel.getSpongeStartCoords();
        int ABPosition = exModel.getABPosition(mQuadTag);
        
        RDMat24 nodalCoordsConverted = mNodalCoords; // s, z
        
        if (!exModel.isCartesian()) {
            for (int i = 0; i < 4; i++) {
                RDCol2 sz = nodalCoordsConverted.col(i);
                nodalCoordsConverted.col(i) = Geodesy::rtheta(sz);
            }
            nodalCoordsConverted.colwise().reverse(); // theta, r
            SpongeStartCoords.colwise().reverse();
        }
            
        if (ABPosition == 1) {
            nodalCoordsConverted.row(0) = RDRow4::Constant(SpongeStartCoords(0)); // s or theta
        } else if (ABPosition == 2) {
            nodalCoordsConverted.row(1) = RDRow4::Constant(SpongeStartCoords(1)); // z or r
        } else if (ABPosition == 3) {
            nodalCoordsConverted.row(0) = RDRow4::Constant(SpongeStartCoords(0));
            nodalCoordsConverted.row(1) = RDRow4::Constant(SpongeStartCoords(1));
        }
        
        if (!exModel.isCartesian()) {
            for (int i = 0; i < 4; i++) {
                RDCol2 sz;
                sz(0) = nodalCoordsConverted(1, i) * sin(nodalCoordsConverted(0, i));
                sz(1) = nodalCoordsConverted(1, i) * cos(nodalCoordsConverted(0, i));
                nodalCoordsConverted.col(i) = sz;
            }
        }
        
        mCopyCoords = nodalCoordsConverted;
    } else {
        mCopyCoords = mNodalCoords;
    }
    
    mABCRightSide = exModel.getSideRightB(mQuadTag);
    mABCLowerSide = exModel.getSideLowerB(mQuadTag);
    mOnAbsBoundary = exModel.hasStaceyABC() && ((mABCRightSide >= 0) || (mABCLowerSide >= 0));
}

Quad::~Quad() {
    delete mMapping;
    delete mMaterial;
    if (mRelabelling) {
        delete mRelabelling;
    }
}

void Quad::addVolumetric3D(const std::vector<Volumetric3D *> &m3D,
    double srcLat, double srcLon, double srcDep, double phi2D) {
    mMaterial->addVolumetric3D(m3D, srcLat, srcLon, srcDep, phi2D);
}

void Quad::addGeometric3D(const std::vector<Geometric3D *> &g3D, 
    double srcLat, double srcLon, double srcDep, double phi2D) {
    if (g3D.size() > 0) {
        mRelabelling = new Relabelling(this);
        mRelabelling->addUndulation(g3D, srcLat, srcLon, srcDep, phi2D);
    }    
}

void Quad::setOceanLoad3D(const OceanLoad3D &o3D, 
    double srcLat, double srcLon, double srcDep, double phi2D) {
    if (!mOnSurface) {
        return;
    }
    if (mIsFluid) {
        throw std::runtime_error("Quad::setOceanLoad3D || Adding ocean load to fluid surface.");
    }
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            bool surface = mOnSurface && (
                (mSurfaceSide == 0 && jpol == 0) || 
                (mSurfaceSide == 1 && ipol == nPol) || 
                (mSurfaceSide == 2 && jpol == nPol) ||
                (mSurfaceSide == 3 && ipol == 0));
            if (surface) {
                int nr_read = getPointNr(ipol, jpol);
                int ipnt = ipol * nPntEdge + jpol;
                const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
                const RDMatX3 &rtpS = computeGeocentricGlobal(srcLat, srcLon, srcDep, xieta, nr_read, phi2D);
                for (int alpha = 0; alpha < nr_read; alpha++) {
                    double t = rtpS(alpha, 1);
                    double p = rtpS(alpha, 2);
                    mOceanDepth[ipnt](alpha) = o3D.getOceanDepth(t, p); 
                }
            }    
        }
    }
}

bool Quad::hasRelabelling() const {
    if (!mRelabelling) {
        return false;
    }
    return !(mRelabelling->isZero());
}

double Quad::getDeltaT(bool isCartesian) const {
    // courant number
    double courant = getCourant();
    
    // vmax on slices
    RDColX vmax = mMaterial->getVMax();
    
    // hmin on slices
    RDColX hmin = getHminSlices(isCartesian);
    
    // dt on slices
    RDColX dt_slices = courant * hmin.schur(vmax.array().pow(-1.).matrix());
    double dt_min_org = dt_slices.minCoeff();
    
    // we have checked this in relabelling
    return dt_min_org;
}

void Quad::setupGLLPoints(std::vector<GLLPoint *> &gllPoints, const IMatPP &myPointTags,
    double distTol, RDCol2 &Vref_range, RDCol2 &U0_range, const ABCParameters *ABCPar) {
    // compute mass on points
    const arPP_RDColX &mass = mMaterial->computeElementalMass();

    // local properties for absorbing boundaries
    RDMatXN Vp3D = RDMatXN::Zero(mNr, nPE);
    RDMatXN Vs3D = RDMatXN::Zero(mNr, nPE);
    RDMatXN Rho3D = RDMatXN::Zero(mNr, nPE);

    if (mOnAbsBoundary || mIsSpongeQuad) {
        Vp3D = mMaterial->getProperty("vp", -1);
        Vs3D = mMaterial->getProperty("vs", -1);
        Rho3D = mMaterial->getProperty("rho", -1);
    }
    
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            int ipnt = ipol * nPntEdge + jpol;
            // get point tag
            int pointTag = myPointTags(ipol, jpol);
            
            // setup properties
            bool axial = mIsAxial && ipol == 0;
            bool surface = mOnSurface && (
                (mSurfaceSide == 0 && jpol == 0) || 
                (mSurfaceSide == 1 && ipol == nPol) || 
                (mSurfaceSide == 2 && jpol == nPol) ||
                (mSurfaceSide == 3 && ipol == 0));
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            const RDCol2 &crds = mapping(xieta);
            gllPoints[pointTag]->setup(mPointNr(ipol, jpol), axial, surface, crds, distTol);
            
            // add mass 
            if (mIsFluid) {
                gllPoints[pointTag]->addMassFluid(mass[ipnt]);
            } else {
                gllPoints[pointTag]->addMassSolid(mass[ipnt]);
            }
            
            //////// boundary terms ////////
            bool sfbry = mOnSFBoundary && (
                (mSFSide == 0 && jpol == 0) || 
                (mSFSide == 1 && ipol == nPol) || 
                (mSFSide == 2 && jpol == nPol) ||
                (mSFSide == 3 && ipol == 0));
            
            // normal
            if (sfbry) {
                RDMatX3 normal = computeNormal(mSFSide, ipol, jpol);
                // inverse in solid domain
                if (!mIsFluid) {
                    normal *= -1.;
                }
                // both solid and fluid elements contribute, so divided by 2
                gllPoints[pointTag]->addSFNormal(normal * .5);
            }
            
            if (surface) {
                // surface area
                const RDMatX3 &normal = computeNormal(mSurfaceSide, ipol, jpol);
                gllPoints[pointTag]->addSurfNormal(normal);
                gllPoints[pointTag]->setOceanDepth(mOceanDepth[ipnt]);    
            }

            if (mIsSpongeQuad || mOnAbsBoundary) {
                RDColX gamma = RDColX::Zero(mPointNr(ipol, jpol), 1);
                RDMatX3 ABCnormal = RDMatX3::Zero(mPointNr(ipol, jpol),3);
                
                // Stacey Condition
                bool rbry = mOnAbsBoundary && (
                    (mABCRightSide == 0 && jpol == 0) || 
                    (mABCRightSide == 1 && ipol == nPol) || 
                    (mABCRightSide == 2 && jpol == nPol) ||
                    (mABCRightSide == 3 && ipol == 0));
                                 
                bool lbry = mOnAbsBoundary && (
                    (mABCLowerSide == 0 && jpol == 0) || 
                    (mABCLowerSide == 1 && ipol == nPol) || 
                    (mABCLowerSide == 2 && jpol == nPol) ||
                    (mABCLowerSide == 3 && ipol == 0));
                
                if (!rbry && !lbry && !mIsSpongeQuad) {
                    continue;
                }
                
                if (rbry && !lbry) {
                    ABCnormal = computeNormal(mABCRightSide, ipol, jpol);
                } else if (lbry) {
                    ABCnormal = computeNormal(mABCLowerSide, ipol, jpol);
                }
                
                // Kosloff&Kosloff Sponge Boundary
                if (mIsSpongeQuad) {
                    
                    RDCol2 convCrds;
                    if (!ABCPar->isCartesian) {
                        convCrds = Geodesy::rtheta(crds);
                    } else {
                        convCrds = crds;
                    }
                    
                    double sr_frac = (ABCPar->outerCorner[0] - convCrds[0]) / (ABCPar->outerCorner[0] - ABCPar->innerCorner[0]);
                    double ztheta_frac = (ABCPar->outerCorner[1] - convCrds[1]) / (ABCPar->outerCorner[1] - ABCPar->innerCorner[1]);
                    double frac = std::min({sr_frac,ztheta_frac});

                    if (rbry && frac > tinyDouble) {
                        throw std::runtime_error("Quad::setupGLLPoints || Right ABC - mesh inconsistency.");
                    }
                    if (lbry && frac > tinyDouble) {
                        throw std::runtime_error("Quad::setupGLLPoints || Lower ABC - mesh inconsistency.");
                    }
                    
                    RDColX U0;
                    if (ABCPar->U_type == 0) {
                        U0 = RDColX::Constant(mNr, 1, ABCPar->U);
                    } else if (ABCPar->U_type == 1) { // empirical
                        
                        RDColX vs, vp; // get 3D velocities
                        if (Vp3D.rows() > 1) {
                            vp = Vp3D.col(ipnt);
                            vs = Vs3D.col(ipnt);
                        } else {
                            vp = RDColX::Constant(mNr, 1, Vp3D(0, ipnt));
                            vs = RDColX::Constant(mNr, 1, Vs3D(0, ipnt));
                        }
                        
                        double width;
                        if ((ztheta_frac > sr_frac)) {
                            width = ABCPar->width[0];
                        } else {
                            width = ABCPar->width[1];
                            if (!ABCPar->isCartesian) { // handle width given in degrees
                                width *= convCrds[0];
                            }
                        }
                        
                        if (mIsFluid) { // calculate U0
                            U0 = 1.76 * exp(-0.08 * width / (vp.array() * ABCPar->T)) 
                                      / ABCPar->T;
                        } else {
                            U0 = 2.20 * exp(-0.08 * width / (vp.array() * ABCPar->T)) 
                                      * (vs.array() / vp.array()) * (vs.array() / vp.array())
                                      / ABCPar->T;
                        }
                    } else if (ABCPar->U_type == 2) {
                        bool set = false;
                        double Uz;
                        for (int i = 0; i < ABCPar->Us.size(); i+=2) {
                           if (crds[1] > ABCPar->depths[i]) {
                               Uz = ABCPar->Us[i];
                               set = true;
                               break;
                           }
                        }
                        if (set == false) Uz = ABCPar->Us.back();
                        U0 = RDColX::Constant(mNr, 1, Uz);
                    }
                    
                    gamma = U0 * (1 - sin(0.5 * pi * frac) * sin(0.5 * pi * frac));

                    Vref_range[0]=std::min({Vref_range[0], Vp3D.col(ipnt).minCoeff()});
                    Vref_range[1]=std::max({Vref_range[1], Vp3D.col(ipnt).maxCoeff()});
                    U0_range[0]=std::min({U0_range[0], U0.minCoeff()});
                    U0_range[1]=std::max({U0_range[1], U0.maxCoeff()});
                }
            
                gllPoints[pointTag]->setABC(gamma, Rho3D.col(ipnt), Vp3D.col(ipnt), Vs3D.col(ipnt));
                gllPoints[pointTag]->addABCNormal(ABCnormal);
            }
        }
    }
}

int Quad::release(Domain &domain, const IMatPP &myPointTags, const AttBuilder *attBuild, bool recordingWF) const {
    if (mIsFluid) {
        return releaseFluid(domain, myPointTags, recordingWF);
    } else {
        return releaseSolid(domain, myPointTags, attBuild, recordingWF);
    }
}

int Quad::releaseSolid(Domain &domain, const IMatPP &myPointTags, const AttBuilder *attBuild, bool recordingWF) const {
    bool elem1D = mMaterial->isSolidPar1D(attBuild != 0);
    if (hasRelabelling()) {
        elem1D = elem1D && mRelabelling->isPar1D();
    }
    std::array<Point *, nPntElem> points;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            points[ipol * nPntEdge + jpol] = domain.getPoint(myPointTags(ipol, jpol));
        }
    }
    Gradient *grad = createGraident();
    PRT *prt = mRelabelling ? mRelabelling->createPRT(elem1D) : 0;
    Elastic *elas = mMaterial->createElastic(elem1D, attBuild);
    Element *elem = new SolidElement(grad, prt, points, elas);
    if (recordingWF) {
        RDMatPP WF_interpfact;
        getWFRweights(WF_interpfact);
        RDMatXN dr = mRelabelling ? mRelabelling->getDeltaR() : RDMatXN::Zero(mNr, nPntElem);
        elem->addWFRweights(WF_interpfact, mapping(RDCol2::Zero()), dr);
    }
    return domain.addElement(elem);
}

int Quad::releaseFluid(Domain &domain, const IMatPP &myPointTags, bool recordingWF) const {
    bool elem1D = mMaterial->isFluidPar1D();
    if (hasRelabelling()) {
        elem1D = elem1D && mRelabelling->isPar1D();
    }
    std::array<Point *, nPntElem> points;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            points[ipol * nPntEdge + jpol] = domain.getPoint(myPointTags(ipol, jpol));
        }
    }
    Gradient *grad = createGraident();
    PRT *prt = mRelabelling ? mRelabelling->createPRT(elem1D) : 0;
    Acoustic *acous = mMaterial->createAcoustic(elem1D);
    Element *elem = new FluidElement(grad, prt, points, acous);
    if (recordingWF) {
        RDMatPP WF_interpfact;
        getWFRweights(WF_interpfact);
        RDMatXN dr = mRelabelling ? mRelabelling->getDeltaR() : RDMatXN::Zero(mNr, nPntElem);
        elem->addWFRweights(WF_interpfact, mapping(RDCol2::Zero()), dr);
    }
    return domain.addElement(elem);
}

RDCol2 Quad::mapping(const RDCol2 &xieta) const {
    return mMapping->mapping(mNodalCoords, xieta, mCurvedOuter);
}

RDCol2 Quad::mapping(const RDCol2 &xieta, const RDMat24 nodalcoords) const {
    return mMapping->mapping(nodalcoords, xieta, mCurvedOuter);
}

RDMat22 Quad::jacobian(const RDCol2 &xieta) const {
    return mMapping->jacobian(mNodalCoords, xieta, mCurvedOuter);
}

double Quad::detJacobian(const RDCol2 &xieta) const {
    return mMapping->detJacobian(mNodalCoords, xieta, mCurvedOuter);
}

bool Quad::invMapping(const RDCol2 &sz, RDCol2 &xieta) const {
    return mMapping->invMapping(mNodalCoords, sz, mCurvedOuter, xieta);
}

RDRow4 Quad::computeWeightsCG4() const {
    RDRow4 weights_cg4;
    RDMatPP ifact;
    XMath::structuredUseFirstRow(mIntegralFactor, ifact);
    weights_cg4(0) = (ifact(0, 0) + ifact(0, 1)
                   + ifact(1, 0) + ifact(1, 1)
            + 0.5 * (ifact(0, 2) + ifact(1, 2)
                   + ifact(2, 0) + ifact(2, 1))
            + 0.25 * ifact(2, 2)) / ifact(1, 1);
            
    weights_cg4(1) = (ifact(0, 3) + ifact(0, 4)
                   + ifact(1, 3) + ifact(1, 4)
            + 0.5 * (ifact(0, 2) + ifact(1, 2)
                   + ifact(2, 3) + ifact(2, 4))
            + 0.25 * ifact(2, 2)) / ifact(1, 3);
            
    weights_cg4(2) = (ifact(3, 0) + ifact(3, 1)
                   + ifact(4, 0) + ifact(4, 1)
            + 0.5 * (ifact(2, 0) + ifact(2, 1)
                   + ifact(3, 2) + ifact(4, 2))
            + 0.25 * ifact(2, 2)) / ifact(3, 1);
            
    weights_cg4(3) = (ifact(3, 3) + ifact(3, 4)
                   + ifact(4, 3) + ifact(4, 4)
            + 0.5 * (ifact(2, 3) + ifact(2, 4)
                   + ifact(3, 2) + ifact(4, 2))
            + 0.25 * ifact(2, 2)) / ifact(3, 3);            
    return weights_cg4;
}

RDMatX3 Quad::computeGeocentricGlobal(double srcLat, double srcLon, double srcDep,
    const RDCol2 &xieta, int npnt, double phi2D) const {
    RDMatX3 rtpG_Nr(npnt, 3);
    RDCol3 rtpS;
    Geodesy::rtheta(mapping(xieta, mCopyCoords), rtpS(0), rtpS(1));
    double dphi = 2. * pi / npnt;
    for (int i = 0; i < npnt; i++) {
        rtpS(2) = phi2D < -DBL_MAX * .9 ? dphi * i : phi2D;
        rtpG_Nr.row(i) = Geodesy::rotateSrc2Glob(rtpS, srcLat, srcLon, srcDep).transpose();
    }
    return rtpG_Nr;
}

RDMatX3 Quad::computeCartesian(const RDCol2 &xieta, int npnt, double phi2D) const {
    RDMatX3 zsp_Nr(npnt, 3);
    RDCol2 sz = mapping(xieta, mCopyCoords);
    double dphi = 2. * pi / npnt;
    for (int i = 0; i < npnt; i++) {
        zsp_Nr(i,0) = sz(1);
        zsp_Nr(i,1) = sz(0);
        zsp_Nr(i,2) = phi2D < -DBL_MAX * .9 ? dphi * i : phi2D;
    }
    return zsp_Nr;
}

double Quad::computeCenterRadius(bool cartesian) const {
    // xi = eta = 0
    RDCol2 crds = mapping(RDCol2::Zero(), mCopyCoords);
    if (cartesian) {
        return crds(1);
    } else { 
        return crds.norm();
    }
}

void Quad::computeGradientScalar(const vec_CDMatPP &u, vec_ar3_CDMatPP &u_i) const {
    // create Gradient object
    RDMatPP dsdxii, dsdeta, dzdxii, dzdeta, inv_s;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            const RDMat22 &J = jacobian(xieta);
            double detJ = J.determinant();
            dsdxii(ipol, jpol) = J(0, 0) / detJ;
            dsdeta(ipol, jpol) = -J(0, 1) / detJ;
            dzdxii(ipol, jpol) = -J(1, 0) / detJ;
            dzdeta(ipol, jpol) = J(1, 1) / detJ;
            double s = mapping(xieta)(0);
            inv_s(ipol, jpol) = (mIsAxial && ipol == 0) ? 0. : 1. / s;
        }
    }
    PreloopGradient grad(dsdxii, dsdeta, dzdxii, dzdeta, inv_s, mIsAxial);
    // compute gradient
    grad.gradScalar(u, u_i, mNr / 2, mNr % 2 == 0);
}

Gradient *Quad::createGraident() const {
    RDMatPP dsdxii, dsdeta, dzdxii, dzdeta, inv_s;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            const RDMat22 &J = jacobian(xieta);
            double detJ = J.determinant();
            dsdxii(ipol, jpol) = J(0, 0) / detJ;
            dsdeta(ipol, jpol) = -J(0, 1) / detJ;
            dzdxii(ipol, jpol) = -J(1, 0) / detJ;
            dzdeta(ipol, jpol) = J(1, 1) / detJ;
            double s = mapping(xieta)(0);
            inv_s(ipol, jpol) = (mIsAxial && ipol == 0) ? 0. : 1. / s;
        }
    }
    return new Gradient(dsdxii, dsdeta, dzdxii, dzdeta, inv_s, mIsAxial);  
}

void Quad::formIntegralFactor() {
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            const RDCol2 &weights = SpectralConstants::getWeights(ipol, jpol, mIsAxial);
            double wxi_weta = weights(0) * weights(1);
            double detJ = detJacobian(xieta);
            double s = mapping(xieta)(0);
            if (mIsAxial) {
                if (ipol == 0) {
                    const RDMat22 &J = jacobian(xieta);
                    mIntegralFactor(ipol * nPntEdge + jpol) = wxi_weta * J(0, 0) * detJ; 
                } else {
                    mIntegralFactor(ipol * nPntEdge + jpol) = wxi_weta * s / (1. + xieta(0)) * detJ; 
                }
            } else {
                mIntegralFactor(ipol * nPntEdge + jpol) = wxi_weta * s * detJ;  
            }
        }
    }
}

void Quad::formNrField(const NrField &nrf, double distTol) {
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            // interpolate
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            RDCol2 crds = mapping(xieta);
            // erase numerical error
            double offset = std::min((double)tinySingle, distTol / 1e3);
            crds(0) = round(crds(0) / offset) * offset;
            crds(1) = round(crds(1) / offset) * offset;
            mPointNr(ipol, jpol) = nrf.getNrAtPoint(crds);
            
            // upper limit
            double spacing = Mapping::interpolate(mNodalAveGLLSpacing, xieta);
            double circ = 2. * pi * crds(0);
            int upper = (int)(circ / spacing);
            if (upper < 3) upper = 3; // axis
            mPointNr(ipol, jpol) = std::min(mPointNr(ipol, jpol), upper);
            
            ////////// deal with even numbers ////////
            bool forceOdd = false;
            if (mPointNr(ipol, jpol) % 2 == 0) {
                if (mIsAxial) {
                    // axis
                    mPointNr(ipol, jpol)++;
                    forceOdd = true;
                } else if (mNearAxisNodes.maxCoeff() >= 0) {
                    // near axis 
                    for (int i = 0; i < 4; i++) {
                        int n0 = mNearAxisNodes(i);
                        int n1 = mNearAxisNodes(Mapping::period0123(i + 1));
                        if (n0 >= 0 && n1 >= 0) {
                            // on side
                            bool onSide = 
                            (n0 == 0 && jpol == 0) || 
                            (n0 == 1 && ipol == nPol) || 
                            (n0 == 2 && jpol == nPol) ||
                            (n0 == 3 && ipol == 0);
                            if (onSide) {
                                mPointNr(ipol, jpol)++;
                                forceOdd = true;
                                break;
                            }
                        } else if (n0 >= 0 && n1 < 0) {
                            // on point
                            bool onPoint = 
                            (n0 == 0 && ipol == 0 && jpol == 0) || 
                            (n0 == 1 && ipol == nPol && jpol == 0) || 
                            (n0 == 2 && ipol == nPol && jpol == nPol) ||
                            (n0 == 3 && ipol == 0 && jpol == nPol);
                            if (onPoint) {
                                mPointNr(ipol, jpol)++;
                                forceOdd = true;
                                break;
                            }
                        }
                    }
                }
            }
            
            // luck number
            if (nrf.useLuckyNumber()) {
                mPointNr(ipol, jpol) = PreloopFFTW::nextLuckyNumber(mPointNr(ipol, jpol), forceOdd);
            }
        }
    }
    mNr = mPointNr.maxCoeff();
}

double Quad::getCourant() const {
    // 1D reference
    double vmaxRef = mMaterial->getVMaxRef();
    double hminRef = DBL_MAX;
    for (int i = 0; i < 4; i++) {
        double sideLen = (mNodalCoords.col(i) - mNodalCoords.col(Mapping::period0123(i + 1))).norm();
        hminRef = std::min(hminRef, sideLen);
    }
    // dt = courant * hmin / vmax
    return mDeltaTRef * vmaxRef / hminRef;
}

RDColX Quad::getHminSlices(bool isCartesian) const {
    // point coordinates
    std::vector<RDCol2> coordsRef;
    for (int ipol = 0; ipol <= nPol; ipol++) {
        for (int jpol = 0; jpol <= nPol; jpol++) {
            const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
            coordsRef.push_back(mapping(xieta));
        }
    }
    // hmin of reference mesh
    int nslices = getNr();
    RDColX hmin = RDColX::Constant(nslices, XMath::findClosestDist(coordsRef));
    // hmin of deformed mesh
    if (hasRelabelling()) {
        const RDMatXN &deltaR = mRelabelling->getDeltaR();
        for (int islice = 0; islice < nslices; islice++) {
            std::vector<RDCol2> coordsSlice = coordsRef;
            for (int ipol = 0; ipol <= nPol; ipol++) {
                for (int jpol = 0; jpol <= nPol; jpol++) {
                    int ipnt = ipol * nPntEdge + jpol;
                    if (isCartesian) {
                        coordsSlice[ipnt](1) += deltaR(islice, ipnt);
                    } else {
                        double theta = Geodesy::theta(coordsRef[ipnt]);
                        coordsSlice[ipnt](0) += deltaR(islice, ipnt) * sin(theta);
                        coordsSlice[ipnt](1) += deltaR(islice, ipnt) * cos(theta);
                    }
                    
                }
            }
            hmin(islice) = XMath::findClosestDist(coordsSlice);
        }
    }
    return hmin;
}

RDMatX3 Quad::computeNormal(int side, int ipol, int jpol) const {
    RDMatX3 normal;
    if (mMapping->getType() == Mapping::MappingTypes::Spherical && (side + mCurvedOuter) % 2 == 0) {
        normal = computeNormalSpherical(side, ipol, jpol);
    } else if (mMapping->getType() == Mapping::MappingTypes::SemiSpherical && side == mCurvedOuter) {
        normal = computeNormalSpherical(side, ipol, jpol);
    } else {
        normal = computeNormalLinear(side, ipol, jpol);
    }
    return normal;
}

RDMatX3 Quad::computeNormalLinear(int side, int ipol, int jpol) const {
    RDMatX3 normal(mPointNr(ipol, jpol), 3);
    
    const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
    const RDCol2 &crds = mapping(xieta);

    const RDCol2 &sz0 = mNodalCoords.col(side);
    const RDCol2 &sz1 = mNodalCoords.col(Mapping::period0123(side + 1));
    double r0, r1, theta0, theta1;
    Geodesy::rtheta(sz0, r0, theta0);
    Geodesy::rtheta(sz1, r1, theta1);
    
    double rsf = .5 * (r0 + r1); // average radius
    double dl = rsf * std::abs(theta1 - theta0); // average side length in theta-direction

    double sint = crds(0) / crds.norm();
    double cost = crds(1) / crds.norm();
    double s = rsf * sint; // scaling in phi-direction
    
    if (abs(sz0(0) - sz1(0)) < tinyDouble) { // cartesian side normal
        normal.col(0).fill(sz1(1)-sz0(1));//(sz1(1)-sz0(1));
        normal.col(1).fill(0.);
        normal.col(2).fill(0.);
        
        if (hasRelabelling()) {
            normal.col(0) *= mRelabelling->getStaceyNormalWeight(ipol, jpol, side);
        }
    } else if (abs(sz0(1) - sz1(1)) < tinyDouble) { // cartesian top/bottom normal
        if (hasRelabelling()) {
            normal = (sz0(0) - sz1(0)) * mRelabelling->getSFNormalRTZ(ipol, jpol);
        } else {
            normal.col(0).fill(0.);
            normal.col(1).fill(0.);
            normal.col(2).fill(sz1(0)-sz0(0));
        }
    } else if (abs(theta1 - theta0) < tinyDouble) { // non-cartesian side normal
        normal.col(0).fill(sz0(1)-sz1(1));
        normal.col(1).fill(0.);
        normal.col(2).fill(sz1(0)-sz0(0));
        if (hasRelabelling()) {
            normal.col(0) *= mRelabelling->getStaceyNormalWeight(ipol, jpol, side);
        }
    } else {
        throw std::runtime_error("Quad::computeNormal || Invalid element shape.");
    }
    
    const RDCol2 &weights = 0.5 * SpectralConstants::getWeights(ipol, jpol, mIsAxial);
    double wsf = (side == 0 || side == 2) ? weights(0) : weights(1);
    if (mIsAxial) {
        if (ipol == 0) {
            const RDMat22 &J = jacobian(xieta);
            normal *= wsf * J(0, 0);
        } else {
            normal *= wsf / (1. + xieta(0)) * crds(0);
        }
    } else {
        normal *= wsf * crds(0);
    }

    return normal;
}

RDMatX3 Quad::computeNormalSpherical(int side, int ipol, int jpol) const {
    const RDCol2 &xieta = SpectralConstants::getXiEta(ipol, jpol, mIsAxial);
    const RDCol2 &crds = mapping(xieta);

    double r0, r1, theta0, theta1;
    Geodesy::rtheta(mNodalCoords.col(side), r0, theta0);
    Geodesy::rtheta(mNodalCoords.col(Mapping::period0123(side + 1)), r1, theta1);
    
    double rsf = .5 * (r0 + r1); // average radius
    double dl = rsf * std::abs(theta1 - theta0); // average side length in theta-direction

    double sint = crds(0) / crds.norm();
    double cost = crds(1) / crds.norm();
    double s = rsf * sint; // scaling in phi-direction
    
    RDMat33 Q = RDMat33::Zero(3, 3);
    Q(0, 0) = cost;
    Q(0, 2) = sint;
    Q(1, 1) = 1.;
    Q(2, 0) = -sint;
    Q(2, 2) = cost;
    
    RDMatX3 nRTZ(mPointNr(ipol, jpol), 3);
    if (hasRelabelling()) {
        nRTZ = mRelabelling->getSFNormalRTZ(ipol, jpol);
    } else {
        nRTZ.col(0).fill(0.);
        nRTZ.col(1).fill(0.);
        nRTZ.col(2).fill(1.);
    }
    
    RDMatX3 normal = nRTZ * Q.transpose();
    const RDCol2 &weights = 0.5 * SpectralConstants::getWeights(ipol, jpol, mIsAxial);
    double wsf = (side == 0 || side == 2) ? weights(0) : weights(1);
    if (mIsAxial) {
        if (ipol == 0) {
            const RDMat22 &J = jacobian(xieta);
            normal *= wsf * J(0, 0) * dl;
        } else {
            normal *= wsf / (1. + xieta(0)) * s * dl;
        }
    } else {
        normal *= wsf * s * dl;
    }
    if (side != mCurvedOuter) {
        normal *= -1.;
    }
    return normal;
}

RDRowN Quad::getUndulationOnSlice(double phi) const {
    if (hasRelabelling()) {
        return XMath::computeFourierAtPhi(mRelabelling->getDeltaR(), phi);
    } else {
        return RDRowN::Zero();
    }
}

RDRowN Quad::getMaterialOnSlice(const std::string &parName, int refType, double phi) const {
    return XMath::computeFourierAtPhi(mMaterial->getProperty(parName, refType), phi);
}

void Quad::getSpatialRange(double &s_max, double &s_min, double &z_max, double &z_min) const {
    s_max = mNodalCoords.row(0).maxCoeff();
    s_min = mNodalCoords.row(0).minCoeff();
    z_max = mNodalCoords.row(1).maxCoeff();
    z_min = mNodalCoords.row(1).minCoeff();
}

bool Quad::nearMe(double s, double z) const {
    if (s > mNodalCoords.row(0).maxCoeff() + tinySingle || s < mNodalCoords.row(0).minCoeff() - tinySingle) {
        return false;
    }
    if (z > mNodalCoords.row(1).maxCoeff() + tinySingle || z < mNodalCoords.row(1).minCoeff() - tinySingle) {
        return false;
    }
    return true;
}

void Quad::getWFRweights(RDMatPP &interpFact) const {
    RDColP interpXi, interpEta;
    XMath::interpLagrange(0, nPntEdge,
        mIsAxial ? SpectralConstants::getP_GLJ().data():
        SpectralConstants::getP_GLL().data(), interpXi.data());
    XMath::interpLagrange(0, nPntEdge,
        SpectralConstants::getP_GLL().data(), interpEta.data());
    interpFact = interpXi * interpEta.transpose();
    if (mIsFluid) {
        for (int ipol = 0; ipol <= nPol; ipol++) {
            for (int jpol = 0; jpol <= nPol; jpol++) {
                int ipnt = ipol * nPntEdge + jpol;
                interpFact(ipol, jpol) /= mIntegralFactor(ipnt);    
            }
        }
    }
}

int Quad::edgeAtRadius(double radius, double distTol, bool upper) const {
    double r0 = mNodalCoords.col(0).norm();
    double r1 = mNodalCoords.col(1).norm();
    double r2 = mNodalCoords.col(2).norm();
    double r3 = mNodalCoords.col(3).norm();
    double rc = (r0 + r1 + r2 + r3) / 4;
    if (upper != (rc > radius)) {
        return -1;
    }
    
    if (std::abs(r0 - radius) < distTol && std::abs(r1 - radius) < distTol) {
        return 0;
    }
    if (std::abs(r1 - radius) < distTol && std::abs(r2 - radius) < distTol) {
        return 1;
    }
    if (std::abs(r2 - radius) < distTol && std::abs(r3 - radius) < distTol) {
        return 2;
    }
    if (std::abs(r3 - radius) < distTol && std::abs(r0 - radius) < distTol) {
        return 3;
    }
    return -1;
}


