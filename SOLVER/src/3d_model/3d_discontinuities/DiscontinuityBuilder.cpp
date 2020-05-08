// DiscontinuityBuilder.cpp
// created by Claudia on 2-Mar_2020
// General Interface handling: -> read continuous vertical discontinuities from volumetric 3D models and hande them as geometric EMC
//                             -> effectively flatten 3D models so they are read correctly in the presence of relabelling
//                             -> define geometric EMC models as discontinuities, so there is no interpolation across them when reading volumettric models

#include "Geometric3D.h"
#include "Geometric3D_EMC.h"
#include "ExodusModel.h"
#include "Volumetric3D.h"
#include <sstream>
#include "XMPI.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Parameters.h"
#include "eigenp.h"
#include "DiscontinuityBuilder.h"

int DiscontinuityBuilder::mNGeom = 0;
int DiscontinuityBuilder::mNVol = 0;
int DiscontinuityBuilder::mNGeomFlat = 0;
int DiscontinuityBuilder::mNVolFlat = 0;
int DiscontinuityBuilder::mNGeomDisc = 0;
int DiscontinuityBuilder::mNVolDisc = 0;
    
void DiscontinuityBuilder::buildInparam(std::vector<Volumetric3D *> Vmodels, std::vector<Geometric3D *> &Gmodels, 
    const Parameters &par, const ExodusModel *exModel, int verb) {

    int ndiscs = par.getValue<int>("MODEL_3D_DISCONTINUITIES_NUM");
    int nsize = par.getSize("MODEL_3D_DISCONTINUITIES_LIST");
    if (ndiscs > nsize) {
        throw std::runtime_error("DiscontinuityBuilder::buildInparam || "
            "Not enough model names provided in MODEL_3D_DISCONTINUITIES_LIST, ||"
            "MODEL_3D_DISCONTINUITIES_NUM = " + boost::lexical_cast<std::string>(ndiscs) + 
            ", but only " + boost::lexical_cast<std::string>(nsize) + " provided.");
    }
    
    const std::string source = "DiscontinuityBuilder::buildInparam";
    
    std::vector<Geometric3D *> flattening;
    std::vector<Geometric3D *> discontinuities;
    
    for (int idisc = 0; idisc < ndiscs; idisc++) {
        // split model name and parameters
        std::string dstr = par.getValue<std::string>("MODEL_3D_DISCONTINUITIES_LIST", idisc);
        std::vector<std::string> strs = Parameters::splitString(dstr, "$");
        std::string name(strs[0]);
        
        int imodel;
        bool isDiscontinuity, causesFlattening;
        Parameters::castValue(imodel, strs[1], source);
        Parameters::castValue(isDiscontinuity, strs[2], source);
        
        imodel -= 1;
        
        if (boost::iequals(name, "from_volumetric")) {
            std::vector<std::string> params(strs.begin() + 3, strs.end());
            
            Geometric3D *m;
            m = new Geometric3D_EMC();
            buildFromVolumetric3D(m, Vmodels[imodel], params, exModel, verb);
            Gmodels.push_back(m);
            
            causesFlattening = true;
            flattening.push_back(m);
            mNVolFlat++;
            if (isDiscontinuity) {
                discontinuities.push_back(m);
                mNVolDisc++;
            }
            mNVol++;
        } else if (boost::iequals(name, "from_geometric")) {
            Parameters::castValue(causesFlattening, strs[3], source);
            if (isDiscontinuity) {
                discontinuities.push_back(Gmodels[imodel]);
                mNGeomDisc++;
            }
            if (causesFlattening) {
                flattening.push_back(Gmodels[imodel]);
                mNGeomFlat++;
            }
            mNGeom++;
        } else {
            throw std::runtime_error("DiscontinuityBuilder::buildInparam || "
                "Unknown discontinuity name " + name + ".");
        }
        
   }
   
   for (auto &model: Vmodels) {
       model->setFlattening(flattening);
       model->setDiscontinuities(discontinuities);
   }
   
   if (verb) XMPI::cout << verbose();
}

void DiscontinuityBuilder::buildFromVolumetric3D(Geometric3D *&Gmodel, const Volumetric3D *Vmodel, 
    const std::vector<std::string> params, const ExodusModel *exModel, int verbose) {
    
    const std::string source = "DiscontinuityBuilder::buildFromVolumetric3D";
    
    if (!Vmodel->isEMC()) {
        throw std::runtime_error("DiscontinuityBuilder::buildFromVolumetric3D || "
                "Auto-generation of discontinuities is only implemented for EMC-type models.");
    }
    
    double search_lower = exModel->getROuter();
    double search_upper = 0;
    std::string compType = "greater";
    bool from_bottom = false;
    
    double cutoff, range;
    bool layer_from_mesh, range_absolute;
    
    Parameters::castValue(cutoff, params[0], source);
    Parameters::castValue(layer_from_mesh, params[1], source);
    Parameters::castValue(range, params[2], source);
    Parameters::castValue(range_absolute, params[3], source);    
    
    try {
        int ipar = 4;
        Parameters::castValue(compType, params.at(ipar++), source);
        Parameters::castValue(search_lower, params.at(ipar++), source);
        search_lower *= 1000;
        Parameters::castValue(search_upper, params.at(ipar++), source);
        search_upper *= 1000;
        Parameters::castValue(from_bottom, params.at(ipar++), source);
    } catch (std::out_of_range) {
        // nothing
    }
    
    RDColX lat, lon;
    RDMatXX depth;
    bool cartesian;
    Vmodel->findDiscontinuity(lat, lon, depth, cartesian, cutoff, search_upper, search_lower, compType, from_bottom);
    
    double layer, upper, lower;
    if (layer_from_mesh) {
        bool found = exModel->findDiscontinuity(layer, Vmodel->getVarName(), cutoff, search_upper, 
        search_lower, compType, from_bottom);
        if (!found) {
            throw std::runtime_error("DiscontinuityBuilder::buildFromVolumetric3D || "
                "Discontinuity not found in Exodus Model.");
        }
    } else {
        layer = depth.mean();
        exModel->findClosestMeshLine(layer);
    }
    depth = layer - depth.array();
    
    if (range_absolute) {
        upper = layer - range * 1000;
        lower = layer + range * 1000;
    } else {
        double maxdev = depth.array().abs().maxCoeff();
        upper = layer - range * maxdev;
        lower = layer + range * maxdev;
    }
    upper = std::max(Geodesy::getROuter() - exModel->upperEdge(), upper);
    lower = std::min(Geodesy::getROuter() - exModel->lowerEdge(), lower);
    
    Gmodel->initializeFromVolumetric(layer, upper, lower, lat, lon, depth, cartesian);
    
    double srcLat, srcLon, srcDep;
    Vmodel->getSourceLocation(srcLat, srcLon, srcDep);
    Gmodel->setSourceLocation(srcLat, srcLon, srcDep);
    
    if (verbose) XMPI::cout << Gmodel->verbose();
}


std::string DiscontinuityBuilder::verbose() {
    std::stringstream ss;
    ss << "\n====================== Discontinuities =====================" << std::endl;
    ss << "  Vol. Interfaces     =   "<< mNVol << std::endl;
    ss << "    Discont. from Vol.  =   " << mNVolDisc << std::endl;
    ss << "    Flatten. from Vol.  =   " << mNVolFlat << std::endl;
    ss << "  Geom. Interfaces    =   " << mNGeom << std::endl;
    ss << "    Discont. from Geom. =   " << mNGeomDisc << std::endl;
    ss << "    Flatten. from Geom. =   " << mNGeomFlat << std::endl;
    ss << "====================== Discontinuities =====================\n" << std::endl;
    return ss.str();
}