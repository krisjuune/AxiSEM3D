// Volumetric3D.cpp
// created by Kuangdai on 16-May-2016 
// base class of volumetric 3D models, including 
// velocity, density, elasticity, attenuation

#include "Volumetric3D.h"

#include "XMPI.h"
#include "Parameters.h"
#include "AutoGeometricParams.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

/////////////////////////////// user-defined models here
#include "Volumetric3D_s20rts.h"
#include "Volumetric3D_s40rts.h"
#include "Volumetric3D_crust1.h"
#include "Volumetric3D_bubble.h"
#include "Volumetric3D_cylinder.h"
#include "Volumetric3D_block.h"
#include "Volumetric3D_EMC.h"
#include "Volumetric3D_SEG.h"
/////////////////////////////// user-defined models here

void Volumetric3D::buildInparam(std::vector<Volumetric3D *> &models,
    const Parameters &par, const ExodusModel *exModel, std::vector<AutoGeometricParams *> &Vol2GeoModels,
    double srcLat, double srcLon, double srcDep, int verbose) {
    
    // clear the container
    for (const auto &m: models) {
        delete m;    
    }    
    models.clear();
    
    // check size
    int nmodels = par.getValue<int>("MODEL_3D_VOLUMETRIC_NUM");
    int nsize = par.getSize("MODEL_3D_VOLUMETRIC_LIST");
    if (nmodels > nsize) {
        throw std::runtime_error("Volumetric3D::buildInparam || "
            "Not enough model names provided in MODEL_3D_VOLUMETRIC_LIST, ||"
            "MODEL_3D_VOLUMETRIC_NUM = " + boost::lexical_cast<std::string>(nmodels) + 
            ", but only " + boost::lexical_cast<std::string>(nsize) + " provided.");
    }
    
    for (int imodel = 0; imodel < nmodels; imodel++) {
        
        // split model name and parameters
        std::string mstr = par.getValue<std::string>("MODEL_3D_VOLUMETRIC_LIST", imodel);
        std::vector<std::string> strs = Parameters::splitString(mstr, "$");
        std::string name(strs[0]);
        std::vector<std::string> params(strs.begin() + 1, strs.end());
        
        // create model
        Volumetric3D *m;
        if (boost::iequals(name, "s20rts")) {
            m = new Volumetric3D_s20rts();
        } else if (boost::iequals(name, "s40rts")) {
            m = new Volumetric3D_s40rts();
        } else if (boost::iequals(name, "crust1")) {
            m = new Volumetric3D_crust1();
        } else if (boost::iequals(name, "bubble")) {
            m = new Volumetric3D_bubble(); 
        } else if (boost::iequals(name, "cylinder")) {
            m = new Volumetric3D_cylinder();
        } else if (boost::iequals(name, "block")) {
            m = new Volumetric3D_block();        
        } else if (boost::iequals(name, "emc")) {
            m = new Volumetric3D_EMC();
        } else if (boost::iequals(name, "SEG_C3NA")) {
            m = new Volumetric3D_SEG();

        ///////////////////////////////
        // user-defined models here
        /////////////////////////////// 
            
        } else {
            throw std::runtime_error("Volumetric3D::buildInparam || "
                "Unknown volumetric model name " + name + ".");
        }
        
        // initialize
        m->setSourceLocation(srcLat, srcLon, srcDep);
        m->setupExodusModel(exModel);
        m->initialize(params);
        models.push_back(m);

        // SEG-model autogenerates ocean floor relabelling
        m->modelBathymetry(Vol2GeoModels);

        // verbose
        if (verbose) {
            XMPI::cout << m->verbose();
        }
    }
}
