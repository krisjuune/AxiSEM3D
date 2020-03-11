// Geometric3D.cpp
// created by Kuangdai on 4-Jun-2016 
// base class of geometric 3D models, such as 
// topography, ellipticity and undulating Moho

#include "Geometric3D.h"

#include "XMPI.h"
#include "Parameters.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

/////////////////////////////// user-defined models here
#include "Ellipticity.h"
#include "Geometric3D_crust1.h"
#include "Geometric3D_EMC.h"
/////////////////////////////// user-defined models here

void Geometric3D::buildInparam(std::vector<Geometric3D *> &models,
        const Parameters &par, double srcLat, double srcLon, double srcDep, int verbose) {
    // clear the container
    for (const auto &m: models) {
        delete m;    
    }    
    models.clear();
    
    // check size
    int nmodels = par.getValue<int>("MODEL_3D_GEOMETRIC_NUM");
    int nsize = par.getSize("MODEL_3D_GEOMETRIC_LIST");
    if (nmodels > nsize) {
        throw std::runtime_error("Geometric3D::buildInparam || "
            "Not enough model names provided in MODEL_3D_GEOMETRIC_LIST, ||"
            "MODEL_3D_GEOMETRIC_NUM = " + boost::lexical_cast<std::string>(nmodels) + 
            ", but only " + boost::lexical_cast<std::string>(nsize) + " provided.");
    }
    
    int nShiftSize = par.getSize("MODEL_3D_GEOMETRIC_SHIFT");
    
    for (int imodel = 0; imodel < nmodels; imodel++) {
        
        // split model name and parameters
        std::string mstr = par.getValue<std::string>("MODEL_3D_GEOMETRIC_LIST", imodel);
        std::vector<std::string> strs = Parameters::splitString(mstr, "$");
        std::string name(strs[0]);
        std::vector<std::string> params(strs.begin() + 1, strs.end());
        
        // create model
        Geometric3D *m;
        if (boost::iequals(name, "crust1")) {
            m = new Geometric3D_crust1();
            
        } else if (boost::iequals(name, "emc")) {
            m = new Geometric3D_EMC();
            
        /////////////////////////////// 
        // user-defined models here
        /////////////////////////////// 
            
        } else {
            throw std::runtime_error("Geometric3D::buildInparam || "
                "Unknown 3D geometric model: " + name + ".");
        }
        
        // initialize
        m->initialize(params);
        m->setSourceLocation(srcLat, srcLon, srcDep);        
        bool hasShift;
        std::vector<std::string> shiftstrs;
        if (imodel >= nShiftSize) {
            hasShift = false;
        } else {
            std::string shiftstr = par.getValue<std::string>("MODEL_3D_GEOMETRIC_SHIFT", imodel);
            shiftstrs = Parameters::splitString(shiftstr, "$");
            Parameters::castValue(hasShift, shiftstrs[0], "Geometric3D::buildInparam");
        }
        
        if (hasShift) {
            if (!(shiftstrs.size() == 3)) {
                throw std::runtime_error("Geometric3D::buildInparam || "
                    "2 coordinate values required for shifting 3D model.");
            }
            double x, y, z;
            Parameters::castValue(x, shiftstrs[1], "Geometric3D::buildInparam");
            Parameters::castValue(y, shiftstrs[2], "Geometric3D::buildInparam");
            m->applyShift(x, y);
        }
        models.push_back(m);
        
        // verbose
        if (verbose) {
            XMPI::cout << m->verbose();
        }
    }
    
    // ELLIPTICITY
    std::string emode = par.getValue<std::string>("MODEL_3D_ELLIPTICITY_MODE");
    if (boost::iequals(emode, "full")) {
        Geometric3D *m = new Ellipticity();
        models.push_back(m);
        if (verbose) {
            XMPI::cout << m->verbose();
        }
    }
}
