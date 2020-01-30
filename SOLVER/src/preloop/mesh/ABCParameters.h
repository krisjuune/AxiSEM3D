// ABCParameters.h
// created by Claudia on 27-Mar-2018
// AxiSEM3D Mesh

#include "eigenp.h"
#include "Parameters.h"
#include <math.h>

class ExodusModel;

class ABCParameters {

public:
    bool isCartesian;
    int n, U_type;
    double U, T;
    RDCol2 outerCorner, innerCorner, width;
    std::vector<double> depths;
    std::vector<double> Us;
    std::string s_type;

    ABCParameters(const Parameters &par, const ExodusModel *exModel) {
    
    isCartesian = exModel->isCartesian();
    T = 2 * par.getValue<double>("SOURCE_STF_HALF_DURATION");

    n = exModel->getNumAbsElements();
    
    innerCorner = exModel->getSpongeStartCoords();
    outerCorner = exModel->getMeshBoundaryCoords();
    
    if (!isCartesian) {
        innerCorner = innerCorner.colwise().reverse();
        outerCorner = outerCorner.colwise().reverse();
    }
    
    width = outerCorner - innerCorner;
    width = width.cwiseAbs();
    
    // sponge parametrisation type
    s_type = par.getValue<std::string>("ABC_SPONGE_BOUNDARIES_TYPE");
    if (boost::iequals(s_type, "constant")) {
        U_type = 0;
        U = par.getValue<double>("ABC_SPONGE_BOUNDARIES_FACTOR");
    } else if (boost::iequals(s_type, "empirical")) {
        U_type = 1;
    } else if (boost::iequals(s_type, "vertical_profile")) {
        U_type = 2;
        std::string mstr = par.getValue<std::string>("ABC_SPONGE_BOUNDARIES_FACTOR");
        std::vector<std::string> strs = Parameters::splitString(mstr, "$");
        for (int i = 0; i < strs.size(); i+=2) {
            Us.push_back(boost::lexical_cast<double>(strs[i]));
            if (strs.size() > i+1) {
                depths.push_back(exModel->getROuter() - 1000 * boost::lexical_cast<double>(strs[i+1]));
            }
        }
    } else {
        throw std::runtime_error("ABCParameters::ABCParameters || "
            "Unknown ABC absorbtion factor format " + s_type + ".");
    }
    
    };

};
