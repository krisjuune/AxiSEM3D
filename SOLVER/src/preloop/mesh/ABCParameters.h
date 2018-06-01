// ABCParameters.h
// created by Claudia on 27-Mar-2018
// AxiSEM3D Mesh

#include "eigenp.h"
#include "Parameters.h"

class ExodusModel;

class ABCParameters {

public:
    int n;
    double Hmax, width, Vp_min, Ufac;
    RDCol2 boundaries;

    ABCParameters(const Parameters &par, const ExodusModel *exModel) {
    n = par.getValue<int>("ABC_ELEMENTS");
    Ufac = par.getValue<double>("ABC_FACTOR");

    if (Ufac <= 0) {Ufac = 1200 / (n + 8);}

    Hmax = exModel->getHmax();
    boundaries = exModel->getBoundaries();

    width = n * Hmax;
    };

};
