// ABCParameters.h
// created by Claudia on 27-Mar-2018
// AxiSEM3D Mesh

#include "eigenp.h"
#include "Parameters.h"
#include <math.h>

class ExodusModel;

class ABCParameters {

public:
    int n;
    double Hmax, width, Ufac;
    RDCol2 boundaries;

    ABCParameters(const Parameters &par, const ExodusModel *exModel) {
    n = exModel->getNumAbsElements();
    Ufac = par.getValue<double>("SPONGE_ABSORBING_FACTOR");
    Hmax = exModel->getHmax();
    boundaries = exModel->getBoundaries();

    double T = par.getValue<double>("SOURCE_STF_HALF_DURATION");
    double v_max = exModel->getABVmax();
    if (Ufac <= 0) {Ufac = 5 * 0.17 * (27 / pow(T, 0.45) + 27.82) / (n + 7);}

    width = n * Hmax;
    };

};
