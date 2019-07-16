// NrField.cpp
// created by Kuangdai on 13-May-2016 
// base class of nr integer field

#include "NrField.h"
#include "ConstNrField.h"
#include "EmpNrField.h"
#include "WisdomNrField.h"
#include "UserNrField.h"
#include "ExodusModel.h"

#include "Parameters.h"
#include "XMPI.h"
#include <boost/algorithm/string.hpp>

bool NrField::mHasLowOrderExt = false;
RDCol2 NrField::mInnerBoundaries = RDCol2::Zero();
int NrField::mExtNu = 0;

void NrField::buildInparam(NrField *&nrf, const Parameters &par, const ExodusModel *exModel, int verbose) {
    if (nrf) {
        delete nrf;
    }
    
    std::string type = par.getValue<std::string>("NU_TYPE");
    bool useLucky = par.getValue<bool>("FFTW_LUCKY_NUMBER");
    
    mHasLowOrderExt = (exModel->hasExtension() & !exModel->hasSpongeABC());
    mInnerBoundaries = exModel->getInnerBoundaries();
    mExtNu = par.getValue<int>("EXTENSION_ORDER");
    
    if (boost::iequals(type, "constant")) {
        int nu = par.getValue<int>("NU_CONST");
        nrf = new ConstNrField(useLucky, nu);
    } else if (boost::iequals(type, "empirical")) {
        int nu_ref = par.getValue<int>("NU_EMP_REF");
        int nu_min = par.getValue<int>("NU_EMP_MIN");
        bool scaleS = par.getValue<bool>("NU_EMP_SCALE_AXIS");
        bool scaleT = par.getValue<bool>("NU_EMP_SCALE_THETA");
        bool scaleD = par.getValue<bool>("NU_EMP_SCALE_DEPTH");
        double powS = par.getValue<double>("NU_EMP_POW_AXIS");
        double factPI = par.getValue<double>("NU_EMP_FACTOR_PI");
        double startT = par.getValue<double>("NU_EMP_THETA_START") * degree;
        double powT = par.getValue<double>("NU_EMP_POW_THETA");
        double factD0 = par.getValue<double>("NU_EMP_FACTOR_SURF");
        double startD = par.getValue<double>("NU_EMP_DEPTH_START") * 1e3;
        double endD = par.getValue<double>("NU_EMP_DEPTH_END") * 1e3;
        nrf = new EmpNrField(useLucky, nu_ref, nu_min, scaleS, scaleT, scaleD, 
            powS, factPI, startT, powT, factD0, startD, endD);
    } else if (boost::iequals(type, "wisdom")) {
        std::string fname = Parameters::sInputDirectory + "/" + par.getValue<std::string>("NU_WISDOM_REUSE_INPUT");
        double factor = par.getValue<double>("NU_WISDOM_REUSE_FACTOR");
        if (factor <= tinyDouble) {
            factor = 1.0;
        }
        nrf = new WisdomNrField(useLucky, fname, factor);
    } else if (boost::iequals(type, "user-defined")) {
        int nsize = par.getSize("NU_USER_PARAMETER_LIST");
        std::vector<double> params;
        for (int ipar = 0; ipar < nsize; ipar++) {
            double param = par.getValue<double>("NU_USER_PARAMETER_LIST", ipar);
            params.push_back(param);
        }
        nrf = new UserNrField(useLucky, params);
    } else {
        throw std::runtime_error("NrField::build || "
            "Invalid parameter, keyword = NU_TYPE.");
    }
    
    if (verbose) {
        XMPI::cout << nrf->verbose();
    }
}

int NrField::getNrAtPoint(const RDCol2 &coords) const {
    int nr = getNrAtPointInternal(coords);
    if (mHasLowOrderExt && (coords(0) > mInnerBoundaries(0) || coords(1) < mInnerBoundaries(1))) {
        nr = std::min({2 * mExtNu + 1, nr});
    }
    return nr;
}