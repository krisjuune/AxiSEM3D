// DiscontinuityBuilder.h
// created by Claudia on 2-Mar_2020

#include <string>
#include <vector>
#include <eigenp.h>

class Volumetric3D;
class Geometric3D;
class ExodusModel;
class Parameters;

class DiscontinuityBuilder {
public:
    
    static void buildInparam(std::vector<Volumetric3D *> Vmodels, std::vector<Geometric3D *> &Gmodels, 
        const Parameters &par, const ExodusModel *exModel, int verbose);
    
private:
    static void buildFromVolumetric3D(Geometric3D *&Gmodel, const Volumetric3D *Vmodel, 
        const std::vector<std::string> params, const ExodusModel *exModel, int verbose);
    static std::string verbose();
    
    static int mNGeom, mNVol, mNGeomFlat, mNVolFlat, mNGeomDisc, mNVolDisc;
};