// Mass.h
// created by Kuangdai on 3-Apr-2016 
// base class of mass

#pragma once

#include "eigenc.h"

class ABC {
public:
    ~ABC() {};
    ABC(const RColX gamma, const RColX v_rho_p, const RColX v_rho_s, const int nr);
    void setNormal(RMatX3 const normal);

    CMatX3 StaceyTraction(const CMatX3 veloc) const;
    CColX StaceyTraction(const CColX veloc) const;

    void applyKosloffDamping(CMatX3 &accel, const CMatX3 veloc, const CMatX3 displ);
    void applyKosloffDamping(CColX &accel, const CColX veloc, const CColX displ);
    
private:    
    int mNr, mNu;
    double mWeight;
    RColX mGamma;
    RMatX3 mNormal; 
    RColX mVp_Rho, mVs_Rho;
    bool mStacey, mKosloff, mStacey3D, mKosloff3D;
};
