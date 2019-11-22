// SolidPoint.h
// created by Kuangdai on 3-Apr-2016 
// solid gll points 

#pragma once

class Mass;
class ABC;
#include "Point.h"

class SolidPoint: public Point {
    friend class SolidFluidPoint;
    
public:
    SolidPoint(int nr, bool axial, const RDCol2 &crds, Mass *mass, ABC *ABC);
    ~SolidPoint();
    
    // update in time domain by Newmark
    void updateNewmark(double dt);
    
    // check stability
    bool stable() const {return mDispl.allFinite();};
    
    // reset to zero 
    void resetZero(); 
    
    // randomize disp and stiff
    void randomDispl(Real factor = one, int seed = -1, int max_order = -1);
    void randomStiff(Real factor = one, int seed = -1, int max_order = -1);
    
    // verbose
    std::string verbose() const;
    
    // measure cost 
    double measure(int count);
    
    // test mass
    void test();
    
    // communication size
    int sizeComm() const {return mStiff.size();};
    
    // communication
    void feedBuffer(CColX &buffer, int &row);
    void extractBuffer(CColX &buffer, int &row);
    
    ///////////// solid-only /////////////   
    // scatter displ to element
    void scatterDisplToElement(vec_ar3_CMatPP &displ, int ipol, int jpol, int maxNu) const;
    
    // gather stiff from element
    void gatherStiffFromElement(const vec_ar3_CMatPP &stiff, int ipol, int jpol);
    
    // add to stiff, used by source
    void addToStiff(const CMatX3 &source);
    
    // wisdom
    void learnWisdom(Real cutoff, int nPeaks);
    int getNuWisdom() const;
    
    // get displacement
    const CMatX3 &getDispFourierSolid() const {return mDispl;};
    
    void applyABC();
    
private:
    
    // mask 
    void maskField(CMatX3 &field);

    // fields
    CMatX3 mDispl;
    CMatX3 mVeloc;
    CMatX3 mAccel;
    CMatX3 mStiff;
    
    // mass
    Mass *mMass;
    
    // ABC factor
    ABC *mABC = 0;

    // max disp norm for wisddom learning
    CMatX3 mLastDisplWisdom = CMatX3::Zero(mNu + 1, 3);
    RRow3 mLastMaxDisplWisdom = -RRow3::Ones();
    Eigen::Matrix<int, 1, 3> mNuWisdom = Eigen::Matrix<int, 1, 3>::Zero();
    Eigen::Matrix<int, 1, 3> mPeaksWisdom = Eigen::Matrix<int, 1, 3>::Zero();
    Eigen::Matrix<bool, 1, 3> mApproachingPeak = Eigen::Matrix<bool, 1, 3>::Ones();
};

