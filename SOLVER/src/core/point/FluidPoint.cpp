// FluidPoint.cpp
// created by Kuangdai on 4-Apr-2016 
// fluid gll points 

#include "FluidPoint.h"
#include "Mass.h"
#include "MultilevelTimer.h"
#include "ABC.h"

FluidPoint::FluidPoint(int nr, bool axial, const RDCol2 &crds, Mass *mass, const bool fluidSurf, 
    ABC *ABC): Point(nr, axial, crds), mMass(mass), mFluidSurf(fluidSurf), mABC(ABC) {
    mDispl = CColX::Zero(mNu + 1, 1);
    mVeloc = CColX::Zero(mNu + 1, 1);
    mAccel = CColX::Zero(mNu + 1, 1);
    mStiff = CColX::Zero(mNu + 1, 1);
    mMass->checkCompatibility(nr);
    //mNuWisdom = mNu;
}

FluidPoint::~FluidPoint() {
    delete mMass;
}

void FluidPoint::updateNewmark(Real dt) {
    // calculate new acceleration

    if(mFluidSurf) {
      resetZero();
      return;
    }

    // mask stiff
    maskField(mStiff);
      // compute accel inplace
    mMass->computeAccel(mStiff);
      // mask accel (masking must be called twice if mass is 3D)
    maskField(mStiff);
    // update dt
    double half_dt = half * dt;
    double half_dt_dt = half_dt * dt;
    
    // update velocity and old acceleration
    mVeloc += (Real)half_dt * (mAccel + mStiff);
    mAccel = mStiff;
    // update displacement
    mDispl += (Real)dt * mVeloc + (Real)half_dt_dt * mAccel;

    // zero stiffness for next time step
    mStiff.setZero();
}

void FluidPoint::applyABC() {
    if (mABC) {
        mStiff -= mABC->StaceyTraction(mVeloc);
        mABC->applyKosloffDamping(mAccel, mVeloc, mDispl);
    }
}

void FluidPoint::resetZero() {
    mStiff.setZero();
    mDispl.setZero();
    mVeloc.setZero();
    mAccel.setZero();
}

void FluidPoint::randomDispl(Real factor, int seed, int max_order) {
    if (seed >= 0) {
        std::srand(seed);
    }
    if (max_order < 0 || max_order > mNu) {
        mDispl.setRandom(); 
    } else {
        mDispl.topRows(max_order + 1).setRandom();
    }
    mDispl *= factor;
    maskField(mDispl);
}

void FluidPoint::randomStiff(Real factor, int seed, int max_order) {
    if (seed >= 0) {
        std::srand(seed);
    }
    if (max_order < 0 || max_order > mNu) {
        mStiff.setRandom(); 
    } else {
        mStiff.topRows(max_order + 1).setRandom();
    }
    mStiff *= factor;
    maskField(mStiff);
}

std::string FluidPoint::verbose() const {
    return "FluidPoint$" + mMass->verbose();
}

double FluidPoint::measure(int count) {
    Real dt = .1;
    randomStiff((Real)1e6); 
    MyBoostTimer timer;
    timer.start();
    for (int i = 0; i < count; i++) {
        maskField(mStiff);
        mMass->computeAccel(mStiff);
        maskField(mStiff);
        Real half_dt = half * dt;
        Real half_dt_dt = half_dt * dt;
        mVeloc += half_dt * (mAccel + mStiff);
        mAccel = mStiff;
        mDispl += dt * mVeloc + half_dt_dt * mAccel; 
        mVeloc.setZero();
    }
    double elapsed_time = timer.elapsed();
    resetZero();
    return elapsed_time / count;
}

void FluidPoint::test() {
    // mass matrix
    int totalDim = mNu + 1;
    RMatXX M = RMatXX::Zero(totalDim, totalDim);
    
    resetZero();
    for (int alpha = 0; alpha <= mNu; alpha++) {
        if (mNr % 2 == 0 && alpha == mNu) {continue;}
        // delta function
        if (axial() && alpha > 0) {continue;}
        mStiff(alpha) = one;
        if (alpha == 0) {mStiff(alpha) = two;}
                
        // compute stiff 
        updateNewmark(1.);
                
        // positive-definite
        Real sr = mAccel(alpha).real();
        if (sr <= zero) {
            // add code here to debug
            throw std::runtime_error("FluidPoint::test || "
                "Mass matrix is not positive definite.");  
        }
                    
        // store mass
        int row = alpha;
        for (int alpha1 = 0; alpha1 <= mNu; alpha1++) {
            if (mNr % 2 == 0 && alpha1 == mNu) {continue;}
            if (axial() && alpha > 0) {continue;}
            int col = alpha1;
            M(row, col) = mAccel(alpha1).real();
        }
                
        // restore zero 
        mStiff(alpha) = czero;
    }
    resetZero();

    // test self-adjointness 
    Real maxM = M.array().abs().maxCoeff();
    Real tole = maxM * tinyReal;
    for (int i = 0; i < totalDim; i++) {
        for (int j = i + 1; j < totalDim; j++) {
            Real diff = std::abs(M(i, j) - M(j, i));
            if (diff > tole) {
                // add code here to debug
                throw std::runtime_error("FluidPoint::test || "
                    "Mass matrix is not self-adjoint."); 
            }
        }
    }
}

void FluidPoint::feedBuffer(CColX &buffer, int &row) {
    int rows = mStiff.rows();
    buffer.block(row, 0, rows, 1) = mStiff;
    row += rows;
}

void FluidPoint::extractBuffer(CColX &buffer, int &row) {
    int rows = mStiff.rows();
    mStiff += buffer.block(row, 0, rows, 1);
    row += rows;
}

void FluidPoint::scatterDisplToElement(vec_CMatPP &displ, int ipol, int jpol, int maxNu) const {
    // lower orders
    int nyquist = (int)(mNr % 2 == 0);
    for (int alpha = 0; alpha <= mNu - nyquist; alpha++) {
        displ[alpha](ipol, jpol) = mDispl(alpha);
    }
    // mask Nyquist
    if (nyquist) {
        displ[mNu](ipol, jpol) = czero;
    }
    // mask higher orders
    for (int alpha = mNu + 1; alpha <= maxNu; alpha++) {
        displ[alpha](ipol, jpol) = czero;
    }
}

void FluidPoint::gatherStiffFromElement(const vec_CMatPP &stiff, int ipol, int jpol) {
    // lower orders
    int nyquist = (int)(mNr % 2 == 0);
    for (int alpha = 0; alpha <= mNu - nyquist; alpha++) {
        mStiff(alpha) -= stiff[alpha](ipol, jpol);
    }
    // mask Nyquist
    if (nyquist) {
        mStiff(mNu) = czero;
    }
}

void FluidPoint::addToStiff(const CMatX3 &source) {
    // make sure the length of "source" does not exceed mNu + 1
    int rows = std::min(mNu + 1, int(source.rows()));
    mStiff.topRows(rows) -= source.col(0).topRows(rows);
}

void FluidPoint::maskField(CColX &field) {
    field.row(0).imag().setZero();
    // axial boundary condition
    if (mAxial) {
        field.bottomRows(mNu).setZero();
    }
    // mask Nyquist
    if (mNr % 2 == 0) {
        field(mNu) = czero;
    }
}

void FluidPoint::learnWisdom(Real cutoff, int nPeaks) {
    // using only maximum of displacement
    if (nPeaks == 0) {
        // L2 norm
        Real L2norm = mDispl.squaredNorm();
        // Hilbert norm
        Real h2norm = L2norm - .5 * mDispl.row(0).squaredNorm();
        if (h2norm <= mLastMaxDisplWisdom) {
            return;
        }
        mLastMaxDisplWisdom = h2norm;
        // try smaller orders
        Real tol = h2norm * cutoff * cutoff;
        for (int newNu = 0; newNu < mNu; newNu++) {
            Real diff = L2norm - mDispl.topRows(newNu + 1).squaredNorm();
            if (diff <= tol) {
                mNuWisdom = newNu;
                return;
            }
        }
        mNuWisdom = mNu;
    
    // using N peaks of displacement
    } else {
    
        if (mPeaksWisdom >= nPeaks) {
            return;
        }
        // L2 norm
        Real L2norm = mDispl.squaredNorm();
        // Hilbert norm
        Real h2norm = L2norm - .5 * mDispl.row(0).squaredNorm();
    
        if (h2norm <= mLastMaxDisplWisdom) {
            if (mApproachingPeak) {
                // try smaller orders
                Real tol = mLastMaxDisplWisdom * cutoff * cutoff;
                bool found = false;
                for (int newNu = 0; newNu < mNu; newNu++) {
                    Real diff = mLastDisplWisdom.squaredNorm() - mLastDisplWisdom.topRows(newNu + 1).squaredNorm();
                    if (diff <= tol) {
                        mNuWisdom = std::max({newNu, mNuWisdom});
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    mNuWisdom = mNu;
                }
                // reset peak-search 
                mApproachingPeak = false;    
                mPeaksWisdom += 1;
            }
        } else {
            mApproachingPeak = h2norm > tinyDouble;
        }
        mLastDisplWisdom = mDispl;
        mLastMaxDisplWisdom = h2norm;
        
    }
}


