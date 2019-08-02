// AutoGeometricParams.h
// created by Claudia on 11-Oct-2018

#include "eigenp.h"

class AutoGeometricParams {

public:
    RDMatXX mUndulation;
    double mRLayer, mRUpper, mRLower;
    RDColX mX, mY;

    AutoGeometricParams(double z_mesh, RDColX x, RDColX y): mX(x), mY(y), mRLayer(z_mesh) {
        mUndulation = RDMatXX::Zero(mX.rows(),mY.rows());
    };

    void initialize() {
        double maxdiff = std::max({mUndulation.maxCoeff() * mUndulation.maxCoeff(), mUndulation.minCoeff() * mUndulation.minCoeff()});
        maxdiff = std::sqrt(maxdiff);
        mRLower = std::min({mRLayer + 3 * maxdiff, mMeshLimits(1)});
        mRUpper = std::max({mRLayer - 3 * maxdiff, mMeshLimits(0)});
    };

    void setMeshLimits(RDCol2 meshLimits) {mMeshLimits = meshLimits;};
    void setDepth(int xi, int yi, double z) {
        if ((xi >= mX.rows()) || (yi >= mY.rows())) {
            throw std::runtime_error("AutoGeometricParams::setDepth || Point out of bounds.");
        }
        mUndulation(xi, yi) = z - mRLayer;
    }

private:
    RDCol2 mMeshLimits;

};
