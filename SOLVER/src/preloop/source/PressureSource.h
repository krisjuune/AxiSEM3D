// PressureSource.h
// created by Claudia on 6-Apr-2018
// axial earthquake source

#pragma once
#include "Source.h"

class PressureSource: public Source {
public:
    // input: CMTSOLUTION components
    PressureSource(double depth = 0., double lat = 0., double lon = 0.,
        double M0 = 0.);

    std::string verbose() const;

protected:
    void computeSourceFourier(const Quad &myQuad, const RDColP &interpFactZ,
        arPP_CMatX3 &fouriers) const;

private:
    // store: Cartesian, paper components
    double mM0;
};
