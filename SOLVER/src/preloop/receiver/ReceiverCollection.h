// Receiver.h
// created by Kuangdai on 1-Jun-2016 
// receiver collections

#pragma once

#include <vector>
#include <string>

class Domain;
class Mesh;
class Parameters;
class Receiver;
class PointwiseIO;

class ReceiverCollection {
public:
    ReceiverCollection(const std::string &fileRec, bool geographic, 
        double srcLat, double srcLon, double srcDep, int duplicated, 
        double saveSurfRadius, double saveSurfDistMin, double saveSurfDistMax, 
        bool saveSurfUpper, bool cartesian);
    ~ReceiverCollection();
    
    void release(Domain &domain, const Mesh &mesh, bool depthInRef); 
    
    std::string verbose() const;
    
        
    static void buildInparam(ReceiverCollection *&rec, const Parameters &par,
        double srcLat, double srcLon, double srcDep, int totalStepsSTF,
        bool cartesian, int verbose);
private:
    
    // receivers
    std::vector<Receiver *> mReceivers;
    
    // input
    std::string mInputFile;
    bool mGeographic, mCartesian;

    // options
    int mTotalRecordSteps = 0;
    int mRecordInterval = 1;
    int mBufferSize = 1000;
    std::string mComponents = "RTZ";
    
    // IO
    std::vector<PointwiseIO *> mPointwiseIO;

    // for verbose
    int mWidthName;
    int mWidthNetwork;
    
    // surface wavefield
    double mSaveSurfaceAtRadius = -1.;
    double mSaveSurfaceDistMin = 0.;
    double mSaveSurfaceDistMax = 180.;
    bool mSaveSurfaceFromUpper = false;
    bool mAssemble = true;
    
    // source location
    double mSrcLat, mSrcLon, mSrcDep;
};

