#ifndef ACTSSeededCKFTrackingProc_h
#define ACTSSeededCKFTrackingProc_h 1

#include <EVENT/TrackerHit.h>

#include <UTIL/CellIDDecoder.h>

#include <Acts/Definitions/Units.hpp>

#include "ACTSProcBase.hxx"
#include "GeometryIdSelector.hxx"

/**
 * This code performs a true pattern recognition by looping over all MC
 * particles and adding all hits associated to them onto a prototrack. This is
 * then fitted and output.
 */
class ACTSSeededCKFTrackingProc : public ACTSProcBase {
 public:
  virtual marlin::Processor *newProcessor() {
    return new ACTSSeededCKFTrackingProc;
  }

  ACTSSeededCKFTrackingProc(const ACTSSeededCKFTrackingProc &) = delete;
  ACTSSeededCKFTrackingProc &operator=(const ACTSSeededCKFTrackingProc &) =
      delete;
  ACTSSeededCKFTrackingProc();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

 private:
  /** Call to get collections
   */
  LCCollection *getCollection(const std::string &, LCEvent *);

 protected:
  // Collection names for (in/out)put
  std::vector<std::string> _inputTrackerHitCollections;
  std::string _outputSeedCollection;
  std::string _outputTrackCollection;

  // Run settings
  bool _runCKF = true;
  bool _propagateBackward = false;

  // Extrapolation to calo settings
  float _caloFaceR = 1857; //mm
  float _caloFaceZ = 2307; //mm

  // Seed finding configuration
  float _seedFinding_rMax = 150;
  float _seedFinding_deltaRMin = 5;
  float _seedFinding_deltaRMax = 80;
  float _seedFinding_deltaRMinTop = 0;
  float _seedFinding_deltaRMaxTop = 0;
  float _seedFinding_deltaRMinBottom = 0;
  float _seedFinding_deltaRMaxBottom = 0;
  float _seedFinding_collisionRegion = 75;
  float _seedFinding_zMax = 600;
  float _seedFinding_sigmaScattering = 50;
  float _seedFinding_radLengthPerSeed = 0.1;
  float _seedFinding_minPt = 500;
  float _seedFinding_impactMax = 3 * Acts::UnitConstants::mm;

  StringVec _seedFinding_zBinEdges;
  int _zTopBinLen = 1;
  int _zBottomBinLen = 1;
  int _phiTopBinLen = 1;
  int _phiBottomBinLen = 1;

  // Track fit parameters
  double _initialTrackError_pos;
  double _initialTrackError_phi;
  double _initialTrackError_relP;
  double _initialTrackError_lambda;
  double _initialTrackError_time =
      100 * Acts::UnitConstants::ns;  // No Marlin default

  double _CKF_chi2CutOff = 15;
  int32_t _CKF_numMeasurementsCutOff = 10;

  // Seeding configuration
  std::vector<std::string> _seedingLayers;
  ACTSTracking::GeometryIdSelector _seedGeometrySelection;

  uint32_t _fitFails;
};

#endif
