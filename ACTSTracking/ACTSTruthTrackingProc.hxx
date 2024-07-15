#ifndef ACTSTruthTrackingProc_h
#define ACTSTruthTrackingProc_h 1

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>

#include <UTIL/CellIDDecoder.h>

#include <Acts/Definitions/Units.hpp>
#include "Acts/EventData/ParticleHypothesis.hpp"

#include "ACTSProcBase.hxx"

/**
 * This code performs a true pattern recognition by looping over all MC
 * particles and adding all hits associated to them onto a prototrack. This is
 * then fitted and output.
 */
class ACTSTruthTrackingProc : public ACTSProcBase {
 public:
  virtual marlin::Processor* newProcessor() {
    return new ACTSTruthTrackingProc;
  }

  ACTSTruthTrackingProc(const ACTSTruthTrackingProc&) = delete;
  ACTSTruthTrackingProc& operator=(const ACTSTruthTrackingProc&) = delete;
  ACTSTruthTrackingProc();

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent* evt);

  virtual void check(LCEvent* evt);

  /** Called after data processing for clean up.
   */
  virtual void end();

 private:
  /** Call to get collections
   */
  LCCollection* getCollection(const std::string&, LCEvent*);

 protected:
  // Encoder
  std::shared_ptr<UTIL::BitField64> _encoder;

  // Get the subdetector ID from a hit
  int getSubdetector(const lcio::TrackerHit* hit);

  // Get the layer ID from a hit
  int getLayer(const lcio::TrackerHit* hit);

  // Remove hits in the same layer of the same subdetector
  void removeHitsSameLayer(const std::vector<lcio::TrackerHit*>&,
                           std::vector<lcio::TrackerHit*>&);

  // Collection names for (in/out)put
  std::vector<std::string> _inputTrackerHitCollections;
  std::vector<std::string> _inputTrackerHitRelationCollections;
  std::string _inputParticleCollection;
  std::string _outputTrackCollection;

  // Run and event counters
  uint32_t _eventNumber;
  uint32_t _runNumber;

  // Extrapolation to calo settings
  float _caloFaceR = 1857; //mm
  float _caloFaceZ = 2307; //mm

  // Track fit parameters
  double _initialTrackError_d0 = 20 * Acts::UnitConstants::um;  // Marlin: 1.e3
  double _initialTrackError_phi =
      1 * Acts::UnitConstants::degree;    // Marlin default: 1.e1
  double _initialTrackError_relP = 0.01;  // Marlin default: 1.e-2
  double _initialTrackError_lambda =
      1 * Acts::UnitConstants::degree;  // Marlin (tanlambda) default: 1.e1
  double _initialTrackError_z0 =
      100 * Acts::UnitConstants::um;  // Marlin default: 1.e3
  double _initialTrackError_time =
      1 * Acts::UnitConstants::ns;  // No Marlin default
  // double _maxChi2perHit;

  uint32_t _fitFails;
};

#endif
