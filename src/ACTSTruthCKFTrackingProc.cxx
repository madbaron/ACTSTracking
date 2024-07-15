#include "ACTSTruthCKFTrackingProc.hxx"

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTrackerConf.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"

using namespace Acts::UnitLiterals;

#include "Helpers.hxx"
#include "MeasurementCalibrator.hxx"
#include "SourceLink.hxx"

using TrackFinderOptions =
    Acts::CombinatorialKalmanFilterOptions<ACTSTracking::SourceLinkAccessor::Iterator,
                                           Acts::VectorMultiTrajectory>;
/*
using TrackFinderResult = Acts::Result<
    Acts::CombinatorialKalmanFilterResult<ACTSTracking::SourceLink>>;

using TrackFinderResultContainer = std::vector<TrackFinderResult>;
*/
ACTSTruthCKFTrackingProc aACTSTruthCKFTrackingProc;

ACTSTruthCKFTrackingProc::ACTSTruthCKFTrackingProc()
    : ACTSProcBase("ACTSTruthCKFTrackingProc") {
  // modify processor description
  _description =
      "Fit tracks using the Combinatorial Kalman Filter algorithm with MC "
      "particles as seeds.";

  // CKF configuration
  registerProcessorParameter("CKF_Chi2CutOff",
                             "Maximum local chi2 contribution.",
                             _CKF_chi2CutOff, _CKF_chi2CutOff);

  registerProcessorParameter(
      "CKF_NumMeasurementsCutOff",
      "Maximum number of associated measurements on a single surface.",
      _CKF_numMeasurementsCutOff, _CKF_numMeasurementsCutOff);

  // Input collections - mc particles, tracker hits and the relationships
  // between them
  registerInputCollections(LCIO::TRACKERHITPLANE, "TrackerHitCollectionNames",
                           "Name of the TrackerHit input collections.",
                           _inputTrackerHitCollections, {});

  registerInputCollection(
      LCIO::MCPARTICLE, "MCParticleCollectionName",
      "Name of the MCParticle input collection (used for seeding).",
      _inputParticleCollection, std::string("MCParticle"));

  // Output collections - tracks and relations
  registerOutputCollection(LCIO::TRACK, "TrackCollectionName",
                           "Name of track output collection.",
                           _outputTrackCollection, std::string("TruthTracks"));

  // Extrapolation to calo surface
  registerProcessorParameter("CaloFace_Radius",
                             "ECAL Inner Radius (mm).",
                             _caloFaceR, _caloFaceR);
  
  registerProcessorParameter("CaloFace_Z",
                             "ECAL half length (mm).",
                             _caloFaceZ, _caloFaceZ);
}

void ACTSTruthCKFTrackingProc::init() {
  ACTSProcBase::init();

  // Reset counters
  _runNumber = 0;
  _eventNumber = 0;
  _fitFails = 0;
}

void ACTSTruthCKFTrackingProc::processRunHeader(LCRunHeader*) { _runNumber++; }

void ACTSTruthCKFTrackingProc::processEvent(LCEvent* evt) {
  // Construct a perigee surface as the target surface
  std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(
          Acts::Vector3{0., 0., 0.});

  //
  // Make seeds using truth particles
  std::vector<Acts::BoundTrackParameters> seeds;

  // Get the collection of MC particles
  LCCollection* particleCollection =
      getCollection(_inputParticleCollection, evt);
  if (particleCollection == nullptr) return;

  for (uint32_t idxP = 0; idxP < particleCollection->getNumberOfElements();
       ++idxP) {
    const MCParticle* mcParticle =
        static_cast<const MCParticle*>(particleCollection->getElementAt(idxP));

    // Tracks are made by stable charged particles from generation
    if (mcParticle->isCreatedInSimulation() ||
        mcParticle->getGeneratorStatus() != 1 || mcParticle->getCharge() == 0)
      continue;

    // Create initial parameters
    double px = mcParticle->getMomentum()[0];
    double py = mcParticle->getMomentum()[1];
    double pz = mcParticle->getMomentum()[2];
    double pt = sqrt(px * px + py * py);
    double p = sqrt(px * px + py * py + pz * pz);

    Acts::BoundVector params = Acts::BoundVector::Zero();
    // position/time
    params[Acts::eBoundLoc0] = 0;
    params[Acts::eBoundLoc1] = 0;
    params[Acts::eBoundTime] = mcParticle->getTime();
    // direction angles phi,theta
    params[Acts::eBoundPhi] = atan2(py, px);
    params[Acts::eBoundTheta] = atan2(pt, pz);
    // absolute momentum vector
    params[Acts::eBoundQOverP] = mcParticle->getCharge() / p;

    // build the track covariance matrix using the smearing sigmas
    Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();
    cov(Acts::eBoundLoc0, Acts::eBoundLoc0) =
        std::pow(_initialTrackError_d0, 2);
    cov(Acts::eBoundLoc1, Acts::eBoundLoc1) =
        std::pow(_initialTrackError_z0, 2);
    cov(Acts::eBoundTime, Acts::eBoundTime) = std::pow(1_ns, 2);
    cov(Acts::eBoundPhi, Acts::eBoundPhi) = std::pow(_initialTrackError_phi, 2);
    cov(Acts::eBoundTheta, Acts::eBoundTheta) =
        std::pow(_initialTrackError_lambda, 2);
    cov(Acts::eBoundQOverP, Acts::eBoundQOverP) =
        std::pow(_initialTrackError_relP * p / (p * p), 2);

    Acts::BoundTrackParameters seed(perigeeSurface, params, cov,
                                    ACTSTracking::convertParticle(mcParticle));
    seeds.push_back(seed);
  }

  //
  // Prepare input hits in ACTS format

  // Loop over each hit collections and get a single vector with hits
  // from all of the subdetectors. Also include the Acts GeoId in
  // the vector. It will be important for the sort to speed up the
  // population of the final SourceLink multiset.
  std::vector<std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit*>>
      sortedHits;
  for (const std::string& collection : _inputTrackerHitCollections) {
    // Get the collection of tracker hits
    LCCollection* trackerHitCollection = getCollection(collection, evt);
    if (trackerHitCollection == nullptr) continue;

    for (uint32_t idxHit = 0;
         idxHit < trackerHitCollection->getNumberOfElements(); idxHit++) {
      EVENT::TrackerHit* hit = static_cast<EVENT::TrackerHit*>(
          trackerHitCollection->getElementAt(idxHit));

      sortedHits.push_back(
          std::make_pair(geoIDMappingTool()->getGeometryID(hit), hit));
    }
  }

  // Sort by GeoID
  std::sort(
      sortedHits.begin(), sortedHits.end(),
      [](const std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit*>& hit0,
         const std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit*>& hit1)
          -> bool { return hit0.first < hit1.first; });

  // Turn the LCIO TrackerHit's into Acts objects
  // Assuems that the hits are ssorted by the GeoID
  ACTSTracking::SourceLinkContainer sourceLinks;
  ACTSTracking::MeasurementContainer measurements;

  sourceLinks.reserve(sortedHits.size());
  for (std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit*>& hit :
       sortedHits) {
    // Convert to Acts hit
    const Acts::Surface* surface = trackingGeometry()->findSurface(hit.first);

    const double* lcioglobalpos = hit.second->getPosition();
    Acts::Vector3 globalPos = {lcioglobalpos[0], lcioglobalpos[1],
                               lcioglobalpos[2]};
    Acts::Result<Acts::Vector2> lpResult =
        surface->globalToLocal(geometryContext(), globalPos, {0, 0, 0}, 0.5_um);
    if (!lpResult.ok())
      throw std::runtime_error(
          "Global to local transformation did not succeed.");

    Acts::Vector2 loc = lpResult.value();

    Acts::SquareMatrix2 localCov = Acts::SquareMatrix2::Zero();
    const EVENT::TrackerHitPlane* hitplane =
        dynamic_cast<const EVENT::TrackerHitPlane*>(hit.second);
    if (hitplane) {
      localCov(0, 0) = std::pow(hitplane->getdU() * Acts::UnitConstants::mm, 2);
      localCov(1, 1) = std::pow(hitplane->getdV() * Acts::UnitConstants::mm, 2);
    } else {
      throw std::runtime_error("Currently only support TrackerHitPlane.");
    }

    ACTSTracking::SourceLink sourceLink(surface->geometryId(),
                                        measurements.size(), hit.second);
    Acts::SourceLink src_wrap { sourceLink };
    Acts::Measurement meas = Acts::makeMeasurement(
        src_wrap, loc, localCov, Acts::eBoundLoc0, Acts::eBoundLoc1);

    measurements.push_back(meas);
    sourceLinks.emplace_hint(sourceLinks.end(), sourceLink);
  }

  //
  // Caches
  Acts::MagneticFieldContext magFieldContext = Acts::MagneticFieldContext();
  Acts::MagneticFieldProvider::Cache magCache =
      magneticField()->makeCache(magFieldContext);

  //
  // Initialize track finder
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;
  using Stepper = Acts::EigenStepper<>;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using CKF = Acts::CombinatorialKalmanFilter<Propagator, Acts::VectorMultiTrajectory>;

  // Configurations
  Navigator::Config navigatorCfg{trackingGeometry()};
  navigatorCfg.resolvePassive = false;
  navigatorCfg.resolveMaterial = true;
  navigatorCfg.resolveSensitive = true;

  // construct all components for the fitter
  Stepper stepper(magneticField());
  Navigator navigator(navigatorCfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  CKF trackFinder(std::move(propagator));

  // Set the options
  Acts::MeasurementSelector::Config measurementSelectorCfg = {
      {Acts::GeometryIdentifier(),
       { {}, { _CKF_chi2CutOff }, { (std::size_t)(_CKF_numMeasurementsCutOff) }}}};

  Acts::PropagatorPlainOptions pOptions;
  pOptions.maxSteps = 10000;

  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;

  Acts::MeasurementSelector measSel { measurementSelectorCfg };
  ACTSTracking::MeasurementCalibrator measCal { measurements };
  Acts::CombinatorialKalmanFilterExtensions<Acts::VectorMultiTrajectory>
      extensions;
  extensions.calibrator.connect<
      &ACTSTracking::MeasurementCalibrator::calibrate>(
      &measCal);
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &kfUpdater);
  extensions.smoother.connect<
      &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
      &kfSmoother);
  extensions.measurementSelector
      .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
          &measSel);

  using ACTSTracking::SourceLinkAccessor;
  SourceLinkAccessor slAccessor;
  slAccessor.container = &sourceLinks;
  Acts::SourceLinkAccessorDelegate<SourceLinkAccessor::Iterator> slAccessorDelegate;
  slAccessorDelegate.connect<&SourceLinkAccessor::range>(&slAccessor);

  // std::unique_ptr<const Acts::Logger>
  // logger=Acts::getDefaultLogger("TrackFitting",
  // Acts::Logging::Level::VERBOSE);
  TrackFinderOptions ckfOptions = TrackFinderOptions(
      geometryContext(), magneticFieldContext(), calibrationContext(),
      slAccessorDelegate,
      extensions, pOptions, perigeeSurface.get());

  //
  // Output

  // Make the output track collection
  LCCollectionVec* trackCollection = new LCCollectionVec(LCIO::TRACK);

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackCollection->setFlag(trkFlag.getFlag());

  //
  // Find the tracks
  using TrackContainer =
      Acts::TrackContainer<Acts::VectorTrackContainer,
                           Acts::VectorMultiTrajectory, std::shared_ptr>;
  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  for (std::size_t iseed = 0; iseed < seeds.size(); ++iseed) {
    auto result = trackFinder.findTracks(seeds.at(iseed), ckfOptions, tracks);
    if (result.ok()) {
      const auto& fitOutput = result.value();
      for (const auto& trackTip : fitOutput)
      {
        std::cout << trackTip.chi2() << std::endl;
        EVENT::Track* track = ACTSTracking::ACTS2Marlin_track(
            trackTip, magneticField(), magCache,
            _caloFaceR, _caloFaceZ, geometryContext(), magneticFieldContext(), trackingGeometry());
        trackCollection->addElement(track);
      }
    } else {
      streamlog_out(WARNING) << "Track fit error: " << result.error() << std::endl;
      _fitFails++;
    }

  }

  // Save the output track collection
  evt->addCollection(trackCollection, _outputTrackCollection);

  // Increment the event number
  _eventNumber++;
}

void ACTSTruthCKFTrackingProc::check(LCEvent*) {
  // nothing to check here - could be used to fill checkplots in reconstruction
  // processor
}

void ACTSTruthCKFTrackingProc::end() {
  streamlog_out(MESSAGE) << " end()  " << name() << " processed "
                         << _eventNumber << " events in " << _runNumber
                         << " runs " << std::endl;
}

LCCollection* ACTSTruthCKFTrackingProc::getCollection(
    const std::string& collectionName, LCEvent* evt) {
  try {
    return evt->getCollection(collectionName);
  } catch (DataNotAvailableException& e) {
    streamlog_out(DEBUG5) << "- cannot get collection. Collection "
                          << collectionName << " is unavailable" << std::endl;
    return nullptr;
  }
}
