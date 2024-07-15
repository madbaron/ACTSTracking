#include "ACTSSeededCKFTrackingProc.hxx"

#include <EVENT/MCParticle.h>
#include <EVENT/SimTrackerHit.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/LCRelationNavigator.h>
#include <UTIL/LCTrackerConf.h>

#include <Acts/EventData/MultiTrajectory.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
// #include <Acts/Seeding/BinFinder.hpp>
// #include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <Acts/Seeding/SeedFinder.hpp>
#include <Acts/Seeding/SpacePointGrid.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

using namespace Acts::UnitLiterals;

#include "Helpers.hxx"
#include "MeasurementCalibrator.hxx"
#include "SeedSpacePoint.hxx"
#include "SourceLink.hxx"

// Track fitting definitions
using TrackFinderOptions =
    Acts::CombinatorialKalmanFilterOptions<ACTSTracking::SourceLinkAccessor::Iterator,
                                        Acts::VectorMultiTrajectory>;
using SSPoint = ACTSTracking::SeedSpacePoint;

ACTSSeededCKFTrackingProc aACTSSeededCKFTrackingProc;

ACTSSeededCKFTrackingProc::ACTSSeededCKFTrackingProc()
    : ACTSProcBase("ACTSSeededCKFTrackingProc") {
  // modify processor description
  _description =
      "Fit tracks using the Combinatorial Kalman Filter algorithm with triplet "
      "seeding.";

  // Settings
  registerProcessorParameter(
      "RunCKF",
      "Run tracking using CKF. False means stop at the seeding stage.", _runCKF,
      _runCKF);

  registerProcessorParameter(
      "InitialTrackError_RelP",
      "Track error estimate, momentum component (relative).",
      _initialTrackError_relP, 0.25);

  registerProcessorParameter("InitialTrackError_Phi",
                             "Track error estimate, phi (radians).",
                             _initialTrackError_phi, 1_degree);

  registerProcessorParameter("InitialTrackError_Lambda",
                             "Track error estimate, lambda (radians).",
                             _initialTrackError_lambda, 1_degree);

  registerProcessorParameter("InitialTrackError_Pos",
                             "Track error estimate, local position (mm).",
                             _initialTrackError_pos, 10_um);

  // Extrapolation to calo surface
  registerProcessorParameter("CaloFace_Radius",
                             "ECAL Inner Radius (mm).",
                             _caloFaceR, _caloFaceR);
  
  registerProcessorParameter("CaloFace_Z",
                             "ECAL half length (mm).",
                             _caloFaceZ, _caloFaceZ);

  // Seeding configurations
  registerProcessorParameter(
      "SeedingLayers",
      "Layers to use for seeding in format \"VolID LayID\", one per line. ID's "
      "are ACTS GeometryID's. * can be used to wildcard.",
      _seedingLayers, {"*", "*"});

  registerProcessorParameter("SeedFinding_RMax",
                             "Maximum radius of hits to consider.",
                             _seedFinding_rMax, _seedFinding_rMax);

  registerProcessorParameter("SeedFinding_DeltaRMin",
                             "Minimum dR between hits in a seed.",
                             _seedFinding_deltaRMin, _seedFinding_deltaRMin);

  registerProcessorParameter("SeedFinding_DeltaRMax",
                             "Maximum dR between hits in a seed.",
                             _seedFinding_deltaRMax, _seedFinding_deltaRMax);

  registerProcessorParameter("SeedFinding_DeltaRMinTop",
                             "Minimum dR between the reference hit and outer ones in a seed.",
                             _seedFinding_deltaRMinTop, 0.f);

  registerProcessorParameter("SeedFinding_DeltaRMaxTop",
                             "Maximum dR between the reference hit and outer ones in a seed.",
                             _seedFinding_deltaRMaxTop, 0.f);

  registerProcessorParameter("SeedFinding_DeltaRMinBottom",
                             "Minimum dR between the reference hit and inner ones in a seed.",
                             _seedFinding_deltaRMinBottom, 0.f);

  registerProcessorParameter("SeedFinding_DeltaRMaxBottom",
                             "Maximum dR between the reference hit and inner ones in a seed.",
                             _seedFinding_deltaRMaxBottom, 0.f);

  registerProcessorParameter("SeedFinding_zTopBinLen",
                             "Number of top bins along Z for seeding",
                             _zTopBinLen, 1);

  registerProcessorParameter("SeedFinding_zBottomBinLen",
                             "Number of bottom bins along Z for seeding",
                             _zBottomBinLen, 1);

  registerProcessorParameter("SeedFinding_phiTopBinLen",
                             "Number of top bins along phi for seeding",
                             _phiTopBinLen, 1);

  registerProcessorParameter("SeedFinding_phiBottomBinLen",
                             "Number of bottom bins along phi for seeding",
                             _phiBottomBinLen, 1);

  registerProcessorParameter("SeedFinding_zBinEdges",
                             "Bins placement along Z for seeding.",
                             _seedFinding_zBinEdges, StringVec(0));

  registerProcessorParameter(
      "SeedFinding_CollisionRegion",
      "Size of the collision region in one direction (assumed symmetric).",
      _seedFinding_collisionRegion, _seedFinding_collisionRegion);

  registerProcessorParameter("SeedFinding_ZMax",
                             "Maximum z of hits hits to consider.",
                             _seedFinding_zMax, _seedFinding_zMax);

  registerProcessorParameter(
      "SeedFinding_RadLengthPerSeed", "Average radiation length per seed.",
      _seedFinding_radLengthPerSeed, _seedFinding_radLengthPerSeed);

  registerProcessorParameter("SeedFinding_SigmaScattering",
                             "Number of sigmas to allow in scattering angle.",
                             _seedFinding_sigmaScattering,
                             _seedFinding_sigmaScattering);

  registerProcessorParameter("SeedFinding_MinPt",
                             "Minimum pT of tracks to seed.",
                             _seedFinding_minPt, _seedFinding_minPt);

  registerProcessorParameter("SeedFinding_ImpactMax",
                             "Maximum d0 of tracks to seed.",
                             _seedFinding_impactMax, _seedFinding_impactMax);

  registerProcessorParameter("PropagateBackward",
                             "Extrapolates tracks towards beamline.",
                             _propagateBackward, _propagateBackward);

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

  // Output collections - tracks and relations
  registerOutputCollection(LCIO::TRACK, "SeedCollectionName",
                           "Name of seed output collection.",
                           _outputSeedCollection, std::string("SeedTracks"));

  registerOutputCollection(LCIO::TRACK, "TrackCollectionName",
                           "Name of track output collection.",
                           _outputTrackCollection, std::string("Tracks"));
}

void ACTSSeededCKFTrackingProc::init() {
  ACTSProcBase::init();

  // Reset counters
  _fitFails = 0;

  // Initialize seeding layers
  std::vector<std::string> seedingLayers;
  std::copy_if(_seedingLayers.begin(), _seedingLayers.end(),
               std::back_inserter(seedingLayers),
               [](const std::string &s) { return !s.empty(); });

  if (seedingLayers.size() % 2 != 0) {
    throw std::runtime_error("SeedingLayers needs an even number of entries");
  }

  std::vector<Acts::GeometryIdentifier> geoSelection;
  for (uint32_t i = 0; i < seedingLayers.size(); i += 2) {
    Acts::GeometryIdentifier geoid;
    if (_seedingLayers[i + 0] != "*")  // volume
      geoid = geoid.setVolume(std::stoi(_seedingLayers[i + 0]));
    if (_seedingLayers[i + 1] != "*")  // layer
      geoid = geoid.setLayer(std::stoi(_seedingLayers[i + 1]));

    geoSelection.push_back(geoid);
  }

  _seedGeometrySelection = ACTSTracking::GeometryIdSelector(geoSelection);

  if (_seedFinding_deltaRMinTop == 0.f) _seedFinding_deltaRMinTop = _seedFinding_deltaRMin;
  if (_seedFinding_deltaRMaxTop == 0.f) _seedFinding_deltaRMaxTop = _seedFinding_deltaRMax;
  if (_seedFinding_deltaRMinBottom == 0.f) _seedFinding_deltaRMinBottom = _seedFinding_deltaRMin;
  if (_seedFinding_deltaRMaxBottom == 0.f) _seedFinding_deltaRMaxBottom = _seedFinding_deltaRMax;
}

void ACTSSeededCKFTrackingProc::processRunHeader(LCRunHeader *) {}

void ACTSSeededCKFTrackingProc::processEvent(LCEvent *evt) {
  //
  // Prepare the output
  // Make the output track collection
  LCCollectionVec *seedCollection = new LCCollectionVec(LCIO::TRACK);
  LCCollectionVec *trackCollection = new LCCollectionVec(LCIO::TRACK);

  // Enable the track collection to point back to hits
  LCFlagImpl trkFlag(0);
  trkFlag.setBit(LCIO::TRBIT_HITS);
  trackCollection->setFlag(trkFlag.getFlag());

  //
  // Prepare input hits in ACTS format

  // Loop over each hit collections and get a single vector with hits
  // from all of the subdetectors. Also include the Acts GeoId in
  // the vector. It will be important for the sort to speed up the
  // population of the final SourceLink multiset.
  std::vector<std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit *>>
      sortedHits;
  for (const std::string &collection : _inputTrackerHitCollections) {
    // Get the collection of tracker hits
    LCCollection *trackerHitCollection = getCollection(collection, evt);
    if (trackerHitCollection == nullptr) continue;

    for (uint32_t idxHit = 0;
         idxHit < trackerHitCollection->getNumberOfElements(); idxHit++) {
      EVENT::TrackerHit *hit = static_cast<EVENT::TrackerHit *>(
          trackerHitCollection->getElementAt(idxHit));

      sortedHits.push_back(
          std::make_pair(geoIDMappingTool()->getGeometryID(hit), hit));
    }
  }

  // Sort by GeoID
  std::sort(
      sortedHits.begin(), sortedHits.end(),
      [](const std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit *> &hit0,
         const std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit *> &hit1)
          -> bool { return hit0.first < hit1.first; });

  // Turn the LCIO TrackerHit's into Acts objects
  // Assuems that the hits are ssorted by the GeoID
  ACTSTracking::SourceLinkContainer sourceLinks;
  ACTSTracking::MeasurementContainer measurements;
  ACTSTracking::SeedSpacePointContainer spacePoints;

  sourceLinks.reserve(sortedHits.size());
  for (std::pair<Acts::GeometryIdentifier, EVENT::TrackerHit *> &hit :
       sortedHits) {
    // Convert to Acts hit
    const Acts::Surface *surface = trackingGeometry()->findSurface(hit.first);
    if (surface == nullptr) throw std::runtime_error("Surface not found");

    const double *lcioglobalpos = hit.second->getPosition();
    Acts::Vector3 globalPos = {lcioglobalpos[0], lcioglobalpos[1],
                               lcioglobalpos[2]};
    Acts::Result<Acts::Vector2> lpResult =
        surface->globalToLocal(geometryContext(), globalPos, {0, 0, 0}, 0.5_um);
    if (!lpResult.ok())
      throw std::runtime_error(
          "Global to local transformation did not succeed.");

    Acts::Vector2 loc = lpResult.value();

    Acts::SquareMatrix2 localCov = Acts::SquareMatrix2::Zero();
    const EVENT::TrackerHitPlane *hitplane =
        dynamic_cast<const EVENT::TrackerHitPlane *>(hit.second);
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

    //
    // Seed selection and conversion to useful coordinates
    if (_seedGeometrySelection.check(surface->geometryId())) {
      Acts::RotationMatrix3 rotLocalToGlobal =
          surface->referenceFrame(geometryContext(), globalPos, {0, 0, 0});

      // Convert to a seed space point
      // the space point requires only the variance of the transverse and
      // longitudinal position. reduce computations by transforming the
      // covariance directly from local to rho/z.
      //
      // compute Jacobian from global coordinates to rho/z
      //
      //         rho = sqrt(x² + y²)
      // drho/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
      //             = 2 * {x,y} / r
      //       dz/dz = 1 (duuh!)
      //
      double x = globalPos[Acts::ePos0];
      double y = globalPos[Acts::ePos1];
      double scale = 2 / std::hypot(x, y);
      Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
      jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
      jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
      jacXyzToRhoZ(1, Acts::ePos2) = 1;
      // compute Jacobian from local coordinates to rho/z
      Acts::ActsMatrix<2, 2> jac =
          jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(Acts::ePos0, Acts::ePos0);
      // compute rho/z variance
      Acts::ActsVector<2> var = (jac * localCov * jac.transpose()).diagonal();

      // Save spacepoint
      spacePoints.push_back(
          ACTSTracking::SeedSpacePoint(globalPos, var[0], var[1], sourceLink));
    }
  }

  streamlog_out(DEBUG0) << "Created " << spacePoints.size() << " space points"
                        << std::endl;

  //
  // Run seeding + tracking algorithms
  //

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
  if (_propagateBackward) {
    pOptions.direction = Acts::Direction::Backward;
  }

  // Construct a perigee surface as the target surface
  std::shared_ptr<Acts::PerigeeSurface> perigeeSurface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(
          Acts::Vector3{0., 0., 0.});

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
  // Finder configuration
  static const Acts::Vector3 zeropos(0, 0, 0);

  Acts::SeedFinderConfig<ACTSTracking::SeedSpacePoint> finderCfg;
  finderCfg.rMax = _seedFinding_rMax;
  finderCfg.deltaRMin = _seedFinding_deltaRMin;
  finderCfg.deltaRMax = _seedFinding_deltaRMax;
  finderCfg.deltaRMinTopSP = _seedFinding_deltaRMinTop;
  finderCfg.deltaRMaxTopSP = _seedFinding_deltaRMaxTop;
  finderCfg.deltaRMinBottomSP = _seedFinding_deltaRMinBottom;
  finderCfg.deltaRMaxBottomSP = _seedFinding_deltaRMaxBottom;
  finderCfg.collisionRegionMin = -_seedFinding_collisionRegion;
  finderCfg.collisionRegionMax = _seedFinding_collisionRegion;
  finderCfg.zMin = -_seedFinding_zMax;
  finderCfg.zMax = _seedFinding_zMax;
  finderCfg.maxSeedsPerSpM = 1;
  finderCfg.cotThetaMax = 7.40627;  // 2.7 eta;
  finderCfg.sigmaScattering = _seedFinding_sigmaScattering;
  finderCfg.radLengthPerSeed = _seedFinding_radLengthPerSeed;
  finderCfg.minPt = _seedFinding_minPt * Acts::UnitConstants::MeV;
  finderCfg.impactMax = _seedFinding_impactMax * Acts::UnitConstants::mm;
  finderCfg.useVariableMiddleSPRange = true;

  Acts::SeedFilterConfig filterCfg;
  filterCfg.maxSeedsPerSpM = finderCfg.maxSeedsPerSpM;

  finderCfg.seedFilter =
      std::make_unique<Acts::SeedFilter<ACTSTracking::SeedSpacePoint>>(
          Acts::SeedFilter<ACTSTracking::SeedSpacePoint>(filterCfg.toInternalUnits()));
  finderCfg = finderCfg.toInternalUnits().calculateDerivedQuantities();

  Acts::SeedFinderOptions finderOpts;
  finderOpts.bFieldInZ = (*magneticField()->getField(zeropos, magCache))[2];   // TODO investigate
  finderOpts.beamPos = {0, 0};
  finderOpts = finderOpts.toInternalUnits();
  finderOpts = finderOpts.calculateDerivedQuantities(finderCfg);

  Acts::CylindricalSpacePointGridConfig gridCfg;
  gridCfg.cotThetaMax = finderCfg.cotThetaMax;
  gridCfg.deltaRMax = finderCfg.deltaRMax;
  gridCfg.minPt = finderCfg.minPt;
  gridCfg.rMax = finderCfg.rMax;
  gridCfg.zMax = finderCfg.zMax;
  gridCfg.zMin = finderCfg.zMin;
  gridCfg.impactMax = finderCfg.impactMax;
  if (_seedFinding_zBinEdges.size() > 0)
  {
    gridCfg.zBinEdges.resize(_seedFinding_zBinEdges.size());
    for (int k = 0; k < _seedFinding_zBinEdges.size(); k++)
    {
      float pos = std::atof(_seedFinding_zBinEdges[k].c_str());
      if (pos >= finderCfg.zMin && pos < finderCfg.zMax)
      {
        gridCfg.zBinEdges[k] = pos;
      }
      else
      {
        streamlog_out(WARNING) << "Wrong parameter SeedFinding_zBinEdges; "
                             << "default used" << std::endl;
        gridCfg.zBinEdges.clear();
        break;
      }
    }
  }

  Acts::CylindricalSpacePointGridOptions gridOpts;
  gridOpts.bFieldInZ = (*magneticField()->getField(zeropos, magCache))[2];

  // Create tools
  auto extractGlobalQuantities = [](const SSPoint& sp, float, float,
                                     float) {
    Acts::Vector3 position { sp.x(), sp.y(), sp.z() };
    Acts::Vector2 covariance { sp.varianceR(), sp.varianceZ() };
    return std::make_tuple(position, covariance, sp.t());
  };

  std::vector<const ACTSTracking::SeedSpacePoint *> spacePointPtrs(
      spacePoints.size(), nullptr);
  std::transform(spacePoints.begin(), spacePoints.end(), spacePointPtrs.begin(),
                 [](const ACTSTracking::SeedSpacePoint &sp) { return &sp; });

  Acts::Extent rRangeSPExtent;

  Acts::CylindricalSpacePointGrid<SSPoint> grid =
      Acts::CylindricalSpacePointGridCreator::createGrid<SSPoint>(
          gridCfg.toInternalUnits(), gridOpts.toInternalUnits());
  Acts::CylindricalSpacePointGridCreator::fillGrid(finderCfg, finderOpts, grid,
      spacePointPtrs.begin(), spacePointPtrs.end(), extractGlobalQuantities,
      rRangeSPExtent);

  const Acts::GridBinFinder<2ul> bottomBinFinder(_phiBottomBinLen, _zBottomBinLen);
  const Acts::GridBinFinder<2ul> topBinFinder(_phiTopBinLen, _zTopBinLen);

  auto spacePointsGrouping = Acts::CylindricalBinnedGroup<SSPoint>(
      std::move(grid), bottomBinFinder, topBinFinder);

  Acts::SeedFinder<SSPoint> finder(finderCfg);
  decltype(finder)::SeedingState state;
  std::vector<Acts::Seed<SSPoint>> seeds;

  state.spacePointData.resize(spacePointPtrs.size(),
      finderCfg.useDetailedDoubleMeasurementInfo);

  float up = Acts::clampValue<float>(
      std::floor(rRangeSPExtent.max(Acts::binR) / 2) * 2);
  const Acts::Range1D<float> rMiddleSPRange(
      std::floor(rRangeSPExtent.min(Acts::binR) / 2) * 2 +
          finderCfg.deltaRMiddleMinSPRange,
      up - finderCfg.deltaRMiddleMaxSPRange);                  // TODO investigate

  std::vector<Acts::BoundTrackParameters> paramseeds;

  for (const auto [bottom, middle, top] : spacePointsGrouping)
  {
    seeds.clear();

    finder.createSeedsForGroup(
        finderOpts, state, spacePointsGrouping.grid(),
        std::back_inserter(seeds), bottom, middle, top, rMiddleSPRange);

    //
    // Loop over seeds and get track parameters
    paramseeds.clear();
    for (const Acts::Seed<SSPoint> &seed : seeds)
    {
      const SSPoint* bottomSP = seed.sp().front();

      const auto& sourceLink = bottomSP->sourceLink();
      const Acts::GeometryIdentifier& geoId = sourceLink.geometryId();
      const Acts::Surface* surface = trackingGeometry()->findSurface(geoId);
      if (surface == nullptr) {
        std::cout << "surface with geoID " << geoId
                  << " is not found in the tracking gemetry";
        continue;
      }

      // Get the magnetic field at the bottom space point
      const Acts::Vector3 seedPos(bottomSP->x(), bottomSP->y(), bottomSP->z());
      Acts::Result<Acts::Vector3> seedField =
          magneticField()->getField(seedPos, magCache);
      if (!seedField.ok()) {
        throw std::runtime_error("Field lookup error: " +
                                 seedField.error().value());
      }

      std::optional<Acts::BoundVector> optParams =
          Acts::estimateTrackParamsFromSeed(geometryContext(),
                                            seed.sp().begin(), seed.sp().end(),
                                            *surface, *seedField, 0.1_T);
      if (!optParams.has_value()) {
        std::cout << "Failed estimation of track parameters for seed."
                  << std::endl;
        continue;
      }

      const Acts::BoundVector &params = optParams.value();

      float charge = std::copysign(1, params[Acts::eBoundQOverP]);
      float p = std::abs(1 / params[Acts::eBoundQOverP]);

      // build the track covariance matrix using the smearing sigmas
      Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();
      cov(Acts::eBoundLoc0, Acts::eBoundLoc0) =
          std::pow(_initialTrackError_pos, 2);
      cov(Acts::eBoundLoc1, Acts::eBoundLoc1) =
          std::pow(_initialTrackError_pos, 2);
      cov(Acts::eBoundTime, Acts::eBoundTime) =
          std::pow(_initialTrackError_time, 2);
      cov(Acts::eBoundPhi, Acts::eBoundPhi) =
          std::pow(_initialTrackError_phi, 2);
      cov(Acts::eBoundTheta, Acts::eBoundTheta) =
          std::pow(_initialTrackError_lambda, 2);
      cov(Acts::eBoundQOverP, Acts::eBoundQOverP) =
          std::pow(_initialTrackError_relP * p / (p * p), 2);

      Acts::BoundTrackParameters paramseed(surface->getSharedPtr(), params,
                                           cov, Acts::ParticleHypothesis::pion());
      paramseeds.push_back(paramseed);

      //
      // Add seed to LCIO collection
      IMPL::TrackImpl *seedtrack = new IMPL::TrackImpl;
      seedCollection->addElement(seedtrack);

      Acts::Vector3 globalPos = surface->localToGlobal(
          geometryContext(),
          {params[Acts::eBoundLoc0], params[Acts::eBoundLoc1]}, {0, 0, 0});

      // state
      Acts::Result<Acts::Vector3> hitField =
          magneticField()->getField(globalPos, magCache);
      if (!hitField.ok()) {
        throw std::runtime_error("Field lookup error: " +
                                 hitField.error().value());
      }

      EVENT::TrackState *seedTrackState = ACTSTracking::ACTS2Marlin_trackState(
          lcio::TrackState::AtFirstHit, paramseed,
          (*hitField)[2] / Acts::UnitConstants::T);

      // hits
      for (const ACTSTracking::SeedSpacePoint *sp : seed.sp()) {
        const ACTSTracking::SourceLink &sourceLink = sp->sourceLink();
        seedtrack->addHit(sourceLink.lciohit());
      }

      seedtrack->trackStates().push_back(seedTrackState);

      streamlog_out(DEBUG) << "Seed Paramemeters" << std::endl
                           << paramseed << std::endl;
    }

    streamlog_out(DEBUG) << "Seeds found: " << std::endl
                         << paramseeds.size() << std::endl;

    //
    // Find the tracks
    if (!_runCKF) continue;

    using TrackContainer = Acts::TrackContainer<Acts::VectorTrackContainer,
                                                Acts::VectorMultiTrajectory,
                                                std::shared_ptr>;
    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tracks(trackContainer, trackStateContainer);

    for (std::size_t iseed = 0; iseed < paramseeds.size(); ++iseed) {

      tracks.clear();

      auto result = trackFinder.findTracks(paramseeds.at(iseed), ckfOptions, tracks);
      if (result.ok()) {
        const auto& fitOutput = result.value();
        for (const TrackContainer::TrackProxy& trackTip : fitOutput)
        {
          //
          // Helpful debug output
          streamlog_out(DEBUG) << "Trajectory Summary" << std::endl;
          streamlog_out(DEBUG)
              << "\tchi2Sum       " << trackTip.chi2() << std::endl;
          streamlog_out(DEBUG)
              << "\tNDF           " << trackTip.nDoF() << std::endl;
          streamlog_out(DEBUG)
              << "\tnHoles        " << trackTip.nHoles() << std::endl;
          streamlog_out(DEBUG)
              << "\tnMeasurements " << trackTip.nMeasurements() << std::endl;
          streamlog_out(DEBUG)
              << "\tnOutliers     " << trackTip.nOutliers() << std::endl;
          streamlog_out(DEBUG)
              << "\tnStates       " << trackTip.nTrackStates() << std::endl;

          // Make track object
          EVENT::Track *track = ACTSTracking::ACTS2Marlin_track(
              trackTip, magneticField(), magCache, _caloFaceR, _caloFaceZ, geometryContext(), magneticFieldContext(), trackingGeometry());

          // Save results
          trackCollection->addElement(track);
        }
      } else {
        streamlog_out(WARNING)
            << "Track fit error: " << result.error() << std::endl;
        _fitFails++;
      }
    }
  }

  // Save the output seed collection
  evt->addCollection(seedCollection, _outputSeedCollection);

  // Save the output track collection
  evt->addCollection(trackCollection, _outputTrackCollection);
}

void ACTSSeededCKFTrackingProc::check(LCEvent *) {
  // nothing to check here - could be used to fill checkplots in reconstruction
  // processor
}

void ACTSSeededCKFTrackingProc::end() {}

LCCollection *ACTSSeededCKFTrackingProc::getCollection(
    const std::string &collectionName, LCEvent *evt) {
  try {
    return evt->getCollection(collectionName);
  } catch (DataNotAvailableException &e) {
    streamlog_out(DEBUG5) << "- cannot get collection. Collection "
                          << collectionName << " is unavailable" << std::endl;
    return nullptr;
  }
}
