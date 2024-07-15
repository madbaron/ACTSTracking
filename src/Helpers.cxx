#include "Helpers.hxx"

#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>

#include <Acts/Surfaces/CylinderSurface.hpp>
#include <Acts/Surfaces/DiscSurface.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/ParticleHypothesis.hpp>

#include <filesystem>

#include "config.h"

#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>


using Stepper = Acts::EigenStepper<>;
using Navigator = Acts::Navigator;
using Propagator = Acts::Propagator<Stepper, Navigator>;

namespace ACTSTracking {

std::string findFile(const std::string& inpath) {
  if (inpath.empty()) return inpath;

  // Already absolute path
  if (inpath[0] == '/') return inpath;

  // relative to cwd
  if (std::filesystem::exists(inpath)) {
    return inpath;
  }

  // relative to absolute paths
  if (std::filesystem::exists(ACTSTRACKING_SOURCEDIR + inpath)) {
    return ACTSTRACKING_SOURCEDIR + inpath;
  }

  if (std::filesystem::exists(ACTSTRACKING_DATADIR + inpath)) {
    return ACTSTRACKING_DATADIR + inpath;
  }

  // nothing was found :(
  return inpath;
}
/*
EVENT::Track* ACTS2Marlin_track(
    const Acts::CombinatorialKalmanFilterResult<ACTSTracking::SourceLink>&
        fitOutput,
    std::size_t trackTip,
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
    Acts::MagneticFieldProvider::Cache& magCache) {
  IMPL::TrackImpl* track = new IMPL::TrackImpl;
  Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(fitOutput.fittedStates,
                                                    trackTip);

  //
  // Fit state
  track->setChi2(trajState.chi2Sum);
  track->setNdf(trajState.NDF);

  // Track state at IP
  static const Acts::Vector3 zeroPos(0, 0, 0);
  Acts::Result<Acts::Vector3> fieldRes =
      magneticField->getField(zeroPos, magCache);
  if (!fieldRes.ok()) {
    throw std::runtime_error("Field lookup error: " + fieldRes.error().value());
  }
  Acts::Vector3 field = *fieldRes;

  const Acts::BoundTrackParameters& params =
      fitOutput.fittedParameters.at(trackTip);
  EVENT::TrackState* trackStateAtIP = ACTSTracking::ACTS2Marlin_trackState(
      EVENT::TrackState::AtIP, params, field[2] / Acts::UnitConstants::T);
  track->trackStates().push_back(trackStateAtIP);

  //
  // Hits and associated track states
  EVENT::TrackerHitVec hitsOnTrack;
  EVENT::TrackStateVec statesOnTrack;
  fitOutput.fittedStates.visitBackwards(
      trackTip, [&](const Acts::MultiTrajectory<
                    ACTSTracking::SourceLink>::ConstTrackStateProxy& state) {
        // No measurement at this state
        if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return true;
        }

        // register all particles that generated this hit
        EVENT::TrackerHit* myHit = state.uncalibrated().lciohit();
        hitsOnTrack.push_back(myHit);

        // Save track state information
        const Acts::Vector3 hitPos(myHit->getPosition()[0],
                                   myHit->getPosition()[1],
                                   myHit->getPosition()[2]);
        Acts::Result<Acts::Vector3> fieldRes =
            magneticField->getField(hitPos, magCache);
        if (!fieldRes.ok()) {
          throw std::runtime_error("Field lookup error: " +
                                   fieldRes.error().value());
        }
        Acts::Vector3 field = *fieldRes;

        EVENT::TrackState* trackState = ACTSTracking::ACTS2Marlin_trackState(
            EVENT::TrackState::AtOther, state.smoothed(),
            state.smoothedCovariance(), field[2] / Acts::UnitConstants::T);
        statesOnTrack.push_back(trackState);

        return true;
      });

  // Reverse hits and states, above creates them backwards
  std::reverse(hitsOnTrack.begin(), hitsOnTrack.end());
  std::reverse(statesOnTrack.begin(), statesOnTrack.end());

  // Save hits
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(
      lcio::LCTrackerCellID::encoding_string());
  EVENT::IntVec& subdetectorHitNumbers = track->subdetectorHitNumbers();
  for (EVENT::TrackerHit* hit : hitsOnTrack) {
    track->addHit(hit);

    uint32_t sysid = decoder(hit)["system"];
    if (subdetectorHitNumbers.size() <= sysid) {
      subdetectorHitNumbers.resize(sysid + 1, 0);
    }
    subdetectorHitNumbers[sysid]++;
  }

  // Save the track states at hits
  if (statesOnTrack.size() > 0) {
    dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.back())
        ->setLocation(EVENT::TrackState::AtLastHit);
    dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.front())
        ->setLocation(EVENT::TrackState::AtFirstHit);
  }

  EVENT::TrackStateVec& myTrackStates = track->trackStates();
  myTrackStates.insert(myTrackStates.end(), statesOnTrack.begin(),
                       statesOnTrack.end());

  return track;
}

EVENT::Track* ACTS2Marlin_track(
    const Acts::KalmanFitterResult<ACTSTracking::SourceLink>& fitOutput,
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
    Acts::MagneticFieldProvider::Cache& magCache) {
  IMPL::TrackImpl* track = new IMPL::TrackImpl;
  Acts::MultiTrajectoryHelpers::TrajectoryState trajState =
      Acts::MultiTrajectoryHelpers::trajectoryState(
          fitOutput.fittedStates, fitOutput.lastMeasurementIndex);

  //
  // Fit state
  track->setChi2(trajState.chi2Sum);
  track->setNdf(trajState.NDF);

  // Track state at IP
  static const Acts::Vector3 zeroPos(0, 0, 0);
  Acts::Result<Acts::Vector3> fieldRes =
      magneticField->getField(zeroPos, magCache);
  if (!fieldRes.ok()) {
    throw std::runtime_error("Field lookup error: " + fieldRes.error().value());
  }
  Acts::Vector3 field = *fieldRes;

  const Acts::BoundTrackParameters& params = fitOutput.fittedParameters.value();
  EVENT::TrackState* trackStateAtIP = ACTSTracking::ACTS2Marlin_trackState(
      EVENT::TrackState::AtIP, params, field[2] / Acts::UnitConstants::T);
  track->trackStates().push_back(trackStateAtIP);

  //
  // Hits and associated track states
  EVENT::TrackerHitVec hitsOnTrack;
  EVENT::TrackStateVec statesOnTrack;
  fitOutput.fittedStates.visitBackwards(
      fitOutput.lastMeasurementIndex,
      [&](const Acts::MultiTrajectory<
          ACTSTracking::SourceLink>::ConstTrackStateProxy& state) {
        // No measurement at this state
        if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return true;
        }

        // Save hit information
        EVENT::TrackerHit* myHit = state.uncalibrated().lciohit();
        hitsOnTrack.push_back(myHit);

        // Save track state information
        const Acts::Vector3 hitPos(myHit->getPosition()[0],
                                   myHit->getPosition()[1],
                                   myHit->getPosition()[2]);
        Acts::Result<Acts::Vector3> fieldRes =
            magneticField->getField(hitPos, magCache);
        if (!fieldRes.ok()) {
          throw std::runtime_error("Field lookup error: " +
                                   fieldRes.error().value());
        }
        Acts::Vector3 field = *fieldRes;

        EVENT::TrackState* trackState = ACTSTracking::ACTS2Marlin_trackState(
            EVENT::TrackState::AtOther, state.smoothed(),
            state.smoothedCovariance(), hitPos[2] / Acts::UnitConstants::T);
        statesOnTrack.push_back(trackState);

        return true;
      });

  // Reverse hits and states, above creates them backwards
  std::reverse(hitsOnTrack.begin(), hitsOnTrack.end());
  std::reverse(statesOnTrack.begin(), statesOnTrack.end());

  // Save hits
  for (EVENT::TrackerHit* hit : hitsOnTrack) {
    track->addHit(hit);
  }

  // Save the track states at hits
  if (statesOnTrack.size() > 0) {
    dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.back())
        ->setLocation(EVENT::TrackState::AtLastHit);
    dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.front())
        ->setLocation(EVENT::TrackState::AtFirstHit);
  }

  EVENT::TrackStateVec& myTrackStates = track->trackStates();
  myTrackStates.insert(myTrackStates.end(), statesOnTrack.begin(),
                       statesOnTrack.end());

  return track;
}

*/

// Conversion with extrapolation to calorimeter face
EVENT::Track* ACTS2Marlin_track(
    const TrackResult& fitter_res,
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
    Acts::MagneticFieldProvider::Cache& magCache, double caloFaceR, double caloFaceZ, 
    Acts::GeometryContext geoContext, Acts::MagneticFieldContext magFieldContext, 
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeo)
{
  IMPL::TrackImpl* track = new IMPL::TrackImpl;

  track->setChi2(fitter_res.chi2());
  track->setNdf(fitter_res.nDoF());
  track->setNholes(fitter_res.nHoles());

  const Acts::Vector3 zeroPos(0, 0, 0);
  Acts::Result<Acts::Vector3> fieldRes = magneticField->getField(zeroPos, magCache);
  if (!fieldRes.ok()) {
    throw std::runtime_error("Field lookup error: " + fieldRes.error().value());
  }
  Acts::Vector3 field = *fieldRes;

  const Acts::BoundVector& params = fitter_res.parameters();
  const Acts::BoundMatrix& covariance = fitter_res.covariance();
  EVENT::TrackState* trackStateAtIP = ACTSTracking::ACTS2Marlin_trackState(
      EVENT::TrackState::AtIP, params, covariance, field[2] / Acts::UnitConstants::T);
  track->trackStates().push_back(trackStateAtIP);

  EVENT::TrackerHitVec hitsOnTrack;
  EVENT::TrackStateVec statesOnTrack;

  for (const auto& trk_state : fitter_res.trackStatesReversed())
  {
    if (!trk_state.hasUncalibratedSourceLink()) continue;

    auto sl = trk_state.getUncalibratedSourceLink()
                       .get<ACTSTracking::SourceLink>();
    EVENT::TrackerHit* curr_hit = sl.lciohit();
    hitsOnTrack.push_back(curr_hit);

    const Acts::Vector3 hitPos(curr_hit->getPosition()[0],
                               curr_hit->getPosition()[1],
                               curr_hit->getPosition()[2]);

    Acts::Result<Acts::Vector3> fieldRes =
        magneticField->getField(hitPos, magCache);
    if (!fieldRes.ok()) {
      throw std::runtime_error("Field lookup error: " +
                               fieldRes.error().value());
    }
    Acts::Vector3 field = *fieldRes;

    EVENT::TrackState* trackState = ACTSTracking::ACTS2Marlin_trackState(
            EVENT::TrackState::AtOther, trk_state.smoothed(),
            trk_state.smoothedCovariance(), field[2] / Acts::UnitConstants::T);
    statesOnTrack.push_back(trackState);
  }

  // Create the CaloSurface
  auto caloCylinder = std::make_shared<Acts::CylinderBounds>(caloFaceR, caloFaceZ);
  auto caloSurface = Acts::Surface::makeShared<Acts::CylinderSurface>(Acts::Transform3::Identity(), caloCylinder);
  
  // Define the circle dimensions (circles at both ends of the cylinder)
  Acts::Translation3 circlePosition1(0, 0, -caloFaceZ); // circle at -z end
  Acts::Translation3 circlePosition2(0, 0, caloFaceZ);  // circle at +z end

  // Create the circle surfaces
  auto circleSurface1 = Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3(circlePosition1), 0. ,caloFaceR);
  auto circleSurface2 = Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3(circlePosition2), 0., caloFaceR);

  // define start parameters - swap this out with some smart call
  double d0 = params[Acts::eBoundLoc0];
  double z0 = params[Acts::eBoundLoc1];
  double phi = params[Acts::eBoundPhi];
  double theta = params[Acts::eBoundTheta];
  double qoverp = params[Acts::eBoundQOverP];
  double time = params[Acts::eBoundTime];

  Acts::Vector3 pos(d0 * cos(phi), d0 * sin(phi), z0);
  Acts::CurvilinearTrackParameters start(Acts::VectorHelpers::makeVector4(pos, time), phi, theta,
                                         qoverp, covariance, Acts::ParticleHypothesis::pion());

  // Set propagator options
  Acts::PropagatorOptions myCaloPropOptions(geoContext, magFieldContext);
  myCaloPropOptions.pathLimit = 20 * Acts::UnitConstants::m;

  //Let's try with our private propagator
  // Configurations
  Navigator::Config navigatorCfg{trackingGeo};
  navigatorCfg.resolvePassive = false;
  navigatorCfg.resolveMaterial = true;
  navigatorCfg.resolveSensitive = true;

  // construct all components for the fitter
  Stepper stepper(magneticField);
  Navigator navigator(navigatorCfg);
  Propagator mypropagator(std::move(stepper), std::move(navigator));
  auto resultProp = mypropagator.propagate(start, *caloSurface, myCaloPropOptions);
  if (resultProp.ok()) {
    auto end_parameters = resultProp.value().endParameters;
    const Acts::BoundMatrix& atCaloCovariance = *(end_parameters->covariance());

    EVENT::TrackState* trackStateAtCalo = ACTSTracking::ACTS2Marlin_trackState(
        EVENT::TrackState::AtCalorimeter, end_parameters->parameters(),
        atCaloCovariance, field[2] / Acts::UnitConstants::T);
    track->trackStates().push_back(trackStateAtCalo);
  }
  else {
    if (theta > M_PI/2){
        auto resultPropP = mypropagator.propagate(start, *circleSurface1, myCaloPropOptions);
        if (resultPropP.ok()) {
          auto end_parametersP = resultPropP.value().endParameters;
          const Acts::BoundMatrix& atCaloCovariance = *(end_parametersP->covariance());

          EVENT::TrackState* trackStateAtCalo = ACTSTracking::ACTS2Marlin_trackState(
            EVENT::TrackState::AtCalorimeter, end_parametersP->parameters(),
            atCaloCovariance, field[2] / Acts::UnitConstants::T);
          track->trackStates().push_back(trackStateAtCalo);
        }
        else{
          std::cout << "Failed propagation! " << std::endl;
        }
    }
    else{
        auto resultPropM = mypropagator.propagate(start, *circleSurface2, myCaloPropOptions);
        if (resultPropM.ok()) {
          auto end_parametersM = resultPropM.value().endParameters;
          const Acts::BoundMatrix& atCaloCovariance = *(end_parametersM->covariance());

          EVENT::TrackState* trackStateAtCalo = ACTSTracking::ACTS2Marlin_trackState(
            EVENT::TrackState::AtCalorimeter, end_parametersM->parameters(),
            atCaloCovariance, field[2] / Acts::UnitConstants::T);
          track->trackStates().push_back(trackStateAtCalo);
        }
        else{
          std::cout << "Failed propagation!" << std::endl;
        }      
    }
  }
  
  std::reverse(hitsOnTrack.begin(), hitsOnTrack.end());
  std::reverse(statesOnTrack.begin(), statesOnTrack.end());

  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(
      lcio::LCTrackerCellID::encoding_string());
  EVENT::IntVec& subdetectorHitNumbers = track->subdetectorHitNumbers();
  subdetectorHitNumbers.resize(7, 0);
  for (EVENT::TrackerHit* hit : hitsOnTrack) {
    track->addHit(hit);

    uint32_t sysid = decoder(hit)["system"];
    if (subdetectorHitNumbers.size() <= sysid) {
      subdetectorHitNumbers.resize(sysid + 1, 0);
    }
    subdetectorHitNumbers[sysid]++;
  }

  if (statesOnTrack.size() > 0) {
    dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.back())
        ->setLocation(EVENT::TrackState::AtLastHit);
    dynamic_cast<IMPL::TrackStateImpl*>(statesOnTrack.front())
        ->setLocation(EVENT::TrackState::AtFirstHit);
  }

  EVENT::TrackStateVec& myTrackStates = track->trackStates();
  myTrackStates.insert(myTrackStates.end(), statesOnTrack.begin(),
                       statesOnTrack.end());

  return track;
}



EVENT::TrackState* ACTS2Marlin_trackState(
    int location, const Acts::BoundTrackParameters& params, double Bz) {
  return ACTS2Marlin_trackState(location, params.parameters(),
                                params.covariance().value(), Bz);
}

EVENT::TrackState* ACTS2Marlin_trackState(int location,
                                          const Acts::BoundVector& value,
                                          const Acts::BoundMatrix& cov,
                                          double Bz) {
  // Create new object
  IMPL::TrackStateImpl* trackState = new IMPL::TrackStateImpl();

  // Basic properties
  trackState->setLocation(location);

  //
  // Trajectory parameters

  // Central values
  double d0 = value[Acts::eBoundLoc0];
  double z0 = value[Acts::eBoundLoc1];
  double phi = value[Acts::eBoundPhi];
  double theta = value[Acts::eBoundTheta];
  double qoverp = value[Acts::eBoundQOverP];

  double p = 1e3 / qoverp;
  double omega = (0.3 * Bz) / (p * std::sin(theta));
  double lambda = M_PI / 2 - theta;
  double tanlambda = std::tan(lambda);

  trackState->setPhi(phi);
  trackState->setTanLambda(tanlambda);
  trackState->setOmega(omega);
  trackState->setD0(d0);
  trackState->setZ0(z0);

  // Uncertainties (covariance matrix)
  Acts::ActsMatrix<6, 6> jac = Acts::ActsMatrix<6, 6>::Zero();

  jac(0, Acts::eBoundLoc0) = 1;

  jac(1, Acts::eBoundPhi) = 1;

  jac(2, Acts::eBoundTheta) = omega / std::tan(theta);
  jac(2, Acts::eBoundQOverP) = omega / qoverp;

  jac(3, Acts::eBoundLoc1) = 1;

  jac(4, Acts::eBoundTheta) = std::pow(1 / std::cos(lambda), 2);

  Acts::ActsMatrix<6, 6> trcov = (jac * cov * jac.transpose());

  EVENT::FloatVec lcioCov(15, 0);
  lcioCov[0] = trcov(0, 0);
  lcioCov[1] = trcov(0, 1);
  lcioCov[2] = trcov(1, 1);
  lcioCov[3] = trcov(0, 2);
  lcioCov[4] = trcov(1, 2);
  lcioCov[5] = trcov(2, 2);
  lcioCov[6] = trcov(0, 3);
  lcioCov[7] = trcov(1, 3);
  lcioCov[8] = trcov(2, 3);
  lcioCov[9] = trcov(3, 3);
  lcioCov[10] = trcov(0, 4);
  lcioCov[11] = trcov(1, 4);
  lcioCov[12] = trcov(2, 4);
  lcioCov[13] = trcov(3, 4);
  lcioCov[14] = trcov(4, 4);

  trackState->setCovMatrix(lcioCov);

  return trackState;
}

EVENT::LCCollection* getCollection(EVENT::LCEvent* evt,
                                   const std::string& name) {
  if (name.size() == 0) return nullptr;

  try {
    return evt->getCollection(name);
  } catch (const EVENT::DataNotAvailableException& e) {
    // TODO: Reenable output
    // streamlog_out( DEBUG2 ) << "getCollection :  DataNotAvailableException :
    // " << name <<  std::endl ;
    return nullptr;
  }
}

Acts::ParticleHypothesis convertParticle(const EVENT::MCParticle* mcParticle)
{
  switch (mcParticle->getPDG()) {
  case 11:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eElectron};
  case -11:
    return Acts::ParticleHypothesis {Acts::PdgParticle::ePositron};
  case 13:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eMuon};
  case -13:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eAntiMuon};
  case 15:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eTau};
  case -15:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eAntiTau};
  case 22:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eGamma};
  case 111:
    return Acts::ParticleHypothesis {Acts::PdgParticle::ePionZero};
  case 211:
    return Acts::ParticleHypothesis {Acts::PdgParticle::ePionPlus};
  case -211:
    return Acts::ParticleHypothesis {Acts::PdgParticle::ePionMinus};
  case 2112:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eNeutron};
  case -2112:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eAntiNeutron};
  case 2212:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eProton};
  case -2212:
    return Acts::ParticleHypothesis {Acts::PdgParticle::eAntiProton};
  }

  Acts::PdgParticle pdg = Acts::PdgParticle::eInvalid;
  float mass = 0.0f;
  Acts::AnyCharge charge_type { 0.0f };
  return Acts::ParticleHypothesis { pdg, mass, charge_type };
}


}  // namespace ACTSTracking
