#pragma once

#include <EVENT/LCEvent.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/TrackState.h>

#include <Acts/EventData/TrackParameters.hpp>
#include "Acts/EventData/ParticleHypothesis.hpp"
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

#include "SourceLink.hxx"

namespace ACTSTracking {

using TrackResult = Acts::TrackContainer<Acts::VectorTrackContainer,
                                         Acts::VectorMultiTrajectory,
                                         std::shared_ptr>::TrackProxy;

//! Get path to a resource file
/**
 * Get absolute file of a file `inpath` by looking in the following places:
 *  - `inpath` to the current working directory
 *  - `ACTSTRACKING_SOURCEDIR/inpath`
 *  - `ACTSTRACKING_DATADIR/inpath`
 *
 * If the files is not found at any location, then `inpath` is returned.
 * If `path` starts with a /, then it is returned directly.
 *
 * \parm inpath File to find.
 *
 * \return Absolute path to file.
 */
std::string findFile(const std::string& inpath);

//! Convert ACTS CKF result to LCIO track class
/**
 * Converted propertie are:
 *  - goodness of fit (chi2, ndf)
 *  - associated hits
 *  - track states at IP
 *
 * \param fitOutput CKF fit result
 * \param trackTip index of track to convert inside fit result
 * \param magneticField magnetic field at different locations in the detector
 * \param magCache cache to help with magnetic field lookup
 *
 * \return Track with equivalent parameters of the ACTS track
 */
/*
EVENT::Track* ACTS2Marlin_track(
    const Acts::CombinatorialKalmanFilterResult<ACTSTracking::SourceLink>&
        fitOutput,
    std::size_t trackTip,
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
    Acts::MagneticFieldProvider::Cache& magCache);
*/
//! Convert ACTS KF result to LCIO track class
/**
 * Converted propertie are:
 *  - goodness of fit (chi2, ndf)
 *  - associated hits
 *  - track states at IP
 *
 * \param fitOutput KF fit result
 * \param magneticField magnetic field at different locations in the detector
 * \param magCache cache to help with magnetic field lookup
 *
 * \return Track with equivalent parameters of the ACTS track
 */
/*
EVENT::Track* ACTS2Marlin_track(
    const Acts::KalmanFitterResult<ACTSTracking::SourceLink>& fitOutput,
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
    Acts::MagneticFieldProvider::Cache& magCache);
*/

EVENT::Track* ACTS2Marlin_track(
    const TrackResult& fitter_res,
    std::shared_ptr<Acts::MagneticFieldProvider> magneticField,
    Acts::MagneticFieldProvider::Cache& magCache,
    double caloFaceR, double caloFaceZ, Acts::GeometryContext geoContext,
    Acts::MagneticFieldContext magFieldContext,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeo);

//! Convert ACTS track state class to Marlin class
/**
 * \param location Location where the track state is defined (ie: `AtIP`)
 * \param params ACTS track state parameters
 * \params Bz magnetic field at location of track state [Tesla]
 *
 * \return Track state with equivalent parameters of the ACTS track
 */
EVENT::TrackState* ACTS2Marlin_trackState(
    int location, const Acts::BoundTrackParameters& params, double Bz);

EVENT::TrackState* ACTS2Marlin_trackState(int location,
                                          const Acts::BoundVector& value,
                                          const Acts::BoundMatrix& cov,
                                          double Bz);

//! Get collection from `LCEvent` with silent fail
/**
 * \param evt event store
 * \param event collection name
 *
 * \return Collection, if found, `nullptr` otherwise
 */
EVENT::LCCollection* getCollection(EVENT::LCEvent* evt,
                                   const std::string& name);

Acts::ParticleHypothesis convertParticle(const EVENT::MCParticle* mcParticle);

}  // namespace ACTSTracking

