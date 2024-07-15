#include "GeometryIdMappingTool.hxx"

#include <iomanip>

using namespace ACTSTracking;

const int32_t GeometryIdMappingTool::VertexEndCapNegative = -2;
const int32_t GeometryIdMappingTool::VertexBarrel = 1;
const int32_t GeometryIdMappingTool::VertexEndCapPositive = 2;
const int32_t GeometryIdMappingTool::InnerTrackerEndCapNegative = -4;
const int32_t GeometryIdMappingTool::InnerTrackerBarrel = 3;
const int32_t GeometryIdMappingTool::InnerTrackerEndCapPositive = 4;
const int32_t GeometryIdMappingTool::OuterInnerTrackerEndCapNegative = -8;
const int32_t GeometryIdMappingTool::OuterInnerTrackerBarrel = 7;
const int32_t GeometryIdMappingTool::OuterInnerTrackerEndCapPositive = 8;
const int32_t GeometryIdMappingTool::OuterTrackerEndCapNegative = -6;
const int32_t GeometryIdMappingTool::OuterTrackerBarrel = 5;
const int32_t GeometryIdMappingTool::OuterTrackerEndCapPositive = 6;

const std::unordered_map<uint32_t, uint32_t>
    GeometryIdMappingTool::NLad_VertexBarrel = {
        {0, 5}, {1, 5}, {2, 5}, {3, 5}, {4, 5}};

const std::unordered_map<uint32_t, uint32_t>
    GeometryIdMappingTool::NRng_VertexEndCap = {
        {0, 16}, {1, 16}, {2, 16}, {3, 16}, {4, 16}, {5, 16}, {6, 16}, {7, 16}};

const std::unordered_map<uint32_t, uint32_t>
    GeometryIdMappingTool::NLad_InnerTrackerBarrel = {
        {0, 32},
        {1, 32},
};

const std::unordered_map<uint32_t, uint32_t>
    GeometryIdMappingTool::NRng_InnerTrackerEndCap = {
        {0, 26},
};

const std::unordered_map<uint32_t, uint32_t>
    GeometryIdMappingTool::NLad_OuterInnerTrackerBarrel = {{2, 46}};

const std::unordered_map<uint32_t, uint32_t>
    GeometryIdMappingTool::NRng_OuterInnerTrackerEndCap = {
        {1, 26}, {2, 26}, {3, 26}, {4, 26}, {5, 26}, {6, 26}};

const std::unordered_map<uint32_t, uint32_t>
    GeometryIdMappingTool::NLad_OuterTrackerBarrel = {
        {0, 84}, {1, 84}, {2, 84}};

const std::unordered_map<uint32_t, uint32_t>
    GeometryIdMappingTool::NRng_OuterTrackerEndCap = {
        {0, 48}, {1, 48}, {2, 48}, {3, 48}, {4, 48}, {5, 48}, {6, 48}, {7, 48}};

const std::unordered_map<int32_t, uint32_t> GeometryIdMappingTool::VolumeMap = {
    {VertexEndCapNegative, 13},
    {VertexBarrel, 14},
    {VertexEndCapPositive, 15},
    {InnerTrackerEndCapNegative, 8},
    {InnerTrackerBarrel, 17},
    {InnerTrackerEndCapPositive, 18},
    {OuterInnerTrackerEndCapNegative, 3},
    {OuterInnerTrackerBarrel, 20},
    {OuterInnerTrackerEndCapPositive, 21},
    {OuterTrackerEndCapNegative, 23},
    {OuterTrackerBarrel, 24},
    {OuterTrackerEndCapPositive, 25},
};

GeometryIdMappingTool::GeometryIdMappingTool(const std::string& encoderString)
    : _encoderString(encoderString) {}

uint64_t GeometryIdMappingTool::getGeometryID(const lcio::SimTrackerHit* hit) {
  UTIL::CellIDDecoder<lcio::SimTrackerHit> decoder(_encoderString);
  uint32_t systemID = decoder(hit)["system"];
  uint32_t layerID = decoder(hit)["layer"];
  int32_t sideID = decoder(hit)["side"];
  uint32_t ladderID = decoder(hit)["module"];
  uint32_t moduleID = decoder(hit)["sensor"];

  return getGeometryID(systemID, layerID, sideID, ladderID, moduleID);
}

uint64_t GeometryIdMappingTool::getGeometryID(const lcio::TrackerHit* hit) {
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
  uint32_t systemID = decoder(hit)["system"];
  uint32_t layerID = decoder(hit)["layer"];
  int32_t sideID = decoder(hit)["side"];
  uint32_t ladderID = decoder(hit)["module"];
  uint32_t moduleID = decoder(hit)["sensor"];

  return getGeometryID(systemID, layerID, sideID, ladderID, moduleID);
}

uint64_t GeometryIdMappingTool::getGeometryID(uint32_t systemID,
                                              uint32_t layerID, int32_t sideID,
                                              uint32_t ladderID,
                                              uint32_t moduleID) {
  uint64_t geometry_id = 0;

  //
  // Volume ID determination.

  // the outermost layer of InnerTracker is "OuterInnerTracker" in ACTS
  if (systemID == InnerTrackerBarrel && layerID == 2) {
    systemID = OuterInnerTrackerBarrel;
  }
  if (systemID == InnerTrackerEndCapPositive && layerID != 0) {
    systemID = OuterInnerTrackerEndCapPositive;
  }

  // endcap is split in +/- sides by ACTS
  int32_t signSystemID = (sideID < 0) ? -systemID : systemID;

  // Map
  uint64_t volume_id = (VolumeMap.find(signSystemID) != VolumeMap.end())
                           ? VolumeMap.at(signSystemID)
                           : signSystemID;
  geometry_id |= volume_id << (14 * 4);

  // Layer ID is counting within sub detector, with pairings depending on the
  // sub detector
  uint64_t layer_id;
  switch (signSystemID) {
    case VertexEndCapNegative:
      layer_id = 2 * (7 - layerID) + 2;
      break;
    case VertexEndCapPositive:
      layer_id = 2 * (layerID) + 2;
      break;
    case VertexBarrel:
      switch (layerID) {
        case 0:
          layer_id = 2;
        case 1:
          layer_id = 4;
        case 2:
          layer_id = 6;
        case 4:
          layer_id = 8;
        case 6:
          layer_id = 10;
      }
      break;
    case InnerTrackerBarrel:
    case OuterTrackerBarrel: {
      layer_id = 2 * layerID + 2;
    } break;
    case OuterInnerTrackerBarrel: {
      layer_id = 2 * (2 - layerID) + 2;
    } break;

    case InnerTrackerEndCapNegative: {
      layer_id = 2 * (0 - layerID) + 2;
    } break;
    case InnerTrackerEndCapPositive: {
      layer_id = 2 * layerID + 2;
    } break;

      // OuterInner tracker layer counting starts at layer 1
      // as it is layer 1 of the original Inner tracker.
    case OuterInnerTrackerEndCapNegative: {
      layer_id = 2 * (7 - layerID) + 0;
    } break;
    case OuterInnerTrackerEndCapPositive: {
      layer_id = 2 * layerID + 0;
    } break;

    case OuterTrackerEndCapNegative: {
      layer_id = 2 * (3 - layerID) + 2;
    } break;
    case OuterTrackerEndCapPositive: {
      layer_id = 2 * layerID + 2;
    } break;

    default:
      layer_id = layerID;
      break;
  }
  geometry_id |= layer_id << (9 * 4);

  // Module ID counting depends on sub detector and layer
  uint64_t sensitive_id;
  switch (signSystemID) {
    case VertexBarrel:
      sensitive_id = NLad_VertexBarrel.at(layerID) * ladderID + moduleID + 1;
      break;
    case InnerTrackerBarrel:
      sensitive_id =
          NLad_InnerTrackerBarrel.at(layerID) * ladderID + moduleID + 1;
      break;
    case OuterInnerTrackerBarrel:
      sensitive_id =
          NLad_OuterInnerTrackerBarrel.at(layerID) * ladderID + moduleID + 1;
      break;
    case OuterTrackerBarrel:
      sensitive_id =
          NLad_OuterTrackerBarrel.at(layerID) * ladderID + moduleID + 1;
      break;

    case VertexEndCapNegative:
      sensitive_id = 16 - moduleID;
      break;
    case VertexEndCapPositive:
      sensitive_id = moduleID + 1;
      break;

    case InnerTrackerEndCapNegative:
    case InnerTrackerEndCapPositive:
      sensitive_id =
          NRng_InnerTrackerEndCap.at(layerID) * ladderID + moduleID + 1;
      break;

    case OuterInnerTrackerEndCapNegative:
    case OuterInnerTrackerEndCapPositive:
      sensitive_id =
          NRng_OuterInnerTrackerEndCap.at(layerID) * ladderID + moduleID + 1;
      break;

    case OuterTrackerEndCapNegative:
    case OuterTrackerEndCapPositive:
      sensitive_id =
          NRng_OuterTrackerEndCap.at(layerID) * ladderID + moduleID + 1;
      break;

    default:
      sensitive_id = ladderID;
      break;
  }
  geometry_id |= sensitive_id << (2 * 4);

  return geometry_id;
}
