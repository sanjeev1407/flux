
// This file is part of the Finite-element soLver for Unsteady electromagnetiX
// (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of
// Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=================================================================================

#ifndef FLUX_FLUX_2D_HPP_
#define FLUX_FLUX_2D_HPP_

#include <functional>
#include <unordered_map>
#include <vector>

#include <mfem.hpp>

#include "flux/debug_assert.hpp"
#include "flux/precice_adapter.hpp"

using namespace mfem;

// Builds contiguous list of acceptable elements in Mfem for
// volume coupling. Here, acceptable is defined by
// the provided lambda. Returned is a vector of
// integers mapping precice vertex (the vector index)
// to Mfem element indices (the vector value at that index).
//
// NOTE: The std::function is given a location as an mfem::Vector
// which must then be deemed acceptable (true) or not (false).
/*std::vector<int> BuildPreciceToMfemMap(
  ParMesh* mesh, std::function<bool(Vector)> acceptable_check);*/
std::vector<int> BuildPreciceToMfemMap(
    ParMesh* mesh, const double x1, const double x2, const double r);

// Inverts the precice_to_mfem_map and stores it in an
// unordered_map as a lookup table. Using the returned table
// as table[mfem_index] returns the contiguous index used for preCICE.
std::unordered_map<int, int> BuildMfemToPreciceMap(
    const std::vector<int>& precice_to_mfem_map);

// Build list of vertex positions (centroids) for preCICE.
std::vector<double> BuildPreciceVertexPositions(
    ParMesh* mesh, const std::vector<int>& precice_to_mfem_map);

// Starting from L2 0th order space, project to H1 and then
// successively increase order to final order.
void SuccessiveRefinementInterpolation(
    const ParGridFunction& starting_l2, ParGridFunction& final_h1,
    std::vector<ParGridFunction>& working_gf);

// Class for using precice volume coupled data as a coefficient.
// ONLY TO BE USED FOR L2 0th order data.
// If location is in coupling region, returns
// the value from preCICE, otherwise returns
// a default value provided in the constructor.
class VolumeCouplingCoefficient : public Coefficient {
 private:
  const PreciceAdapter& precice_adapter;
  const std::string precice_name;
  const std::unordered_map<int, int>* mfem_to_precice_map;
  Vector data;
  std::function< double(const Vector &)> sigmaFunc;
  
 public:
  VolumeCouplingCoefficient(
      const PreciceAdapter& a_precice_adapter,
      const std::string& a_precice_name,
      const std::unordered_map<int, int>* a_mfem_to_precice_map,
      std::function< double(const Vector &)> a_sigmaFunc)
      : precice_adapter(a_precice_adapter),
        precice_name(a_precice_name),
        mfem_to_precice_map(a_mfem_to_precice_map),
        sigmaFunc(a_sigmaFunc) {
    DEBUG_ASSERT(mfem_to_precice_map != nullptr, global_assert{},
                 DebugLevel::CHEAP{});
    data.SetSize(mfem_to_precice_map->size());
  }

  void UpdateData(void) {
    // To be explicit about the conversion
    double* data_ptr = static_cast<double*>(data);
    precice_adapter.ReadBlockScalarData(precice_name, data_ptr);
  }

  virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip) {
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    auto loc = mfem_to_precice_map->find(T.ElementNo);
    return loc != mfem_to_precice_map->end() ? data(loc->second)
                                             : sigmaFunc(transip);
  }

  virtual ~VolumeCouplingCoefficient(void) = default;
};

using namespace mfem;

#endif  // FLUX_FLUX_2D_HPP_
