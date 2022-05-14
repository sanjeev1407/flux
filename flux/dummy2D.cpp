// This file is part of the Finite-element soLver for Unsteady electromagnetiX
// (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of
// Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//============================================================================

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

#include <mfem.hpp>

#include "flux/debug_assert.hpp"
#include "flux/precice_adapter.hpp"

using namespace mfem;
using namespace std;

namespace details {
  /*std::vector<int> BuildPreciceToMfemMap(
    ParMesh* mesh, std::function<bool(Vector)> acceptable_check);*/
std::vector<int> BuildPreciceToMfemMap(
    ParMesh* mesh, const double x, const double r);
  
std::unordered_map<int, int> BuildMfemToPreciceMap(
    const std::vector<int>& precice_to_mfem_map);

std::vector<double> BuildPreciceVertexPositions(
    ParMesh* mesh, const std::vector<int>& precice_to_mfem_map);
}  // namespace details

int main(int argc, char* argv[]) {
  // Initialize MPI
  MPI_Session mpi(argc, argv);
  int myid = mpi.WorldRank();

  // Parse command-line options.
  const char* mesh_file = "beam-tet.mesh";
  int order = 1;
  const char* precice_name_in = "";
  const char* psolver_name_in = "DummySolver";
  double cutoff_radius = -1.0;
  double cutoff_length = -1.0;
  const char* basename = "dummy";
  int time_steps = 1;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&time_steps, "-t", "--time-steps",
                 "Number of time steps to take");

  // preCICE related
  args.AddOption(&precice_name_in, "-p", "--precice-config",
                 "Name of preCICE configuration file.");
  args.AddOption(&psolver_name_in, "-s", "--solver-name",
                 "Name of solver in preCICE");
  args.AddOption(&cutoff_radius, "-cr", "--cutoff-radius",
                 "Radius below which sigma will be passed via preCICE");
  args.AddOption(&cutoff_length, "-cl", "--cutoff-length",
                 "Length below which sigma will be passed via preCICE");
  args.Parse();

  const auto precice_name = std::string(precice_name_in);
  const auto psolver_name = std::string(psolver_name_in);
  const bool precice_active = precice_name.size() != 0;
  // Make sure if preCICE is being used, that a cutoff radius is specified
  DEBUG_ASSERT(
      precice_active, global_assert{}, DebugLevel::CHEAP{},
      "Dummy2D exists for the sole purpose of being used with preCICE");
  DEBUG_ASSERT(time_steps > 0, global_assert{}, DebugLevel::CHEAP{},
               "Must take atleast one time-step (time_steps > 0)");

  // 3. Read the mesh from the given mesh file. We can handle triangular,
  //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
  //    the same code.
  Mesh* mesh = new Mesh(mesh_file, 1, 1);
  const int dim = mesh->Dimension();
  int sdim = mesh->SpaceDimension();

  // 8. Define a parallel mesh by a partitioning of the serial mesh.
  //    parallel mesh is defined, the serial mesh can be deleted.
  ParMesh* pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;  // Do not need after building the ParMesh

  // Necessary when using Nedelec spaces with order > 1
  // for a tetrahedral mesh
  pmesh->ReorientTetMesh();

  // Defining all the Finite element collections and spaces to be used
  H1_FECollection HGradFEC_order1(1, dim);
  ParFiniteElementSpace HGradFESpace_order1(pmesh, &HGradFEC_order1);

  H1_FECollection HGradFEC(order, dim);
  ParFiniteElementSpace HGradFESpace(pmesh, &HGradFEC);

  ND_FECollection HCurlFEC(order, dim);
  ParFiniteElementSpace HCurlFESpace(pmesh, &HCurlFEC);

  RT_FECollection HDivFEC(order, dim);
  ParFiniteElementSpace HDivFESpace(pmesh, &HDivFEC);

  L2_FECollection L2FEC_order0(0, dim);
  ParFiniteElementSpace L2FESpace_order0(pmesh, &L2FEC_order0);
  ParFiniteElementSpace L2FESpace_order0_vec(pmesh, &L2FEC_order0, dim,
                                             Ordering::byVDIM);

  // Activate precice and add mesh/data

  const std::string precice_mesh_name = "DummyMesh";
  std::vector<int> precice_to_mfem_map;
  /*precice_to_mfem_map = details::BuildPreciceToMfemMap(
    pmesh, [cutoff_radius](Vector x) { return x(1) <= cutoff_radius; });*/
  precice_to_mfem_map = details::BuildPreciceToMfemMap(
        pmesh, cutoff_length, cutoff_radius);
  std::unordered_map<int, int> mfem_to_precice_map =
      details::BuildMfemToPreciceMap(precice_to_mfem_map);

  PreciceAdapter precice_adapter(psolver_name, precice_name, myid,
                                 mpi.WorldSize());
  precice_adapter.AddMesh(precice_mesh_name);
  std::vector<double> centroid_positions =
      details::BuildPreciceVertexPositions(pmesh, precice_to_mfem_map);
  precice_adapter.SetVertexPositions(precice_mesh_name, centroid_positions);
  precice_adapter.AddData("Sigma", precice_mesh_name, DataOperation::WRITE);
  precice_adapter.AddData("JouleHeating", precice_mesh_name,
                          DataOperation::READ);
  precice_adapter.AddData("LorentzForce", precice_mesh_name,
                          DataOperation::READ);

  double time = 0.0;
  auto sigmaFunc = [&time](Vector x) {
    std::array<double, 2> lower_left{{0.0, 0.0}};
    std::array<double, 2> upper_right{{0.5, 0.08}};
    std::array<double, 2> length{
        {upper_right[0] - lower_left[0], upper_right[1] - lower_left[1]}};
    const double amplitude = 1.0e3 * (time + 1.0);
    x(0) = (x(0) - lower_left[0]) / length[0];
    x(1) = (x(1) - lower_left[1]) / length[1];
    return amplitude * std::sqrt(x(0) * x(0) + x(1) * x(1));
  };
  FunctionCoefficient sigma_coef(sigmaFunc);

  // Function coefficient for sigma

  // Projecting the sigma data (possibly cell-center data) on zeroth-order
  // L2 finite element space
  ParGridFunction sigma_L2_order0_gf(&L2FESpace_order0);
  GridFunctionCoefficient sigma_L2_order0_gf_coef(&sigma_L2_order0_gf);

  ParGridFunction joule_heating(&L2FESpace_order0);
  ParGridFunction lorentz_force(&L2FESpace_order0_vec);

  double dt = precice_adapter.Initialize();

  Vector data;

  ParaViewDataCollection paraview_dc(basename, pmesh);
  paraview_dc.SetPrefixPath("Dummy_ParaView");
  paraview_dc.SetLevelsOfDetail(order);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.RegisterField("Sigma", &sigma_L2_order0_gf);
  paraview_dc.RegisterField("JouleHeating", &joule_heating);
  paraview_dc.RegisterField("LorentzForce", &lorentz_force);
  for (int ts = 0; ts < time_steps; ++ts) {
    time = static_cast<double>(ts);
    sigma_L2_order0_gf.ProjectCoefficient(sigma_coef);

    data.SetSize(static_cast<int>(precice_to_mfem_map.size()));
    for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
      data(i) = sigma_L2_order0_gf(precice_to_mfem_map[i]);
    }
    const double* write_ptr = static_cast<const double*>(data);
    precice_adapter.WriteBlockScalarData("Sigma", write_ptr);

    double precice_dt = precice_adapter.Advance(1.0);

    // Read preCICE data
    data.SetSize(static_cast<int>(precice_to_mfem_map.size()));
    double* read_ptr = static_cast<double*>(data);
    precice_adapter.ReadBlockScalarData("JouleHeating", read_ptr);
    for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
      joule_heating(precice_to_mfem_map[i]) = data(i);
    }

    data.SetSize(dim * static_cast<int>(precice_to_mfem_map.size()));
    read_ptr = static_cast<double*>(data);
    precice_adapter.ReadBlockVectorData("LorentzForce", read_ptr);
    for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
      for (int d = 0; d < dim; ++d) {
        lorentz_force(dim * precice_to_mfem_map[i] + d) = data(dim * i + d);
      }
    }

    // Export to preCICE
    paraview_dc.SetCycle(ts);
    paraview_dc.SetTime(static_cast<double>(ts));
    paraview_dc.Save();
  }
  precice_adapter.Finalize();

  delete pmesh;
  return 0;
}

namespace details {
  /*std::vector<int> BuildPreciceToMfemMap(
    ParMesh* mesh, std::function<bool(Vector)> acceptable_check) {*/
  std::vector<int> BuildPreciceToMfemMap(
  ParMesh* mesh, const double x, const double r) {
  std::vector<int> precice_to_mfem_map;
  mfem::Vector center;
  for (int i = 0; i < mesh->GetNE(); ++i) {
    mesh->GetElementCenter(i, center);
    if (center(1) <= r) {
      if ((center(0) >= 0.) && (center(0) <= x)) {
         precice_to_mfem_map.push_back(i);
      }

      /*if (acceptable_check(center)) {
	precice_to_mfem_map.push_back(i);*/
    }
  }
  return precice_to_mfem_map;
}

std::unordered_map<int, int> BuildMfemToPreciceMap(
    const std::vector<int>& precice_to_mfem_map) {
  std::unordered_map<int, int> mfem_to_precice_map;
  for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
    const auto insert_success =
        mfem_to_precice_map.insert({precice_to_mfem_map[i], i}).second;
    DEBUG_ASSERT(
        insert_success, global_assert{}, DebugLevel::CHEAP{},
        "An element index is included twice in the precice_to_mfem_map");
  }
  return mfem_to_precice_map;
}

std::vector<double> BuildPreciceVertexPositions(
    ParMesh* mesh, const std::vector<int>& precice_to_mfem_map) {
  const std::size_t dim = mesh->Dimension();
  mfem::Vector center(static_cast<int>(dim));
  const auto nvert = precice_to_mfem_map.size();
  std::vector<double> contiguous_positions(dim * nvert);
  for (std::size_t i = 0; i < nvert; ++i) {
    mesh->GetElementCenter(precice_to_mfem_map[i], center);
    for (std::size_t d = 0; d < dim; ++d) {
      contiguous_positions[i * dim + d] = center(static_cast<int>(d));
    }
  }
  return contiguous_positions;
}
}  // namespace details
