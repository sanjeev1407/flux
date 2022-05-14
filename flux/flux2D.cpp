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

#include "flux/flux2D.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

#include <mfem.hpp>

#include "flux/debug_assert.hpp"
#include "flux/flux2D_solver.hpp"
#include "flux/mfem_coefficients.hpp"
#include "flux/precice_adapter.hpp"

using namespace mfem;
using namespace std;


// Physical constants
double w = 0.0;  // Angular frequency (rad/s)
double sigmaCoil = 6.0e7;        // Electrical conductivity of coil (copper)
double sigmaAir = 1.0e-6;        // Electrical conductivity of air
double mu = 1.2566e-6;       // Magnetic Permeability
double cls = 0.0;           //cut-off length start
double cle = 0.0;           //cut-off length end
double cr = 0.0;            //cut-off radius
int dim;


int main(int argc, char* argv[]) {
  // Initialize MPI
  MPI_Session mpi(argc, argv);
  int myid = mpi.WorldRank();

  // Parse command-line options.
  const char* mesh_file = "beam-tet.mesh";
  int order = 1;
  const char* precice_name_in = "";
  const char* psolver_name_in = "FLUX";
  double cutoff_radius = -1.0;
  double cutoff_length_start = -1.0;
  double cutoff_length_end = -1.0;
  const char* basename = "flux2D";
  double precice_dt = -1.0;
  int print_interval = 1;
  double P_target = 100000.0;
  double freq = 0.37e6;
  double nb_coils = 6;
  double coil_radius = 0.109;
  Vector coil_loc;

  
  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");

  // preCICE related
  args.AddOption(&precice_name_in, "-p", "--precice-config",
                 "Name of preCICE configuration file.");
  args.AddOption(&psolver_name_in, "-s", "--solver-name",
                 "Name of solver in preCICE");
  args.AddOption(&precice_dt, "-pdt", "--precice_dt",
                 "Window size used for preCICE");
  args.AddOption(&cutoff_radius, "-cr", "--cutoff-radius",
                 "Radius below which sigma will be passed via preCICE");
  args.AddOption(&cutoff_length_start, "-cls", "--cutoff-length-start",
                 "Starting location of cutoff length where sigma will be passed via preCICE");
  args.AddOption(&cutoff_length_end, "-cle", "--cutoff-length-end",
                 "End location od cutoff length where sigma will be passed via preCICE");
  args.AddOption(&print_interval, "-pr", "--print-interval",
                 "Interval after which flux output will be saved");
  args.AddOption(&P_target, "-power", "--P_target",
                 "Target power for the simulation in Watts");
  args.AddOption(&freq, "-f", "--freq",
                 "Inductor frequency in Hz");
  args.AddOption(&nb_coils, "-nc", "--nb_coils",
                 "Number of coils");
  args.AddOption(&coil_radius, "-co_r", "--coil_radius",
                 "Coil radius in metres");
  args.AddOption(&coil_loc, "-co_loc", "--coil_loc",
                 "Coil locations in metres");

  args.Parse();

  if (coil_loc.Size() != nb_coils) {
    cout << "Number of inputs for coil locations does not match with nb_coils" << endl;
    exit(0);
  }

  //  Angular frequency
  w = 2.0*M_PI*freq;

  // cut-off length and radius
  cls = cutoff_length_start;
  cle = cutoff_length_end;
  cr = cutoff_radius;
  
  const auto precice_name = std::string(precice_name_in);
  const auto psolver_name = std::string(psolver_name_in);
  const bool precice_active = precice_name.size() != 0;
  // Make sure if preCICE is being used, that a cutoff radius is specified
  DEBUG_ASSERT(!precice_active || cutoff_radius > 0.0, global_assert{},
               DebugLevel::CHEAP{},
               "Cutoff radius must be specified when preCICE is being used.");
  DEBUG_ASSERT(!precice_active || cutoff_length_start > 0.0, global_assert{},
               DebugLevel::CHEAP{},
               "Cutoff length start point must be specified when preCICE is being used.");
  DEBUG_ASSERT(!precice_active || cutoff_length_end > 0.0, global_assert{},
               DebugLevel::CHEAP{},
               "Cutoff length end point must be specified when preCICE is being used.");
  DEBUG_ASSERT(!precice_active || precice_dt > 0, global_assert{},
               DebugLevel::CHEAP{},
               "Time window length must be given if using preCICE");
   
  // 3. Read the mesh from the given mesh file. We can handle triangular,
  //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
  //    the same code.

  Mesh* mesh = new Mesh(mesh_file, 1, 1);
  dim = mesh->Dimension();
  int sdim = mesh->SpaceDimension();

  // 8. Define a parallel mesh by a partitioning of the serial mesh.
  //    parallel mesh is defined, the serial mesh can be deleted.
  ParMesh* pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;  // Do not need after building the ParMesh

  // Necessary when using Nedelec spaces with order > 1
  // for a tetrahedral mesh
  pmesh->ReorientTetMesh();

  // Defining all the Finite element collections and spaces to be used

  // Generate all orders of H1 from 1 to order. Needed for successive
  // interpolation of L2 sigma to correct order H1 space.
  std::vector<H1_FECollection*> H1_collections(order);
  std::vector<ParFiniteElementSpace*> H1_fes(order);
  std::vector<ParGridFunction> working_gf;
  for (int s = 1; s <= order; ++s) {
    H1_collections.emplace_back(new H1_FECollection(s, dim));
    H1_fes.emplace_back(
        new ParFiniteElementSpace(pmesh, H1_collections.back()));
    working_gf.emplace_back(H1_fes.back());
  }

  ParFiniteElementSpace& HGradFESpace_order1 = *H1_fes.front();
  ParFiniteElementSpace& HGradFESpace = *H1_fes.back();

  ND_FECollection HCurlFEC(order+1, dim);
  ParFiniteElementSpace HCurlFESpace(pmesh, &HCurlFEC);

  RT_FECollection HDivFEC(order, dim);
  ParFiniteElementSpace HDivFESpace(pmesh, &HDivFEC);

  L2_FECollection L2FEC_order0(0, dim);
  ParFiniteElementSpace L2FESpace_order0(pmesh, &L2FEC_order0);
  ParGridFunction tmp_gf_L20(&L2FESpace_order0);

  // Activate precice, if needed
  std::vector<int> precice_to_mfem_map;
  std::unordered_map<int, int> mfem_to_precice_map;
  mfem::Coefficient* sigma_coef = nullptr;
  PreciceAdapter* precice_adapter = nullptr;
  std::string precice_mesh_name = "";
  Vector write_data;
  if (precice_active) {
    precice_mesh_name = "flux_volume_mesh";
    precice_to_mfem_map = BuildPreciceToMfemMap(pmesh, cutoff_length_start, cutoff_length_end, cutoff_radius);
    mfem_to_precice_map = BuildMfemToPreciceMap(precice_to_mfem_map);

    precice_adapter =
        new PreciceAdapter(psolver_name, precice_name, myid, mpi.WorldSize());
    precice_adapter->AddMesh(precice_mesh_name);
    std::vector<double> centroid_positions =
        BuildPreciceVertexPositions(pmesh, precice_to_mfem_map);
    precice_adapter->SetVertexPositions(precice_mesh_name, centroid_positions);
    precice_adapter->AddData("Sigma", precice_mesh_name, DataOperation::READ);
    precice_adapter->AddData("JouleHeating", precice_mesh_name,
                             DataOperation::WRITE);
    precice_adapter->AddData("AxialLorentzForce", precice_mesh_name,
			     DataOperation::WRITE);
    precice_adapter->AddData("RadialLorentzForce", precice_mesh_name,
			     DataOperation::WRITE);

    sigma_coef = new VolumeCouplingCoefficient(*precice_adapter, "Sigma",
                                               &mfem_to_precice_map, sigmaFunc);

  } else {
    sigma_coef = new FunctionCoefficient(sigmaFunc);
  }

  // ===========================================================
  // Initialize as much as possible outside of initial loop.
  // ===========================================================

  // Projecting the sigma data (possibly cell-center data) on zeroth-order
  // L2 finite element space
  ParGridFunction sigma_L2_order0_gf(&L2FESpace_order0);
  GridFunctionCoefficient sigma_L2_order0_gf_coef(&sigma_L2_order0_gf);

  ParGridFunction sigma_gf(&HGradFESpace);
  GridFunctionCoefficient sigma_gf_coef(&sigma_gf);

  // ===========================================================
  // Initialization for E Computation
  // ===========================================================
  ComplexOperator::Convention conv = ComplexOperator::HERMITIAN;

  // Assign the boundary conditions
  Array<int> E_ess_bdr_zero(pmesh->bdr_attributes.Max());

  E_ess_bdr_zero = 1;
  
  HYPRE_Int glob_size_h1 = HGradFESpace.GlobalTrueVSize();

  if (mpi.Root()) {
    cout << "Number of Electric Field unknowns:    " << glob_size_h1 <<  endl;
  }

  //    Set up the parallel complex linear form b(.) which corresponds to the
  //    right-hand side -i < w mu r Jc, v> of the FEM linear system.

  ConstantCoefficient Zero(0.0);

  ParLinearForm* b_re = new ParLinearForm(&HGradFESpace);
  ParLinearForm* b_im = new ParLinearForm(&HGradFESpace);

  ParComplexLinearForm b(&HGradFESpace, b_re, b_im, conv);

  //FunctionCoefficient source(current_source);  
  //b.AddDomainIntegrator(new DomainLFIntegrator(Zero), new DomainLFIntegrator(source));

  // Point current source using delta coefficient
  const double I0 = -100.0;
  const double s = w * mu * coil_radius * I0;
  std::vector<DeltaCoefficient*> delta_current;
  for (int i = 0; i < nb_coils; i++) {
    delta_current.emplace_back(new DeltaCoefficient(coil_loc[i], coil_radius, s));
    b.AddDomainIntegrator(new DomainLFIntegrator(Zero),
			  new DomainLFIntegrator(*delta_current.back()));
  }
  
  
  // Define the solution vector E as a parallel complex finite element grid
  // function corresponding to fespace
  ParComplexGridFunction E(&HGradFESpace);

  // Initialize E to 0 + i0
  E.ProjectCoefficient(Zero, Zero);

  // E = 0 Essential BC
  E.real().ProjectBdrCoefficient(Zero, E_ess_bdr_zero);
  E.imag().ProjectBdrCoefficient(Zero, E_ess_bdr_zero);


  // Set up the parallel sesquilinear form a(.,.) on the finite element
  // space corrsponding to <r grad u, grad v> + <1/r u, v> + i<w mu r sigma u, v>
  FunctionCoefficient r_coef(rFunc);
  FunctionCoefficient r_inv_coef(rinvFunc);

  LambdaCoefficient w_mu_r_sigma_coef(
      [&sigma_gf_coef](mfem::ElementTransformation& T,
                       const mfem::IntegrationPoint& ip) {
        std::array<double, 3> x;
        Vector transip(x.data(), 3);
        T.Transform(ip, transip);
        return w * mu * transip(1) * sigma_gf_coef.Eval(T, ip);
      });

  // Projecting the zeroth-order L2 element data on 1st order H1 space
  ParGridFunction w_mu_r_sigma_gf(&HGradFESpace);
  GridFunctionCoefficient w_mu_r_sigma_gf_coef(&w_mu_r_sigma_gf);

  ParSesquilinearForm* a1 = new ParSesquilinearForm(&HGradFESpace, conv);
  a1->AddDomainIntegrator(new DiffusionIntegrator(r_coef), nullptr);
  a1->AddDomainIntegrator(new MassIntegrator(r_inv_coef), nullptr);
  a1->AddDomainIntegrator(nullptr, new MassIntegrator(w_mu_r_sigma_gf_coef));

  // Set up the parallel bilinear form for the preconditioner
  // corresponding to the operator :
  //<r grad u, grad v> + <1/r u, v> + <w mu r sigma u, v>
  ParBilinearForm* pcOp = new ParBilinearForm(&HGradFESpace);
  pcOp->AddDomainIntegrator(new DiffusionIntegrator(r_coef));
  pcOp->AddDomainIntegrator(new MassIntegrator(r_inv_coef));
  pcOp->AddDomainIntegrator(new MassIntegrator(w_mu_r_sigma_gf_coef));

  Array<int> E_ess_tdof_list;
  

  HGradFESpace.GetEssentialTrueDofs(E_ess_bdr_zero, E_ess_tdof_list);

  OperatorHandle A1;
  Vector B1, X1;

  // ===========================================================
  // End of initialization for E Computation
  // ===========================================================

  // ===========================================================
  // Initialization for Joule Heating
  // ===========================================================
  ParGridFunction W(&HGradFESpace);
  JouleHeatingCoefficient w_coeff(sigma_gf_coef, E);
  FunctionCoefficient two_pi_r_coef(two_pi_rFunc);
  TotalJouleHeatingCoefficient two_pi_r_w_coeff(two_pi_r_coef, sigma_gf_coef,
                                                E);
  ParLinearForm W_total(&HGradFESpace);
  W_total.AddDomainIntegrator(new DomainLFIntegrator(two_pi_r_w_coeff));

  FunctionCoefficient one(unit_coeff);
  ParGridFunction ones(&HGradFESpace);
  ones.ProjectCoefficient(one);

  // ===========================================================
  // End of initialization for Joule Heating
  // ===========================================================

  // ===========================================================
  // Initialization for B term
  // ===========================================================
  HYPRE_Int glob_size_nd = HCurlFESpace.GlobalTrueVSize();
  // grid function for B
  ParGridFunction B_re(&HCurlFESpace);
  ParGridFunction B_im(&HCurlFESpace);

  ParGridFunction grad_E_re_gf(&HCurlFESpace);
  ParGridFunction grad_E_im_gf(&HCurlFESpace);
  ParGridFunction rEre_gf(&HGradFESpace);
  ParGridFunction rEim_gf(&HGradFESpace);
  ParGridFunction grad_rEre_gf(&HCurlFESpace);
  ParGridFunction grad_rEim_gf(&HCurlFESpace);
  ParDiscreteLinearOperator* grad =
      new ParDiscreteLinearOperator(&HGradFESpace, &HCurlFESpace);
  grad->AddDomainInterpolator(new GradientInterpolator());
  grad->Assemble();

  rEreCoefficient r_E_re_coef(r_coef, E);
  rEimCoefficient r_E_im_coef(r_coef, E);

  BreCoefficient B_re_coef(dim, r_inv_coef, grad_E_im_gf, grad_rEim_gf);
  BimCoefficient B_im_coef(dim, r_inv_coef, grad_E_re_gf, grad_rEre_gf);

  // ===========================================================
  // End of initialization for B term
  // ===========================================================

  // ===========================================================
  // Initialization for Average Lorentz force
  // ===========================================================
  ParGridFunction Fz(&HGradFESpace);
  AxialLorentzForceCoefficient Fz_coeff(sigma_gf_coef, E, B_re, B_im);

  ParGridFunction Fr(&HGradFESpace);
  RadialLorentzForceCoefficient Fr_coeff(sigma_gf_coef, E, B_re, B_im);
  // ===========================================================
  // End of initialization for Average Lorentz force
  // ===========================================================

  // ===========================================================
  // Initialization of IO
  // ===========================================================
  ParaViewDataCollection paraview_dc(basename, pmesh);
  paraview_dc.SetPrefixPath("ParaView");
  paraview_dc.SetLevelsOfDetail(order);
  paraview_dc.SetDataFormat(VTKFormat::BINARY);
  paraview_dc.SetHighOrderOutput(true);
  paraview_dc.RegisterField("E_re", &E.real());
  paraview_dc.RegisterField("E_im", &E.imag());
  paraview_dc.RegisterField("B_re", &B_re);
  paraview_dc.RegisterField("B_im", &B_im);
  paraview_dc.RegisterField("W", &W);
  paraview_dc.RegisterField("Fx", &Fz);
  paraview_dc.RegisterField("Fy", &Fr);
  paraview_dc.RegisterField("sigmaL2", &sigma_L2_order0_gf);
  paraview_dc.RegisterField("sigmaHGrad", &sigma_gf);
  // ===========================================================
  // End of initialization of IO
  // ===========================================================

  int precice_iter = 0;  // Used to track number of times called when coupled
  if (precice_active) {
    precice_dt = precice_adapter->Initialize();
  }

  // We need to loop to readily handle successive couplings
  // from Hegel.
  do {
    // If precice active, update quantity for sigma
    if (precice_active) {
      VolumeCouplingCoefficient* coef =
          dynamic_cast<VolumeCouplingCoefficient*>(sigma_coef);
      coef->UpdateData();  // Performs a preCICE ReadBlockScalar
    }
    // Update value for sigma gf and things that depend on it.
    sigma_L2_order0_gf.ProjectCoefficient(*sigma_coef);
    SuccessiveRefinementInterpolation(sigma_L2_order0_gf, sigma_gf, working_gf);
    w_mu_r_sigma_gf.ProjectCoefficient(w_mu_r_sigma_coef);

    //-------------------------------------
    //---------E computation starts--------
    //------------------------------------
    b.Assemble();
 
    a1->Update();
    a1->Assemble();
    a1->Finalize();

    pcOp->Update();
    pcOp->Assemble();
    pcOp->Finalize();

    
    // Inside since sigma might change
    a1->FormLinearSystem(E_ess_tdof_list, E, b, A1, X1, B1, 1);

    //     Define and apply a parallel FGMRES solver for A1 X1 = B1 with a block
    //     diagonal preconditioner based on the appropriate multigrid
    //     preconditioner from hypre
    {
      Array<int> blockTrueOffsets;
      blockTrueOffsets.SetSize(3);
      blockTrueOffsets[0] = 0;
      blockTrueOffsets[1] = A1->Height() / 2;
      blockTrueOffsets[2] = A1->Height() / 2;
      blockTrueOffsets.PartialSum();

      BlockDiagonalPreconditioner BDP(blockTrueOffsets);

      OperatorHandle PCOp;
      pcOp->FormSystemMatrix(E_ess_tdof_list, PCOp);

      auto amg = new HypreBoomerAMG(*PCOp.As<HypreParMatrix>());
      amg->SetPrintLevel(SOLVER_PRINT_LEVEL);

      Operator* pc_r = amg;
      
      //Scale by -1 for Hermitian case, 1 for block symmetric case
      //Note : For highly stretched grids, -1 scaling does not converge
      //For relatively uniform grids, -1 scaling was working fine
      Operator* pc_i = new ScaledOperator(pc_r, 1.0); 
 
      BDP.SetDiagonalBlock(0, pc_r);
      BDP.SetDiagonalBlock(1, pc_i);
      BDP.owns_blocks = 1;


      FGMRESSolver fgmres(MPI_COMM_WORLD);
      fgmres.iterative_mode = true;
      fgmres.SetPreconditioner(BDP);
      fgmres.SetOperator(*A1.Ptr());
      fgmres.SetAbsTol(SOLVER_TOL);
      fgmres.SetMaxIter(SOLVER_MAX_IT);
      fgmres.SetPrintLevel(SOLVER_PRINT_LEVEL);
      fgmres.Mult(B1, X1);

    }

    a1->RecoverFEMSolution(X1, b, E);

    //-------------------------------------------------
    //----------E computation ends--------------------
    //-------------------------------------------------

    //--------------------------------------------------------
    //-------------Computation of joule heating---------------
    //--------------------------------------------------------

    // Reproject since depends on sigma
    W.ProjectCoefficient(w_coeff);
   
    W_total.Update();
    W_total.Assemble();

    double P_total = W_total(ones);

    if (mpi.Root()) {
      cout << "Initial Total power =    " << P_total / 1000 << "   KW" << endl;
    }

    const double gamma = sqrt(P_target / P_total);

    if (mpi.Root()) {
      cout << "Scaling factor =   " << gamma << endl;
    }

    // Scaling the electric field
    E.real() *= gamma;
    E.imag() *= gamma;

    // Scaling the source current
    for (int i = 0; i < nb_coils; i++) {
      double s_scaled = gamma*delta_current[i]->Scale();
      delta_current[i]->SetScale(s_scaled); 
    }
  
    // Computing the final power again with the scaled Electric field
    // Reassembling will use newly updated E field
    W.ProjectCoefficient(w_coeff);
    W_total.Update();
    W_total.Assemble();

    P_total = W_total(ones);

    if (mpi.Root()) {
      cout << "Final Total power =   " << P_total / 1000 << "   KW" << endl;
      }

    //-------------------------------------------------
    //----------B computation starts-------------------
    //-------------------------------------------------

    if (mpi.Root()) {
      cout << "Number of Magnetic Field unknowns:    " << glob_size_nd << endl;
    }

    grad->Mult(E.real(), grad_E_re_gf);
    grad->Mult(E.imag(), grad_E_im_gf);

    rEre_gf.ProjectCoefficient(r_E_re_coef);
    rEim_gf.ProjectCoefficient(r_E_im_coef);

    grad->Mult(rEre_gf, grad_rEre_gf);
    grad->Mult(rEim_gf, grad_rEim_gf);

    B_re.ProjectCoefficient(B_re_coef);
    B_im.ProjectCoefficient(B_im_coef);

    B_re *= (1.0 / w);
    B_im *= (1.0 / w);

    //---------------------------------------------------
    //------------B Computation ends---------------------
    //---------------------------------------------------

    //---------------------------------------------------
    //-----------Computation of Average Lorentz force----
    //---------------------------------------------------

    Fz.ProjectCoefficient(Fz_coeff);
    Fr.ProjectCoefficient(Fr_coeff);

    //------------------------------------------------------------
    //----------Computation of Lorentz force ends-----------------
    //------------------------------------------------------------

    //----------Paraview visualization--------------------------
    
    if (precice_iter%print_interval == 0)
      {
	paraview_dc.SetCycle(precice_iter);
	paraview_dc.SetTime(static_cast<double>(precice_iter*precice_dt));
	paraview_dc.Save();
      }
    ++precice_iter;
    
    if (precice_active) {
      // Write fields.

      // Write Joule Heating
      tmp_gf_L20.ProjectGridFunction(W);      
      write_data.SetSize(static_cast<int>(precice_to_mfem_map.size()));
      for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
        write_data(i) = tmp_gf_L20(precice_to_mfem_map[i]);
      }
      const double* data_ptr1 = static_cast<const double*>(write_data);
      precice_adapter->WriteBlockScalarData("JouleHeating", data_ptr1);

      // Write Axial Lorentz Force
      tmp_gf_L20.ProjectGridFunction(Fz);
      for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
        write_data(i) = tmp_gf_L20(precice_to_mfem_map[i]);	
      }
      const double* data_ptr2 = static_cast<const double*>(write_data);
      precice_adapter->WriteBlockScalarData("AxialLorentzForce", data_ptr2);

   
      // Write Radial Lorentz Force
      tmp_gf_L20.ProjectGridFunction(Fr);
      for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
        write_data(i) = tmp_gf_L20(precice_to_mfem_map[i]);
      }
      const double* data_ptr3 = static_cast<const double*>(write_data);
      precice_adapter->WriteBlockScalarData("RadialLorentzForce", data_ptr3);

    
      // Write Lorentz Force vector
      /*write_data.SetSize(2 * static_cast<int>(precice_to_mfem_map.size()));
      tmp_gf_L20.ProjectGridFunction(Fz);
      for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
        write_data(2 * i) = tmp_gf_L20(precice_to_mfem_map[i]);
      }
      tmp_gf_L20.ProjectGridFunction(Fr);
      for (int i = 0; i < precice_to_mfem_map.size(); ++i) {
        write_data(2 * i + 1) = tmp_gf_L20(precice_to_mfem_map[i]);
      }

      data_ptr = static_cast<const double*>(write_data);
      precice_adapter->WriteBlockVectorData("LorentzForce", data_ptr);*/

      precice_adapter->Advance(precice_dt);
  
    }
    
   } while ((precice_active && precice_adapter->IsCouplingOngoing()));

  if (precice_active) {
    precice_adapter->Finalize();
  }

  // Free memory
  delete b_re;
  delete b_im;
  delete a1;
  delete pcOp;
  delete grad;
  delete sigma_coef;
  delete precice_adapter;

  for (int s = 0; s < H1_collections.size(); ++s) {
    delete H1_collections[s];
    delete H1_fes[s];
  }
  H1_collections.clear();
  H1_fes.clear();

  delete pmesh;

  return 0;
}

double current_source(const Vector& x) {
  double Js;
  double r = x[1];
  double z = x[0];
  double J0 = -100.0;

  //Plasmatron_zero coil with rectangular c/s
  if (r >= 12.7e-3 && r <= 19.0e-3)
    {
      if (z >= 116.2e-3 && z <= 132.2e-3)
	{
	  Js = J0;
	}
      else
	{
	  Js = 0.0;
	}
    }
  else
    {
      Js = 0.0;
    }


  
  // Circular coil source with radius of coil = R and
  // radius of c/s of coil ri located at x = xi
  /*double R = 0.109;
  double ri = 0.01;
  double x1 = 0.127;
  double x2 = 0.177;
  double x3 = 0.227;
  double x4 = 0.277;
  double x5 = 0.327;
  double x6 = 0.377;

  double J0 = -1.0e11;
  double r = x[1];
  double Js;


  if (r <= (R+ri) && r >= (R-ri))
    {
      if (x[0] <= (x1 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x1 - sqrt(ri*ri -
  (r-R)*(r-R))))
   {
     Js = J0;
   }
      else if (x[0] <= (x2 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x2 -
  sqrt(ri*ri - (r-R)*(r-R))))
   {
     Js = J0;
   }
      else if (x[0] <= (x3 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x3 -
  sqrt(ri*ri - (r-R)*(r-R))))
   {
     Js = J0;
   }
      else if (x[0] <= (x4 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x4 -
  sqrt(ri*ri - (r-R)*(r-R))))
   {
     Js = J0;
   }
      else if (x[0] <= (x5 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x5 -
  sqrt(ri*ri - (r-R)*(r-R))))
   {
     Js = J0;
   }
      else if (x[0] <= (x6 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x6 -
  sqrt(ri*ri - (r-R)*(r-R))))
   {
     Js = J0;
   }
      else
   {
     Js = 0.0;
   }

    }
  else
    {
      Js = 0.0;

      }*/

  // Circular coil source with rectangular c/s
  /*double r1 = 0.4;
  double r2 = 0.5;
  double x1 = -0.1;
  double x2 = 0.1;
  double J0 = -10000.0;
  double r = x[1];
  double Js;

  if (r <= r2 && r >= r1)
    {
      if (x[0] <= x2 && x[0] >= x1)
   {
     Js = J0;
   }
      else
   {
     Js = 0.0;
   }
    }
  else
    {
      Js = 0.0;
      }*/

  // No current source
  //Js = 0.0;
  
  return w * mu * r * Js;
}


double voltage_bc(const Vector& x) { return 200.0; }

double unit_coeff(const Vector& x) {
  double r = x[1];
  return (r <= cr && (x[0] >= cls && x[0] <= cle )) ? 1.0 : 0.0;
}

double sigma_conductor(const Vector& x) { return sigmaCoil; }

double sigma_non_conductor(const Vector& x) { return sigmaAir; }

// Rewrote this as a "LambdaCoefficient" taken
// from CHyPS. Can then capture sigma and
// reuse the same one used in the other equation.
// double w_mu_r_sigmaFunc(const Vector& x) {
//   double r = x[1];
//   double z = x[0];

//   return w * mu * r * sig;
// }

double sigmaFunc(const Vector& x) {
  double r = x[1];
  double z = x[0];

  return (r <= cr && (x[0] >= cls && x[0] <= cle )) ? 100.0 : sigmaAir;
  
}

double rFunc(const Vector& x) { return x[1]; }

double rinvFunc(const Vector& x) { return 1.0 / x[1]; }

double two_pi_rFunc(const Vector& x) { return 2.0 * M_PI * x[1]; }


std::vector<int> BuildPreciceToMfemMap(ParMesh* mesh, const double x1, const double x2, const double r) {
  std::vector<int> precice_to_mfem_map;
  mfem::Vector center;
  for (int i = 0; i < mesh->GetNE(); ++i) {
    mesh->GetElementCenter(i, center);
    if (center(1) <= r) {
      if ((center(0) >= x1) && (center(0) <= x2)) {
         precice_to_mfem_map.push_back(i);
      }
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

void SuccessiveRefinementInterpolation(
    const ParGridFunction& starting_l2, ParGridFunction& final_h1,
    std::vector<ParGridFunction>& working_gf) {
  // Project L2 to order 1 H1
  GridFunctionCoefficient starting_gfc(&starting_l2);
  working_gf[0].ProjectDiscCoefficient(starting_gfc,
                                       GridFunction::AvgType::ARITHMETIC);

  const int total_levels = static_cast<int>(working_gf.size());
  // Assumes whole mesh is of same order
  const int order = final_h1.FESpace()->GetOrder(0);

  DEBUG_ASSERT(
      order <= total_levels, global_assert{}, DebugLevel::CHEAP{},
      "working_gf vector is too small. Does not contain adequate levels of H1 "
      "refinement");
  if (order == 1) {
    final_h1 = working_gf[0];
  } else {
    for (int s = 1; s < order - 1; ++s) {
      GridFunctionCoefficient h1_gfc(&working_gf[s - 1]);
      working_gf[s].ProjectCoefficient(h1_gfc);
    }
    GridFunctionCoefficient h1_gfc(&working_gf[order - 2]);
    final_h1.ProjectCoefficient(h1_gfc);
  }
}
