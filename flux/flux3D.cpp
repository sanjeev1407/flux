// This file is part of the Finite-element soLver for Unsteady electromagnetiX (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=====================================================================================================

#include "flux3D_solver.hpp"
#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace mfem;
using namespace mfem::common;
using namespace mfem::electromagnetics;

  
//Physical constants
double w = 2.0*3.14*0.37e6; //Angular frequency (rad/s)
double sigma = 6.0e7; //Electrical conductivity
double sigmaAir = 1.0e-3;
double mu = 1.2566e-6; //Magnetic Permeability
double P_target = 150000.0; //Desired total power(Watts)
int dim;

/*Array<double> xc(501);
Array<double> yc(101);
double sigma_dist[101][501];
double xc_min, yc_min, dx, dy;*/

Array<double> xc(97);
Array<double> yc(100);
double sigma_dist[100][97];
double xc_min, yc_min, dx, dy;

//electron properties array
Array<double> e_xc(97);
Array<double> e_yc(100);
double Pe[100][97];
double ne[100][97];
double e_xc_min, e_yc_min, e_dx, e_dy;

int main(int argc, char *argv[])
{
   //Initialize MPI
   MPI_Session mpi(argc, argv);
   int myid = mpi.WorldRank();

   //Parse command-line options.
   const char *mesh_file = "beam-tet.mesh";
   int order = 2;
   const char *basename= "emsolve";
   
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree).");
   args.Parse();

   //Reading and storing sigma data from a file
   /*int nb_x, nb_y;   

   ifstream infile;
   infile.open("sigma_distribution.dat");
   infile >> nb_x >> nb_y;
   
   for (int j = 0; j < nb_y; j++)
     {
    for (int i = 0; i <nb_x; i++)
      {
	infile >> xc[i] >> yc[j] >> sigma_dist[j][i];
      }
     }

   infile.close();

   xc_min = xc[0];
   yc_min = yc[0];
   dx = xc[1] - xc[0];
   dy = yc[1] - yc[0];*/

   int nb_x, nb_y;   

   ifstream infile;
   infile.open("sigma_coolfluid_uniform_grid.dat");
   infile >> nb_x >> nb_y;
   
   for (int i = 0; i < nb_x; i++)
     {
    for (int j = 0; j < nb_y; j++)
      {
	infile >> xc[i] >> yc[j] >> sigma_dist[j][i];
      }
     }

   infile.close();

   xc_min = xc[0];
   yc_min = yc[0];
   dx = xc[1] - xc[0];
   dy = yc[1] - yc[0];

   //Reading electron properties file
   int e_nb_x, e_nb_y;   

   ifstream infile2;
   infile2.open("electron_properties_coolfluid_uniform_grid.dat");
   infile2 >> e_nb_x >> e_nb_y;
   
   for (int i = 0; i < e_nb_x; i++)
     {
    for (int j = 0; j < e_nb_y; j++)
      {
	infile2 >> e_xc[i] >> e_yc[j] >> Pe[j][i] >> ne[j][i];
      }
     }

   infile2.close();

   e_xc_min = e_xc[0];
   e_yc_min = e_yc[0];
   e_dx = e_xc[1] - e_xc[0];
   e_dy = e_yc[1] - e_yc[0];

   // 3. Here material properties are assigned to mesh attributes.
   std::map<int, std::function<double(const Vector &)>> sigmaMap;

   sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(1, sigma_profile));
   sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(2, sigma_conductor));
   sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(3, sigma_conductor));
   sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(4, sigma_conductor));
   sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(5, sigma_conductor));
   sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(6, sigma_conductor));
   sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(7, sigma_conductor));
   sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(8, sigma_non_conductor));
   /*sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(1, sigma_conductor));
     sigmaMap.insert(pair<int, std::function<double(const Vector &)>>(2, sigma_non_conductor));*/
   
   MeshDependentCoefficient *sigma_coefficient  = new MeshDependentCoefficient(sigmaMap);
   
   
   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   dim = mesh->Dimension();
   int sdim = mesh->SpaceDimension();
   
   // 8. Define a parallel mesh by a partitioning of the serial mesh.
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);

   //Sometimes it is needed with meshes having tetrahedrals
   pmesh->ReorientTetMesh();
   
   //Defining all the Finite element collections and spaces to be used
   H1_FECollection HGradFEC(order, dim);
   ParFiniteElementSpace  HGradFESpace(pmesh, &HGradFEC);
   
   ND_FECollection HCurlFEC(order, dim);
   ParFiniteElementSpace HCurlFESpace(pmesh, &HCurlFEC);

   RT_FECollection HDivFEC(order, dim);
   ParFiniteElementSpace  HDivFESpace(pmesh, &HDivFEC);

   L2_FECollection L2FEC(order, dim);
   ParFiniteElementSpace  L2FESpace(pmesh, &L2FEC);
    

   //------------------------------------------------
   //------------Solving for phi(Optional)-----------
   //------------------------------------------------
   
   // Assign the essential boundary conditions
   /*Array<int> phi_ess_bdr_zero(mesh->bdr_attributes.Max());
   Array<int> phi_ess_bdr_non_zero(mesh->bdr_attributes.Max());
   
   // For a solid cylinder problem we have 3 surfaces
   // 1)Z=0 plane, 2)Z=Lplane, 3)circular surface 
   //Electrostatic potential BC 
   phi_ess_bdr_zero = 0;
   phi_ess_bdr_zero[1] = 1;
   
   phi_ess_bdr_non_zero = 0;
   phi_ess_bdr_non_zero[0] = 1;
   
   
   // The terminology is TrueVSize is the unique (non-redundant) number of dofs
   HYPRE_Int glob_size_h1 = HGradFESpace.GlobalTrueVSize();

   if (mpi.Root())
     {
       cout << "Number of Electrostatic potential unknowns:     " << glob_size_h1 << endl;
     }

   // grid function for phi
   ParGridFunction phi(&HGradFESpace);

   // Electrostatic potential bc on surface
   phi = 0.0;
   
   FunctionCoefficient voltage(voltage_bc);
   phi.ProjectBdrCoefficient(voltage,phi_ess_bdr_non_zero);
   
   ConstantCoefficient zero_voltage(0.0);
   phi.ProjectBdrCoefficient(zero_voltage,phi_ess_bdr_zero);
      
   Array<int> phi_ess_tdof_list;
   Array<int> phi_ess_bdr(mesh->bdr_attributes.Max());
   
   for (int i=0; i<mesh->bdr_attributes.Max(); i++)
     {
       phi_ess_bdr[i] = phi_ess_bdr_zero[i] + phi_ess_bdr_non_zero[i];
     }
   
   //Eliminate the dofs corresponding to the boundary faces
   HGradFESpace.GetEssentialTrueDofs(phi_ess_bdr, phi_ess_tdof_list);
   
   //Form linear system for potential : M(sigma) laplacian phi = 0
   // RHS = 0
   ParGridFunction *v0 = new ParGridFunction(&HGradFESpace);
   *v0 = 0.0;
   
   //Laplacian bilinear operator
   ParBilinearForm *a0 = new ParBilinearForm(&HGradFESpace);
   a0->AddDomainIntegrator(new DiffusionIntegrator(*sigma_coefficient));
   a0->Assemble();
   
   HypreParMatrix *A0 = new HypreParMatrix;
   Vector *X0 = new Vector;
   Vector *B0 = new Vector;
   
   a0->FormLinearSystem(phi_ess_tdof_list,phi,*v0,*A0,*X0,*B0);
   
   //Definition of linear solver and preconditioner
   HypreSolver * amg_a0;
   HyprePCG    * pcg_a0;
   
   amg_a0 = new HypreBoomerAMG(*A0);
   pcg_a0 = new HyprePCG(*A0);
   pcg_a0->SetTol(SOLVER_TOL);
   pcg_a0->SetMaxIter(SOLVER_MAX_IT);
   pcg_a0->SetPrintLevel(SOLVER_PRINT_LEVEL);
   pcg_a0->SetPreconditioner(*amg_a0);
     
   // pcg "Mult" operation is a solve
   // X0 = A0^(-1) * B0
   pcg_a0->Mult(*B0, *X0);

   a0->RecoverFEMSolution(*X0,*v0,phi);*/

   //--------------------------------------
   //---------phi computation ends---------
   //--------------------------------------


   //-------------------------------------
   //---------E computation starts--------
   //------------------------------------

   ComplexOperator::Convention conv = ComplexOperator::HERMITIAN;
   
   // Assign the boundary conditions
   Array<int> E_ess_bdr_zero(mesh->bdr_attributes.Max());

   E_ess_bdr_zero = 1;
   
   HYPRE_Int glob_size_nd = HCurlFESpace.GlobalTrueVSize();

   if (mpi.Root())
     {
       cout << "Number of Electric Field unknowns:    " << glob_size_nd << endl;
     }

    //    Set up the parallel complex linear form b(.) which corresponds to the
    //    right-hand side (- iJc) of the FEM linear system.
    VectorFunctionCoefficient source(3, current_source);
    Vector zero_vec(3); zero_vec = 0.0;
    VectorConstantCoefficient Zero_vec(zero_vec);
    
    ParLinearForm *b_re = new ParLinearForm(&HCurlFESpace);
    ParLinearForm *b_im = new ParLinearForm(&HCurlFESpace);
 
    ParComplexLinearForm b(&HCurlFESpace, b_re, b_im, conv);
    b.AddDomainIntegrator(new VectorFEDomainLFIntegrator(Zero_vec), new VectorFEDomainLFIntegrator(source));
    b.Assemble();
    
    
    //If we want to apply voltage source
    /*ParBilinearForm *m1 = new ParBilinearForm(&HCurlFESpace);
    m1->AddDomainIntegrator(new VectorFEMassIntegrator(*sigma_coefficient));
    m1->Assemble();

    ParDiscreteLinearOperator *grad = new ParDiscreteLinearOperator(&HGradFESpace, &HCurlFESpace);
    grad->AddDomainInterpolator(new GradientInterpolator());
    grad->Assemble();

    //Complex grid function for RHS (i M1(sigma) grad phi)
    ParComplexGridFunction v1(&HCurlFESpace);
    v1.real() = 0.0;
    v1.imag() = 0.0;

    ParGridFunction grad_phi(&HCurlFESpace);

    grad->Mult(phi,grad_phi);
    m1->AddMult(grad_phi,v1.imag(),1.0);

    b.Add(1.0, v1);*/
    
    //     Define the solution vector E as a parallel complex finite element grid
    //     function corresponding to fespace
    ParComplexGridFunction E(&HCurlFESpace);

    //Initialize E to 0 + i0
    E.ProjectCoefficient(Zero_vec, Zero_vec);

    //n x E = 0 Essential BC
    E.real().ProjectBdrCoefficientTangent(Zero_vec,E_ess_bdr_zero);
    E.imag().ProjectBdrCoefficientTangent(Zero_vec,E_ess_bdr_zero);
    
    // Set up the parallel sesquilinear form a(.,.) on the finite element
    // space corrsponding to Curl (1/(w*mu) Curl) + i sigma
    double inv_wmu = 1.0/(w*mu);
    ConstantCoefficient stiffnessCoef(inv_wmu);

    ParSesquilinearForm *a1 = new ParSesquilinearForm(&HCurlFESpace, conv);
    a1->AddDomainIntegrator(new CurlCurlIntegrator(stiffnessCoef),
			   NULL);
    a1->AddDomainIntegrator(NULL,
			   new VectorFEMassIntegrator(*sigma_coefficient));
    a1->Assemble();
     
    // Set up the parallel bilinear form for the preconditioner
    // corresponding to the operator : Curl(1/(w*mu) Curl) + sigma
    ParBilinearForm *pcOp = new ParBilinearForm(&HCurlFESpace);
    pcOp->AddDomainIntegrator(new CurlCurlIntegrator(stiffnessCoef));
    pcOp->AddDomainIntegrator(new VectorFEMassIntegrator(*sigma_coefficient));
    pcOp->Assemble();


    Array<int> E_ess_tdof_list;
    HCurlFESpace.GetEssentialTrueDofs(E_ess_bdr_zero, E_ess_tdof_list);
    
    OperatorHandle A1;
    Vector B1, X1;

    a1->FormLinearSystem(E_ess_tdof_list, E, b, A1, X1, B1);

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

      Operator * pc_r = NULL;
      Operator * pc_i = NULL;

      OperatorHandle PCOp;
      pcOp->FormSystemMatrix(E_ess_tdof_list, PCOp);

      pc_r = new HypreAMS(*PCOp.As<HypreParMatrix>(), &HCurlFESpace);

      pc_i = new ScaledOperator(pc_r, -1.0);

      BDP.SetDiagonalBlock(0, pc_r);
      BDP.SetDiagonalBlock(1, pc_i);
      BDP.owns_blocks = 1;

      FGMRESSolver fgmres(MPI_COMM_WORLD);
      fgmres.SetPreconditioner(BDP);
      fgmres.SetOperator(*A1.Ptr());
      fgmres.SetRelTol(SOLVER_TOL);
      fgmres.SetMaxIter(SOLVER_MAX_IT);
      fgmres.SetPrintLevel(SOLVER_PRINT_LEVEL);
      fgmres.Mult(B1, X1);
    }
    
    a1->RecoverFEMSolution(X1, b, E);

    // If we are using a potential source
    // the total field is E_tot = E_ind - Grad Phi
    // so we need to subtract out Grad Phi
    // E = E - grad (P)
    //grad->AddMult(phi,E.real(),-1.0);

    //-------------------------------------------------
    //----------E computation ends--------------------
    //-------------------------------------------------

    
    //---------------------------------------------------
    //------------sigma grid function--------------------
    //---------------------------------------------------
    /*ParGridFunction conductivity(&HGradFESpace);
    conductivity.ProjectCoefficient(*sigma_coefficient);*/

    
    //------------------------------------------------------
    //------------Es = gradp Pe/ne*qe-----------------------
    //------------------------------------------------------
    /*ParDiscreteLinearOperator *grad = new ParDiscreteLinearOperator(&HGradFESpace, &HCurlFESpace);
    grad->AddDomainInterpolator(new GradientInterpolator());
    grad->Assemble();

    // grid function for Pe
    ParGridFunction Pe_gf(&HGradFESpace);
    
    FunctionCoefficient P_electron(Pe_func);
    Pe_gf.ProjectCoefficient(P_electron);

    ParGridFunction grad_Pe(&HCurlFESpace);
    grad->Mult(Pe_gf, grad_Pe);

    double qe = 1.60217662e-19;
    double inv_qe;
    inv_qe = 1.0/qe;
    grad_Pe *= inv_qe;

    // grid function for ne
    ParGridFunction ne(&HGradFESpace);

    FunctionCoefficient n_electron(ne_func);
    ne.ProjectCoefficient(n_electron);

    ParGridFunction E_Ambipolar(&HCurlFESpace);

    AmbipolarCoefficient Es_coeff(3, grad_Pe, ne);

    E_Ambipolar.ProjectCoefficient(Es_coeff);

    E.real().Add(1.0, E_Ambipolar);*/

    //--------------------------------------------------------
    //-------------Computation of joule heating---------------
    //--------------------------------------------------------
    
    ParGridFunction W(&HGradFESpace);

    JouleHeatingCoefficient w_coeff(*sigma_coefficient, E);

    W.ProjectCoefficient(w_coeff);
    
    ParLinearForm W_total(&HGradFESpace);
    W_total.AddDomainIntegrator(new DomainLFIntegrator(w_coeff));
    W_total.Assemble();

    FunctionCoefficient one(unit_coeff);
    ParGridFunction ones(&HGradFESpace);
    ones.ProjectCoefficient(one);

    double P_total = W_total(ones);

    if (mpi.Root())
      {
	cout << "Initial Total power =    " << P_total/1000 << "   KW" << endl;
      }
    

    
    double gamma;
    gamma = sqrt(P_target/P_total);

    if (mpi.Root())
      {
	cout << "Scaling factor =   " << gamma << endl;
      }
    
    
    
    //Scaling the electric field
    E.real() *= gamma;
    E.imag() *= gamma;

    //Computing the final power again with the scaled Electric field
    ParGridFunction W_final(&HGradFESpace);

    JouleHeatingCoefficient w_coeff_final(*sigma_coefficient, E);

    W_final.ProjectCoefficient(w_coeff_final);
    
    ParLinearForm W_total_final(&HGradFESpace);
    W_total_final.AddDomainIntegrator(new DomainLFIntegrator(w_coeff_final));
    W_total_final.Assemble();

    double P_total_final = W_total_final(ones);

    if (mpi.Root())
      {
	cout << "Final Total power =   " << P_total_final/1000 << "   KW" << endl;

	}
    
    
    
    //-------------------------------------------------
    //----------B computation starts-------------------
    //-------------------------------------------------

    HYPRE_Int glob_size_rt = HDivFESpace.GlobalTrueVSize();
    if (mpi.Root())
      {
	cout << "Number of Magnetic Field unknowns:    " << glob_size_rt << endl;
      }
    
    // grid function for B
    ParGridFunction B_re(&HDivFESpace);
    ParGridFunction B_im(&HDivFESpace);
    
    ParDiscreteLinearOperator *curl = new ParDiscreteLinearOperator(&HCurlFESpace, &HDivFESpace);
    curl->AddDomainInterpolator(new CurlInterpolator);
    curl->Assemble();

    // Compute B_re = - (1/w) Curl(E_im)
    curl->Mult(E.imag(), B_re);
    double inv_w = 1.0/w;
    B_re *= -inv_w;

    // Compute B_im =  (1/w) Curl(E_re)
    curl->Mult(E.real(), B_im);
    B_im *= inv_w;


    //---------------------------------------------------
    //------------B Computation ends---------------------
    //---------------------------------------------------

    //================================================
    //=========Checking Div.(sigma*E) = 0 ============
    //================================================
    ParGridFunction sigmaEre_gf(&HDivFESpace);
    sigmaEreCoefficient sig_Ere(3,*sigma_coefficient, E);
    sigmaEre_gf.ProjectCoefficient(sig_Ere);
    
    ParGridFunction sigmaEim_gf(&HDivFESpace);
    sigmaEimCoefficient sig_Eim(3,*sigma_coefficient, E);
    sigmaEim_gf.ProjectCoefficient(sig_Eim);
    
    
    ParGridFunction div_Ere(&L2FESpace);
    ParGridFunction div_Eim(&L2FESpace);
    
    ParMixedBilinearForm *weakDiv = new ParMixedBilinearForm(&HDivFESpace, &L2FESpace);
    weakDiv->AddDomainIntegrator(new VectorFEDivergenceIntegrator());
    weakDiv->Assemble();

    sigmaEre_gf.Add(1.0, b.imag());
    
    weakDiv->Mult(sigmaEre_gf, div_Ere);
    weakDiv->Mult(sigmaEim_gf, div_Eim);

    //================================================
    //==========End of divergence check===============
    //================================================
    
    //---------------------------------------------------
    //-----------Computation of Lorentz force-----------
    //--------------------------------------------------
    /*ParFiniteElementSpace  trial_fes(pmesh, &HDivFEC);
    ParFiniteElementSpace  test_fes(pmesh, &HDivFEC);

    ParGridFunction x(&test_fes);

    HYPRE_Int trial_size = trial_fes.GlobalTrueVSize();
    HYPRE_Int test_size = test_fes.GlobalTrueVSize();
    
    ConstantCoefficient unity(1.0);
    ParBilinearForm a(&test_fes);
    ParMixedBilinearForm a_mixed(&trial_fes, &test_fes);

    //Need to compute the vector coefficient sigma*E_re
    sigmaEreCoefficient sigma_E_re(3,*sigma_coefficient, E);
    
    a.AddDomainIntegrator(new VectorFEMassIntegrator(unity));
    a_mixed.AddDomainIntegrator(new MixedCrossProductIntegrator(sigma_E_re));
    
    a.Assemble();
    a.Finalize();
    
    a_mixed.Assemble();
    a_mixed.Finalize();
    
    Vector B(test_fes.GetTrueVSize());
    Vector X(test_fes.GetTrueVSize());
    
    HypreParMatrix *mixed = a_mixed.ParallelAssemble();
    
    Vector P(trial_fes.GetTrueVSize());
    B_re.GetTrueDofs(P);

    mixed->Mult(P,B);
    
    delete mixed;
    

    //Define and apply a parallel PCG solver for AX=B with Jacobi preconditioner
    HypreParMatrix *Amat = a.ParallelAssemble();
    HypreDiagScale Jacobi(*Amat);
    HyprePCG pcg(*Amat);
    pcg.SetTol(SOLVER_TOL);
    pcg.SetMaxIter(SOLVER_MAX_IT);
    pcg.SetPrintLevel(SOLVER_PRINT_LEVEL);
    pcg.SetPreconditioner(Jacobi);
    X = 0.0;
    pcg.Mult(B, X);

    delete Amat;

    x.SetFromTrueDofs(X);*/


    //Average Lorentz force
    /*ParGridFunction F(&HCurlFESpace);

    LorentzForceCoefficient F_coeff(3,*sigma_coefficient, E, B_re, B_im);

    F.ProjectCoefficient(F_coeff);*/
    
    
    //------------------------------------------------------------
    //----------Computation of Lorentz force ends-----------------
    //------------------------------------------------------------
    
    
    //----------Paraview visualization--------------------------
    ParaViewDataCollection paraview_dc(basename, pmesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.SetDataFormat(VTKFormat::BINARY);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.RegisterField("E_re", &E.real());
    paraview_dc.RegisterField("E_im", &E.imag());
    //paraview_dc.RegisterField("B_re", &B_re);
    //paraview_dc.RegisterField("B_im", &B_im);
    //paraview_dc.RegisterField("J_poloidal", &E_Ambipolar);
    //paraview_dc.RegisterField("grad_Pe", &grad_Pe);
    //paraview_dc.RegisterField("ne", &ne);        
    //paraview_dc.RegisterField("phi", &phi);
    paraview_dc.RegisterField("W", &W_final);    
    //paraview_dc.RegisterField("F", &F);    
    //paraview_dc.RegisterField("sigma", &conductivity);
    paraview_dc.RegisterField("div_E_re", &div_Ere);
    paraview_dc.RegisterField("div_E_im", &div_Eim);
    paraview_dc.Save();
    
   
   return 0;
    
}

namespace mfem
{

namespace electromagnetics
{
  
void current_source(const Vector &x, Vector &Js)
{
  //Azimuthal current source
  /*double offset = 0.0;
  double r = sqrt(x[0]*x[0] + x[1]*x[1]);

    if (r >= 0.5+offset)
    {
    double theta = atan2(x[1], x[0]);
    double J0 = -10000.0;

    Js[2] = 0.0;
    Js[0] = -J0*sin(theta);
    Js[1] = J0*cos(theta);
    }
    else
    {
    Js[2] = 0.0;
    Js[0] = 0.0;
    Js[1] = 0.0;
    }*/

  //Axial current source
  /*double offset = 0.0;
  double r = sqrt(x[0]*x[0] + x[1]*x[1]);
  Js[0] = 0.0;
  Js[1] = 0.0;
  if ( r <= (0.5+offset))
    {
      Js[2] = -10000.0;
    }
  else
    {
      Js[2] = 0.0;
      }*/

  //Circular coil source with radius of coil = R and
  //radius of c/s of coil ri located at x = xi
  /*double R = 0.5;
  double ri = 0.05;
  double xi = 1.0;

  double J0 = -1.0e4;
  double r = sqrt(x[1]*x[1] + x[2]*x[2]);
  double theta;
  
  if (r <= (R+ri) && r >= (R-ri))
    {
      if (x[0] <= (xi + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (xi - sqrt(ri*ri - (r-R)*(r-R))))
	{
	  theta = atan2(x[2], x[1]);
	  Js[1] = -J0*sin(theta);
	  Js[2] = J0*cos(theta);
	}
      else
	{
	  Js[1] = 0.0;
	  Js[2] = 0.0;
	}	
    }
  else
	{
	  Js[1] = 0.0;
	  Js[2] = 0.0;
	}
      
	Js[0] = 0.0;*/
  
  
  //Circular coil source with radius of coil = R and
  //radius of c/s of coil ri located at x = xi
  double R = 0.11;
  double ri = 0.01;
  double x1 = 0.127;
  double x2 = 0.177;
  double x3 = 0.227;
  double x4 = 0.277;
  double x5 = 0.327;
  double x6 = 0.377;
  
  double J0 = -1.0e11;
  double r = sqrt(x[1]*x[1] + x[2]*x[2]);
  double theta;
  
  if (r <= (R+ri) && r >= (R-ri))
    {
      if (x[0] <= (x1 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x1 - sqrt(ri*ri - (r-R)*(r-R))))
	{
	  theta = atan2(x[2], x[1]);
	  Js[1] = -J0*sin(theta);
	  Js[2] = J0*cos(theta);
	}
      else if (x[0] <= (x2 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x2 - sqrt(ri*ri - (r-R)*(r-R))))
	{

	  theta = atan2(x[2], x[1]);
	  Js[1] = -J0*sin(theta);
	  Js[2] = J0*cos(theta);
	}
      else if (x[0] <= (x3 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x3 - sqrt(ri*ri - (r-R)*(r-R))))
	{
	  theta = atan2(x[2], x[1]);
	  Js[1] = -J0*sin(theta);
	  Js[2] = J0*cos(theta);
	}
      else if (x[0] <= (x4 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x4 - sqrt(ri*ri - (r-R)*(r-R))))
	{
	  theta = atan2(x[2], x[1]);
	  Js[1] = -J0*sin(theta);
	  Js[2] = J0*cos(theta);
	}
      else if (x[0] <= (x5 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x5 - sqrt(ri*ri - (r-R)*(r-R))))
	{
	  theta = atan2(x[2], x[1]);
	  Js[1] = -J0*sin(theta);
	  Js[2] = J0*cos(theta);
	}
      else if (x[0] <= (x6 + sqrt(ri*ri - (r-R)*(r-R))) && x[0] >= (x6 - sqrt(ri*ri - (r-R)*(r-R))))
	{
	  theta = atan2(x[2], x[1]);
	  Js[1] = -J0*sin(theta);
	  Js[2] = J0*cos(theta);
	}      
      else
	{
	  Js[1] = 0.0;
	  Js[2] = 0.0;
	}
      
      
    }
  else
    {
      Js[1] = 0.0;
      Js[2] = 0.0;
      }

      Js[0] = 0.0;

//Circular coil source with rectangular c/s
      /*double offset = 0.5;
  double r1 = 0.4 + offset;
  double r2 = 0.5 + offset;
  double x1 = -0.1;
  double x2 = 0.1;
  double J0 = -1000000.0;
  double r = sqrt(x[1]*x[1] + x[2]*x[2]);

  if (r <= r2 && r >= r1)
    {
      if (x[0] <= x2 && x[0] >= x1)
	{
	  double theta = atan2(x[2], x[1]);
	  Js[1] = -J0*sin(theta);
	  Js[2] = J0*cos(theta);
	}
      else
	{
	  Js[1] = 0.0;
	  Js[2] = 0.0;
	}
    }
  else
    {
      Js[1] = 0.0;
      Js[2] = 0.0;
    }

    Js[0] = 0.0;*/

  //No current source
  //Js = 0.0;
  }

double voltage_bc(const Vector &x)
{

  double P = 200.0;

  return P;
}

double unit_coeff(const Vector &x)
{

  double r = sqrt(x[1]*x[1] + x[2]*x[2]);
  double P;
  if (r <= 0.08)
    {
      P = 1.0;
    }
  else
    {
      P = 0.0;
    }

  return P;
}


double sigma_conductor(const Vector &x)
{

  double r0, c;
  r0 = 0.0;
  c = 0.1;
  double sig;
  sig = sigma;

  return sig;
  
}

double sigma_non_conductor(const Vector &x)
{

  double sig;

  sig = sigmaAir;
  
  return sig;
  
}

double sigma_profile(const Vector &x)
{

  double r = sqrt(x[1]*x[1] + x[2]*x[2]);
  double sig1, sig2, sig, a, b, c, d, m_x_L, m_x_R, m_y, dx_norm, dy_norm;
  int left_x, right_x, left_y, right_y;

  
  left_x = int((x[0] - xc_min)/dx);
  right_x = left_x + 1;

  left_y = int((r - yc_min)/dy);
  right_y = left_y + 1;

  //Normalized x & y difference
  dx_norm = (x[0]  - xc[left_x])/dx;
  dy_norm = (r  - yc[left_y])/dy;

  //Linear interpolation
  a      = sigma_dist[left_y][left_x];
  b      = sigma_dist[left_y][right_x];
  c      = sigma_dist[right_y][left_x];
  d      = sigma_dist[right_y][right_x];
  m_x_L      = b - a;
  m_x_R      = d - c;

  sig1 = a + m_x_L*dx_norm;
  sig2 = c + m_x_R*dx_norm;

  m_y = sig2 - sig1;
  sig = sig1 + m_y*dy_norm;

  if (sig < sigmaAir)
    {
      sig = sigmaAir;
      }

    
  return sig;
  
}

double Pe_func(const Vector &x)
{

  double r = sqrt(x[1]*x[1] + x[2]*x[2]);
  double P_e1, P_e2, P_e, a, b, c, d, m_x_L, m_x_R, m_y, dx_norm, dy_norm;
  int left_x, right_x, left_y, right_y;

  if (r <= 0.08)
    {
      left_x = int((x[0] - e_xc_min)/e_dx);
      right_x = left_x + 1;

      left_y = int((r - e_yc_min)/e_dy);
      right_y = left_y + 1;

      //Normalized x & y difference
      dx_norm = (x[0]  - e_xc[left_x])/e_dx;
      dy_norm = (r  - e_yc[left_y])/e_dy;

      //Linear interpolation
      a      = Pe[left_y][left_x];
      b      = Pe[left_y][right_x];
      c      = Pe[right_y][left_x];
      d      = Pe[right_y][right_x];
      m_x_L      = b - a;
      m_x_R      = d - c;

      P_e1 = a + m_x_L*dx_norm;
      P_e2 = c + m_x_R*dx_norm;

      m_y = P_e2 - P_e1;
      P_e = P_e1 + m_y*dy_norm;
    }
  else
    {
      P_e = 0.0;
    }

  /*if (P_e <= 1.0e-10)
    {
      P_e = 1.0e-10;
      }*/
  
  return P_e;
  
}

double ne_func(const Vector &x)
{

  double r = sqrt(x[1]*x[1] + x[2]*x[2]);
  double n_e1, n_e2, n_e, a, b, c, d, m_x_L, m_x_R, m_y, dx_norm, dy_norm;
  int left_x, right_x, left_y, right_y;


  if (r <= 0.08)
    {
      left_x = int((x[0] - e_xc_min)/e_dx);
      right_x = left_x + 1;

      left_y = int((r - e_yc_min)/e_dy);
      right_y = left_y + 1;

      //Normalized x & y difference
      dx_norm = (x[0]  - e_xc[left_x])/e_dx;
      dy_norm = (r  - e_yc[left_y])/e_dy;

      //Linear interpolation
      a      = ne[left_y][left_x];
      b      = ne[left_y][right_x];
      c      = ne[right_y][left_x];
      d      = ne[right_y][right_x];
      m_x_L      = b - a;
      m_x_R      = d - c;

      n_e1 = a + m_x_L*dx_norm;
      n_e2 = c + m_x_R*dx_norm;

      m_y = n_e2 - n_e1;
      n_e = n_e1 + m_y*dy_norm;
    }
  else
    {
      n_e = 1.0e-40;
    }

  if (n_e < 1.0e-40)
    {
      n_e = 1.0e-40;
    }
  
  return n_e;
  
}

} // namespace electromagnetics

} // namespace mfem
