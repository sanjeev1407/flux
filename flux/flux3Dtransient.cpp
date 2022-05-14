// This file is part of the Finite-element soLver for Unsteady electromagnetiX (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=====================================================================================================

#include "electric_solver.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;
using namespace mfem::common;
using namespace mfem::electromagnetics;

static double wj = 0.0;

// Initialize variables used in E_solver.cpp
int electromagnetics::SOLVER_PRINT_LEVEL = 0;
int electromagnetics::STATIC_COND        = 0;

int main(int argc, char *argv[])
{
   // 1. Initialize MPI.
   MPI_Session mpi(argc, argv);
   int myid = mpi.WorldRank();


   // 2. Parse command-line options.
   const char *mesh_file = "conductor_air.msh";
   int order = 2;
   int ode_solver_type = 1;
   double t_final = 1.0;
   double dt = 0.1;
   double mu = 1.25e-6;
   double sigma = 100.0;
   double w = 10.0;
   int gfprint = 0;
   const char *basename = "Electric";
   int debug = 0;
   int vis_steps = 1;
   

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Order (degree) of the finite elements.");
   args.AddOption(&ode_solver_type, "-s", "--ode-solver",
                  "ODE solver: 1 - Backward Euler, 2 - SDIRK2, 3 - SDIRK3\n\t."
                  "\t   22 - Mid-Point, 23 - SDIRK23, 34 - SDIRK34.");
   args.AddOption(&t_final, "-tf", "--t-final",
                  "Final time; start time is 0.");
   args.AddOption(&dt, "-dt", "--time-step",
                  "Time step.");
   args.AddOption(&mu, "-mu", "--permeability",
                  "Magnetic permeability coefficient.");
   args.AddOption(&sigma, "-cnd", "--sigma",
                  "Conductivity coefficient.");
   args.AddOption(&basename, "-k", "--outputfilename",
                  "Name of the visit dump files");
   args.AddOption(&gfprint, "-print", "--print",
                  "Print results (grid functions) to disk.");
   args.AddOption(&debug, "-debug", "--debug",
                  "Print matrices and vectors to disk");
   args.AddOption(&SOLVER_PRINT_LEVEL, "-hl", "--hypre-print-level",
                  "Hypre print level");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.Root())
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (mpi.Root())
   {
      args.PrintOptions(cout);
   }

   wj = w;
   double sigmaAir     = 1.0e-6 * sigma;
   // 3. Here material properties are assigned to mesh attributes.  
   
   std::map<int, double> sigmaMap;
      
   sigmaMap.insert(pair<int, double>(1, sigma));
   sigmaMap.insert(pair<int, double>(2, sigmaAir));
   sigmaMap.insert(pair<int, double>(3, sigmaAir));

   
   // 4. Read the serial mesh from the given mesh file on all processors. We can
   //    handle triangular, quadrilateral, tetrahedral and hexahedral meshes
   //    with the same code.
   Mesh *mesh;
   mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 5. Assign the boundary conditions
   Array<int> ess_bdr_zero(mesh->bdr_attributes.Max());
   Array<int> ess_bdr_cosine(mesh->bdr_attributes.Max());
   Array<int> phi_ess_bdr_zero(mesh->bdr_attributes.Max());
   Array<int> phi_ess_bdr_cosine(mesh->bdr_attributes.Max());


   // For hollow cylinder problem we have 4 surfaces
   // 1)Z=0 plane, 2)Z=1 plane, 3)outer circular plane, 4)Inner circular plane   

   //nxE = 0 essential bc
   ess_bdr_zero = 1;
   //ess_bdr_zero[0] = 1;
   //ess_bdr_zero[1] = 1;
   
   ess_bdr_cosine = 0;
   
   //phi=cos(w*t) on z=0 plane and phi=0 on z=1 plane
   phi_ess_bdr_zero = 0;
   phi_ess_bdr_zero[0] = 1;
   
   phi_ess_bdr_cosine = 0;
   phi_ess_bdr_cosine[1] = 1;
   
   // 6. Define the ODE solver used for time integration
   ODESolver *ode_solver;
   switch (ode_solver_type)
   {
      // Implicit L-stable methods
      case 1:  ode_solver = new BackwardEulerSolver; break;
      case 2:  ode_solver = new SDIRK23Solver(2); break;
      case 3:  ode_solver = new SDIRK33Solver; break;
      // Implicit A-stable methods (not L-stable)
      case 22: ode_solver = new ImplicitMidpointSolver; break;
      case 23: ode_solver = new SDIRK23Solver; break;
      case 34: ode_solver = new SDIRK34Solver; break;
      default:
         if (mpi.Root())
         {
            cout << "Unknown ODE solver type: " << ode_solver_type << '\n';
         }
         delete mesh;
         return 3;
   }

   
   // 8. Define a parallel mesh by a partitioning of the serial mesh. 
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   pmesh->ReorientTetMesh();

   // 12. Define the parallel finite element spaces. We use:
   //
   //     H(curl) for electric field,
   //     H(div) for magnetic flux,
   //     H(grad) for electrostatic potential
   
   
   // ND contains Nedelec "edge-centered" vector finite elements with continuous
   // tangential component.
   ND_FECollection HCurlFEC(order, dim);

   // RT contains Raviart-Thomas "face-centered" vector finite elements with
   // continuous normal component.
   RT_FECollection HDivFEC(order, dim);

   // H1 contains continuous "node-centered" Lagrange finite elements
   H1_FECollection HGradFEC(order, dim);

   
   ParFiniteElementSpace HCurlFESpace(pmesh, &HCurlFEC);
   ParFiniteElementSpace  HDivFESpace(pmesh, &HDivFEC);
   ParFiniteElementSpace  HGradFESpace(pmesh, &HGradFEC);
  
   // The terminology is TrueVSize is the unique (non-redundant) number of dofs
   HYPRE_Int glob_size_nd = HCurlFESpace.GlobalTrueVSize();
   HYPRE_Int glob_size_rt = HDivFESpace.GlobalTrueVSize();
   HYPRE_Int glob_size_h1 = HGradFESpace.GlobalTrueVSize();

   if (mpi.Root())
   {
      cout << "Number of Electric Field unknowns:    " << glob_size_nd << endl;
      cout << "Number of Magnetic Field unknowns:    " << glob_size_rt << endl;
      cout << "Number of Electrostatic potential unknowns:     " << glob_size_h1 << endl;
   }

   int Vsize_nd = HCurlFESpace.GetVSize();
   int Vsize_rt = HDivFESpace.GetVSize();
   int Vsize_h1 = HGradFESpace.GetVSize();
   
   // the big BlockVector stores the fields as
   //    1 E field
   //    2 B field
   //    3 P field
   
   Array<int> true_offset(4);
   true_offset[0] = 0;
   true_offset[1] = true_offset[0] + Vsize_nd;
   true_offset[2] = true_offset[1] + Vsize_rt;
   true_offset[3] = true_offset[2] + Vsize_h1;
   
   // The BlockVector is a large contiguous chunk of memory for storing required
   // data for the hypre vectors, in this case: the E-field HCurl, the B-field HDiv
   // and P field HGrad
   BlockVector F(true_offset);

   // grid functions E, B, P
   ParGridFunction E_gf, B_gf, P_gf;
   E_gf.MakeRef(&HCurlFESpace,F,true_offset[0]);
   B_gf.MakeRef(&HDivFESpace,F, true_offset[1]);
   P_gf.MakeRef(&HGradFESpace,F, true_offset[2]);
   
   // 14. Initialize the E_Diffusion operator
   ElectricDiffusionEOperator oper(true_offset[3], HCurlFESpace,
                                   HDivFESpace, HGradFESpace, ess_bdr_zero,  
                                   ess_bdr_cosine, phi_ess_bdr_zero,
				   phi_ess_bdr_cosine, mu, sigmaMap);
   
   // This function initializes all the fields to zero or some provided IC
   oper.Init(F);
   

   //----------Paraview visualization--------------------------
   ParaViewDataCollection paraview_dc(basename, pmesh);
   paraview_dc.SetPrefixPath("ParaView");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("E", &E_gf);
   paraview_dc.RegisterField("B", &B_gf);
   paraview_dc.RegisterField("P", &P_gf);
   paraview_dc.Save();
   

   // 15. Perform time Integration (looping over the time iterations, ti, with a
   //     time-step dt). The object oper is the ElectricDiffusionOperator which
   //     has a Mult() method and an ImplicitSolve() method which are used by
   //     the time integrators.
   ode_solver->Init(oper);
   double t = 0.0;

   bool last_step = false;
   for (int ti = 1; !last_step; ti++)
   {
      if (t + dt >= t_final - dt/2)
      {
         last_step = true;
      }

      // F is the vector of dofs, t is the current time, and dt is the time step
      // to advance.
      ode_solver->Step(F, t, dt);

      if (debug == 1)
      {
         oper.Debug(basename,t);
      }

      

   if (last_step || (ti % vis_steps) == 0)
      {
         
         if (mpi.Root())
         {
            cout << fixed;
            cout << "step " << setw(6) << ti << ",\tt = " << setw(6)
                 << setprecision(3) << t << endl;
         }

         
      }
   if (last_step)
     {
       paraview_dc.SetCycle(ti);
       paraview_dc.SetTime(t);
       paraview_dc.Save();
     }
         
   }
         
   // 16. Free the used memory.
   delete ode_solver;
   delete pmesh;

   return 0;
}

namespace mfem
{

namespace electromagnetics
{

void e_tan_zero_bc(const Vector &x, Vector &E)
{
   E = 0.0;
}

void e_tan_cosine_bc(const Vector &x, double t, Vector &E)
{
  E[0] = cos(wj * t);
  E[1] = 0.0;
  E[2] = 0.0;

  
}

double cosine_p_bc(const Vector &x, double t)
{
    
    double P = 1000.0;
    
    return P*cos(wj * t);
}

double zero_p_bc(const Vector &x, double t)
{
  return 0.0;
}
  
void e_exact(const Vector &x, double t, Vector &E)
{
   E[0] = 0.0;
   E[1] = 0.0;
   E[2] = 0.0;
}

void b_exact(const Vector &x, double t, Vector &B)
{
   B[0] = 0.0;
   B[1] = 0.0;
   B[2] = 0.0;
}

void current_source(const Vector &x, double t, Vector &Js)
{
  //Azimuthal current source
  /*double r = sqrt(x[0]*x[0] + x[1]*x[1]);
    
  if (r >= 0.5)
    {
      double theta = atan2(x[1], x[0]);
      double J0 = 10000*cos(wj*t);
  
      Js[0] = -J0*sin(theta);
      Js[1] = J0*cos(theta);
      Js[2] = 0.0;
    }
  else
    {
      Js[0] = 0.0;
      Js[1] = 0.0;
      Js[2] = 0.0;
      }*/

  //Axial current source
  /*Js[0] = 0.0;
  Js[1] = 0.0;
  
  if ( r >= 0.5)
    {
      Js[2] = cos(wj*t);
    }
  else
    {
      Js[2] = 0.0;
      }*/

  //Circular coil source with radius of coil = R and
  //radius of c/s of coil r0 located at z=z0
  double R = 0.5;
  double r0 = 0.05;
  double z0 = 0.5;
  double J0 = 1000*cos(wj*t);  
  double r = sqrt(x[0]*x[0] + x[1]*x[1]);

  if (r <= (R+r0) && r >= (R-r0))
    {
      if (x[2] <= (z0 + sqrt(r0*r0 - (r-R)*(r-R))) && x[2] >= (z0 - sqrt(r0*r0 - (r-R)*(r-R))))
	{
	  double theta = atan2(x[1], x[0]);
	  Js[0] = -J0*sin(theta);
	  Js[1] = J0*cos(theta);
	}
      else
	{
	  Js[0] = 0.0;
	  Js[1] = 0.0;
	}
    }
  else
    {
      Js[0] = 0.0;
      Js[1] = 0.0;
    }

  Js[2] = 0.0;
}

} // namespace electromagnetics

} // namespace mfem

