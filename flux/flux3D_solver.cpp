// This file is part of the Finite-element soLver for Unsteady electromagnetiX (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=====================================================================================================

#include "flux3D_solver.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

  
  namespace electromagnetics
  {

    MeshDependentCoefficient::MeshDependentCoefficient(
						       std::map<int, std::function<double(const Vector &)>> &inputMap, double scale)
      : Coefficient()
    {
      // make a copy of the magic attribute-value map for later use
      materialMap = new std::map<int, std::function<double(const Vector &)>>(inputMap);
      scaleFactor = scale;
    }

    MeshDependentCoefficient::MeshDependentCoefficient(
						       MeshDependentCoefficient &cloneMe)
      : Coefficient()
    {
      // make a copy of the magic attribute-value map for later use
      materialMap = new std::map<int, std::function<double(const Vector &)>>(*(cloneMe.materialMap));
      scaleFactor = cloneMe.scaleFactor;
    }

    double MeshDependentCoefficient::Eval(ElementTransformation &T,
					  const IntegrationPoint &ip)
    {
      // given the attribute, extract the coefficient value from the map
      std::map<int, std::function<double(const Vector &)>>::iterator it;
      int thisAtt = T.Attribute;
      double x[3];
      Vector transip(x, 3);
      it = materialMap->find(thisAtt);
      T.Transform(ip, transip);

      if (it != materialMap->end())
	{
	  return it->second(transip);
	}
      else
	{
	  return 0.0; // avoid compile warning
	  std::cerr << "MeshDependentCoefficient attribute " << thisAtt
		    << " not found" << std::endl;
	  mfem_error();
	}
      
    }

    double JouleHeatingCoefficient::Eval(ElementTransformation &T,
					 const IntegrationPoint &ip)
    {
      Vector E_re, E_im;
      double thisSigma_coefficient;
      E.real().GetVectorValue(T, ip, E_re);
      E.imag().GetVectorValue(T, ip, E_im);
      thisSigma_coefficient = sigma_coefficient.Eval(T, ip);
      return 0.5*thisSigma_coefficient*(E_re*E_re + E_im*E_im);
      
    }

    void AmbipolarCoefficient::Eval(Vector &V, ElementTransformation &T,
					 const IntegrationPoint &ip)
    {
      double nb_density_e;
      Vector gradient_P;
      grad_Pe.GetVectorValue(T, ip, gradient_P);
      nb_density_e = ne.GetValue(T, ip);
      
      V[0] = gradient_P[0]/nb_density_e;
      V[1] = gradient_P[1]/nb_density_e;
      V[2] = gradient_P[2]/nb_density_e;
      
    }

    void sigmaEreCoefficient::Eval(Vector &V, ElementTransformation &T,
					 const IntegrationPoint &ip)
    {
      double sig;
      Vector Ereal(3);
      E.real().GetVectorValue(T, ip, Ereal);
      sig = sigma_coefficient.Eval(T, ip);
      
      V[0] = sig*Ereal[0];
      V[1] = sig*Ereal[1];
      V[2] = sig*Ereal[2];
      
    }

    void sigmaEimCoefficient::Eval(Vector &V, ElementTransformation &T,
					 const IntegrationPoint &ip)
    {
      double sig;
      Vector Eimag(3);
      E.imag().GetVectorValue(T, ip, Eimag);
      sig = sigma_coefficient.Eval(T, ip);
      
      V[0] = sig*Eimag[0];
      V[1] = sig*Eimag[1];
      V[2] = sig*Eimag[2];
      
    }

    void LorentzForceCoefficient::Eval(Vector &V, ElementTransformation &T,
					 const IntegrationPoint &ip)
    {
      double sig;
      Vector Ereal, Eimag, Breal, Bimag;
      E.real().GetVectorValue(T, ip, Ereal);
      E.imag().GetVectorValue(T, ip, Eimag);
      B_re.GetVectorValue(T, ip, Breal);
      B_im.GetVectorValue(T, ip, Bimag);      
      sig = sigma_coefficient.Eval(T, ip);
      
      V[0] = 0.5*sig*(Ereal[1]*Breal[2] + Eimag[1]*Bimag[2] - Ereal[2]*Breal[1] - Eimag[2]*Bimag[1]);
      V[1] = - 0.5*sig*(Ereal[0]*Breal[2] + Eimag[0]*Bimag[2] - Ereal[2]*Breal[0] - Eimag[2]*Bimag[0]);
      V[2] = 0.5*sig*(Ereal[0]*Breal[1] + Eimag[0]*Bimag[1] - Ereal[1]*Breal[0] - Eimag[1]*Bimag[0]);
      
      
    }
    
  } // namespace electromagnetics

} // namespace mfem

#endif // MFEM_USE_MPI
