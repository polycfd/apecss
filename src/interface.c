// This source file is part of APECSS, an open-source software toolbox
// for the computation of pressure-driven bubble dynamics and acoustic
// emissions in spherical symmetry.
//
// Copyright (C) 2022 The APECSS Developers
//
// The APECSS Developers are listed in the README.md file available in
// the GitHub repository at https://github.com/polycfd/apecss.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "apecss.h"

// -------------------------------------------------------------------
// SET OPTIONS
// -------------------------------------------------------------------

int apecss_interface_setdefaultoptions(struct APECSS_Interface *Interface)
{
  Interface->sigma = 0.0;
  Interface->LipidCoatingModel = APECSS_LIPIDCOATING_NONE;
  Interface->sigma0 = 0.02;
  Interface->Elasticity = 0.5;
  Interface->Viscosity = 7.5e-9;
  Interface->Rbuck = 1.0e-6;
  Interface->Rrupt = 1.0e-6;

  Interface->get_surfacetension = apecss_interface_surfacetension_clean;
  Interface->get_pressure_surfacetension = apecss_interface_surfacetensionpressure_clean;
  Interface->get_pressurederivative_surfacetension = apecss_interface_surfacetensionpressurederivative_clean;

  Interface->get_pressure_viscous = apecss_interface_pressure_viscous_clean;
  Interface->get_pressurederivative_viscous_expl = apecss_interface_pressurederivative_viscous_cleanexpl;
  Interface->get_pressurederivative_viscous_impl = apecss_interface_pressurederivative_viscous_cleanimpl;

  return (0);
}

int apecss_interface_processoptions(struct APECSS_Interface *Interface)
{
  // Set the appropriate function pointers associated with bubble coating
  if (Interface->LipidCoatingModel & APECSS_LIPIDCOATING_MARMOTTANT)
  {
    if (Interface->LipidCoatingModel & APECSS_LIPIDCOATING_GOMPERTZFUNCTION)
    {
      Interface->get_surfacetension = apecss_interface_surfacetension_gompertzmarmottant;
      Interface->get_pressure_surfacetension = apecss_interface_surfacetensionpressure_gompertzmarmottant;
      Interface->get_pressurederivative_surfacetension = apecss_interface_surfacetensionpressurederivative_gompertzmarmottant;
    }
    else
    {
      Interface->get_surfacetension = apecss_interface_surfacetension_marmottant;
      Interface->get_pressure_surfacetension = apecss_interface_surfacetensionpressure_marmottant;
      Interface->get_pressurederivative_surfacetension = apecss_interface_surfacetensionpressurederivative_marmottant;
    }

    Interface->get_pressure_viscous = apecss_interface_pressure_viscous_marmottant;
    Interface->get_pressurederivative_viscous_expl = apecss_interface_pressurederivative_viscous_marmottantexpl;
    Interface->get_pressurederivative_viscous_impl = apecss_interface_pressurederivative_viscous_marmottantimpl;
  }
  else
  {
    Interface->get_surfacetension = apecss_interface_surfacetension_clean;
    Interface->get_pressure_surfacetension = apecss_interface_surfacetensionpressure_clean;
    Interface->get_pressurederivative_surfacetension = apecss_interface_surfacetensionpressurederivative_clean;

    Interface->get_pressure_viscous = apecss_interface_pressure_viscous_clean;
    Interface->get_pressurederivative_viscous_expl = apecss_interface_pressurederivative_viscous_cleanexpl;
    Interface->get_pressurederivative_viscous_impl = apecss_interface_pressurederivative_viscous_cleanimpl;
  }

  return (0);
}

// -------------------------------------------------------------------
// SURFACE TENSION COEFFICIENT
// -------------------------------------------------------------------

APECSS_FLOAT apecss_interface_surfacetension_clean(APECSS_FLOAT R, struct APECSS_Interface *Interface) { return (Interface->sigma); }

APECSS_FLOAT apecss_interface_surfacetension_marmottant(APECSS_FLOAT R, struct APECSS_Interface *Interface)
{
  APECSS_FLOAT sigma;
  if (R < Interface->Rbuck)
    sigma = 0.0;
  else if (R > Interface->Rrupt)
    sigma = Interface->sigma;
  else
    sigma = Interface->Elasticity * (APECSS_POW2(R) / APECSS_POW2(Interface->Rbuck) - 1.0);

  return (sigma);
}

APECSS_FLOAT apecss_interface_surfacetension_gompertzmarmottant(APECSS_FLOAT R, struct APECSS_Interface *Interface)
{
  return (Interface->sigma * APECSS_EXP(-Interface->GompertzB * APECSS_EXP(Interface->GompertzC * (1.0 - R / Interface->Rbuck))));
}

APECSS_FLOAT apecss_interface_surfacetensionderivative_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface)
{
  APECSS_FLOAT dot_sigma;

  if (R < Interface->Rbuck)
    dot_sigma = 0.0;
  else if (R > Interface->Rrupt)
    dot_sigma = 0.0;
  else
    dot_sigma = 2.0 * Interface->Elasticity * R * U / APECSS_POW2(Interface->Rbuck);

  return (dot_sigma);
}

APECSS_FLOAT apecss_interface_surfacetensionderivative_gompertzmarmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface)
{
  return (apecss_interface_surfacetension_gompertzmarmottant(R, Interface) * U * Interface->GompertzB * Interface->GompertzC *
          APECSS_EXP(Interface->GompertzC * (1.0 - R / Interface->Rbuck)) / Interface->Rbuck);
}

// -------------------------------------------------------------------
// SURFACE TENSION PRESSURE
// -------------------------------------------------------------------

APECSS_FLOAT apecss_interface_surfacetensionpressure_clean(APECSS_FLOAT R, struct APECSS_Interface *Interface) { return (2.0 * Interface->sigma / R); }

APECSS_FLOAT apecss_interface_surfacetensionpressure_marmottant(APECSS_FLOAT R, struct APECSS_Interface *Interface)
{
  return (2.0 * apecss_interface_surfacetension_marmottant(R, Interface) / R);
}

APECSS_FLOAT apecss_interface_surfacetensionpressure_gompertzmarmottant(APECSS_FLOAT R, struct APECSS_Interface *Interface)
{
  return (2.0 * apecss_interface_surfacetension_gompertzmarmottant(R, Interface) / R);
}

APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface)
{
  return (2.0 * Interface->sigma * U / APECSS_POW2(R));
}

APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface)
{
  return (2.0 * (apecss_interface_surfacetension_marmottant(R, Interface) * U / APECSS_POW2(R) -
                 apecss_interface_surfacetensionderivative_marmottant(R, U, Interface) / R));
}

APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_gompertzmarmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface)
{
  return (2.0 * (apecss_interface_surfacetension_gompertzmarmottant(R, Interface) * U / APECSS_POW2(R) -
                 apecss_interface_surfacetensionderivative_gompertzmarmottant(R, U, Interface) / R));
}

// -------------------------------------------------------------------
// VISCOUS PRESSURE
// -------------------------------------------------------------------

APECSS_FLOAT apecss_interface_pressure_viscous_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble) { return (0.0); }

APECSS_FLOAT apecss_interface_pressure_viscous_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (4.0 * Bubble->Interface->Viscosity * U / APECSS_POW2(R));
}

APECSS_FLOAT apecss_interface_pressurederivative_viscous_cleanexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble) { return (0.0); }

APECSS_FLOAT apecss_interface_pressurederivative_viscous_marmottantexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (8.0 * Bubble->Interface->Viscosity * APECSS_POW2(U) / APECSS_POW3(R));
}

APECSS_FLOAT apecss_interface_pressurederivative_viscous_cleanimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble) { return (0.0); }

APECSS_FLOAT apecss_interface_pressurederivative_viscous_marmottantimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble)
{
  return (4.0 * Bubble->Interface->Viscosity / APECSS_POW2(R));
}