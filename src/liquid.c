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

int apecss_liquid_setdefaultoptions(struct APECSS_Liquid *Liquid)
{
  Liquid->Type = APECSS_LIQUID_NEWTONIAN;
  Liquid->EoS = 0;
  Liquid->Gamma = 7.15;
  Liquid->B = 0.0;
  Liquid->b = 0.0;
  Liquid->mu = 0.001;
  Liquid->pref = 1.0e5;
  Liquid->rhoref = 1000.0;
  Liquid->cref = 1500.0;
  Liquid->G = 0.0;
  Liquid->eta = 0.0;
  Liquid->lambda = 0.0;

  Liquid->get_density = apecss_liquid_density_fixed;
  Liquid->get_soundspeed = apecss_liquid_soundspeed_fixed;
  Liquid->get_enthalpy = apecss_liquid_enthalpy_quasiacoustic;
  Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall;
  Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_expl;

  return (0);
}

int apecss_liquid_processoptions(struct APECSS_Liquid *Liquid)
{
  // Set the function pointers that define the equation of state of the liquid
  if (Liquid->EoS)
  {
    if (Liquid->EoS == APECSS_LIQUID_TAIT)
    {
      Liquid->get_density = apecss_liquid_density_tait;
      Liquid->get_soundspeed = apecss_liquid_soundspeed_tait;
      Liquid->get_enthalpy = apecss_liquid_enthalpy_tait;

      // Compute the reference coefficient for the liquid
      Liquid->Kref = Liquid->rhoref / (APECSS_POW(Liquid->pref + Liquid->B, 1.0 / Liquid->Gamma));
    }
    else if (Liquid->EoS == APECSS_LIQUID_NASG)
    {
      Liquid->get_density = apecss_liquid_density_nasg;
      Liquid->get_soundspeed = apecss_liquid_soundspeed_nasg;
      Liquid->get_enthalpy = apecss_liquid_enthalpy_nasg;

      // Compute the reference coefficient for the liquid
      Liquid->Kref = Liquid->rhoref / (APECSS_POW(Liquid->pref + Liquid->B, 1.0 / Liquid->Gamma) * (1.0 - Liquid->b * Liquid->rhoref));
    }
    else
    {
      apecss_erroronscreen(-1, "Unknown equation of state defined for the liquid for the Gilmore model.");
    }

    if (Liquid->Kref <= 0.0) apecss_erroronscreen(-1, "K reference factor of the liquid is unphysical (Kref <= 0).");
  }
  else
  {
    Liquid->get_density = apecss_liquid_density_fixed;
    Liquid->get_soundspeed = apecss_liquid_soundspeed_fixed;
    Liquid->get_enthalpy = apecss_liquid_enthalpy_quasiacoustic;
  }

  // Set the function pointers that define the type of the liquid
  Liquid->get_pressure_viscous = apecss_liquid_pressure_viscous;
  Liquid->get_pressurederivative_viscous_expl = apecss_liquid_pressurederivative_viscous_expl;
  Liquid->get_pressurederivative_viscous_impl = apecss_liquid_pressurederivative_viscous_impl;

  if (Liquid->Type == APECSS_LIQUID_NEWTONIAN)
  {
    Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall;
    Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_expl;
  }
  else if (Liquid->Type == APECSS_LIQUID_KELVINVOIGT)
  {
    Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall_kelvinvoigt;
    Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_explkelvinvoigt;
  }
  else if (Liquid->Type == APECSS_LIQUID_ZENER)
  {
    Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall_zener;
    Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_zener;

    Liquid->get_pressurederivative_viscous_impl = apecss_liquid_pressurederivative_viscous_nonimpl;
  }
  else if (Liquid->Type == APECSS_LIQUID_OLDROYDB)
  {
    Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall_oldroydb;
    Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_exploldroydb;
  }

  return (0);
}

// -------------------------------------------------------------------
// PROPERTIES
// -------------------------------------------------------------------

APECSS_FLOAT apecss_liquid_density_fixed(APECSS_FLOAT p, struct APECSS_Liquid *Liquid) { return (Liquid->rhoref); }

APECSS_FLOAT apecss_liquid_density_tait(APECSS_FLOAT p, struct APECSS_Liquid *Liquid)
{
  return (Liquid->Kref * APECSS_POW((p + Liquid->B), (1.0 / Liquid->Gamma)));
}

APECSS_FLOAT apecss_liquid_density_nasg(APECSS_FLOAT p, struct APECSS_Liquid *Liquid)
{
  APECSS_FLOAT peff = APECSS_POW((p + Liquid->B), (1.0 / Liquid->Gamma));
  return (Liquid->Kref * peff / (1.0 + Liquid->b * Liquid->Kref * peff));
}

APECSS_FLOAT apecss_liquid_soundspeed_fixed(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid) { return (Liquid->cref); }

APECSS_FLOAT apecss_liquid_soundspeed_tait(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid)
{
  return (APECSS_SQRT(Liquid->Gamma * (p + Liquid->B) / rho));
}

APECSS_FLOAT apecss_liquid_soundspeed_nasg(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid)
{
  return (APECSS_SQRT(Liquid->Gamma * (p + Liquid->B) / (rho * (1.0 - Liquid->b * rho))));
}

APECSS_FLOAT apecss_liquid_enthalpy_quasiacoustic(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid) { return (p / rho); }

APECSS_FLOAT apecss_liquid_enthalpy_tait(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid)
{
  return (Liquid->Gamma * (p + Liquid->B) / ((Liquid->Gamma - 1.0) * rho));
}

APECSS_FLOAT apecss_liquid_enthalpy_nasg(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid)
{
  APECSS_FLOAT fact = Liquid->Gamma / (Liquid->Gamma - 1.0);
  return (fact * (p + Liquid->B) / rho + Liquid->b * (p - fact * (p + Liquid->B)));
}

// -------------------------------------------------------------------
// BUBBLE WALL PRESSURE
// -------------------------------------------------------------------

APECSS_FLOAT apecss_liquid_pressure_bubblewall(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble->Interface) -
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble) - Bubble->Interface->get_pressure_viscous(Sol[1], Sol[0], Bubble));
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_kelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble->Interface) -
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble) - Bubble->Interface->get_pressure_viscous(Sol[1], Sol[0], Bubble) -
          4.0 * APECSS_ONETHIRD * Bubble->Liquid->G * (1.0 - APECSS_POW3(Bubble->R0) / APECSS_POW3(Sol[1])));
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble->Interface) + 3.0 * Sol[2] -
          Bubble->Interface->get_pressure_viscous(Sol[1], Sol[0], Bubble));
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_oldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble->Interface) -
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble) - Bubble->Interface->get_pressure_viscous(Sol[1], Sol[0], Bubble) + Sol[2] + Sol[3]);
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_expl(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble->Interface) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) +
          Bubble->Interface->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_explkelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble->Interface) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) +
          Bubble->Interface->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) -
          4.0 * Bubble->Liquid->G * APECSS_POW3(Bubble->R0) * Sol[0] / APECSS_POW4(Sol[1]));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble->Interface) +
          3.0 * apecss_viscoelastic_zenervarsigma_ode(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_exploldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble->Interface) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) +
          Bubble->Interface->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) + apecss_viscoelastic_oldroydb1_ode(Sol, t, Bubble) +
          apecss_viscoelastic_oldroydb2_ode(Sol, t, Bubble));
}

// -------------------------------------------------------------------
// VISCOUS PRESSURE
// -------------------------------------------------------------------

APECSS_FLOAT apecss_liquid_pressure_viscous(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble) { return (4.0 * Bubble->Liquid->mu * U / R); }

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_expl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (4.0 * Bubble->Liquid->mu * APECSS_POW2(U / R));
}

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_impl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble) { return (4.0 * Bubble->Liquid->mu / R); }

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_nonimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble) { return (0.0); }
