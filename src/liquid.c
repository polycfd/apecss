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

int apecss_liquid_setdefaultoptions(struct APECSS_Bubble *Bubble)
{
  Bubble->Liquid = (struct APECSS_Liquid *) malloc(sizeof(struct APECSS_Liquid));

  Bubble->Liquid->Type = APECSS_LIQUID_NEWTONIAN;
  Bubble->Liquid->Gamma = 1.186;
  Bubble->Liquid->B = 0.0;
  Bubble->Liquid->b = 0.0;
  Bubble->Liquid->mu = 0.001;
  Bubble->Liquid->pref = -1.0e10;  // Either set by user input or set to ambient pressure
  Bubble->Liquid->rhoref = 1000.0;
  Bubble->Liquid->cref = 1500.0;
  Bubble->Liquid->G = 0.0;
  Bubble->Liquid->eta = 0.0;
  Bubble->Liquid->lambda = 0.0;

  Bubble->Liquid->get_density = apecss_liquid_density_fixed;
  Bubble->Liquid->get_soundspeed = apecss_liquid_soundspeed_fixed;
  Bubble->Liquid->get_enthalpy = apecss_liquid_enthalpy_quasiacoustic;
  Bubble->Liquid->get_pressure_infinity = apecss_liquid_pressure_infinity_noexcitation;
  Bubble->Liquid->get_pressurederivative_infinity = apecss_liquid_pressurederivative_infinity_noexcitation;
  Bubble->Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall;
  Bubble->Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_expl;
  Bubble->Liquid->get_pressure_viscous = apecss_liquid_pressure_viscous_clean;
  Bubble->Liquid->get_pressurederivative_viscous_expl = apecss_liquid_pressurederivative_viscous_cleanexpl;
  Bubble->Liquid->get_pressurederivative_viscous_impl = apecss_liquid_pressurederivative_viscous_cleanimpl;

  return 0;
}

int apecss_liquid_processoptions(struct APECSS_Bubble *Bubble)
{
  if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE)
  {
    Bubble->Liquid->get_density = apecss_liquid_density_nasg;
    Bubble->Liquid->get_soundspeed = apecss_liquid_soundspeed_nasg;
    Bubble->Liquid->get_enthalpy = apecss_liquid_enthalpy_nasg;
  }
  else
  {
    Bubble->Liquid->get_density = apecss_liquid_density_fixed;
    Bubble->Liquid->get_soundspeed = apecss_liquid_soundspeed_fixed;
    Bubble->Liquid->get_enthalpy = apecss_liquid_enthalpy_quasiacoustic;
  }

  if (Bubble->Liquid->pref < -Bubble->Liquid->B) Bubble->Liquid->pref = Bubble->p0;

  // Set the function pointers for the properties of the liquid at infinity
  if (Bubble->Excitation != NULL)
  {
    if (Bubble->Excitation->type == APECSS_EXCITATION_SIN)
    {
      Bubble->Liquid->get_pressure_infinity = apecss_liquid_pressure_infinity_sinexcitation;
      Bubble->Liquid->get_pressurederivative_infinity = apecss_liquid_pressurederivative_infinity_sinexcitation;
    }
    else
    {
      Bubble->Liquid->get_pressure_infinity = apecss_liquid_pressure_infinity_noexcitation;
      Bubble->Liquid->get_pressurederivative_infinity = apecss_liquid_pressurederivative_infinity_noexcitation;
    }
  }
  else
  {
    Bubble->Liquid->get_pressure_infinity = apecss_liquid_pressure_infinity_noexcitation;
    Bubble->Liquid->get_pressurederivative_infinity = apecss_liquid_pressurederivative_infinity_noexcitation;
  }

  // Set the function pointers that define the type of the liquid
  if (Bubble->Liquid->Type == APECSS_LIQUID_NEWTONIAN)
  {
    Bubble->Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall;
    Bubble->Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_expl;
  }
  else if (Bubble->Liquid->Type == APECSS_LIQUID_KELVINVOIGT)
  {
    Bubble->Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall_kelvinvoigt;
    Bubble->Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_explkelvinvoigt;
  }
  else if (Bubble->Liquid->Type == APECSS_LIQUID_ZENER)
  {
    Bubble->Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall_zener;
    Bubble->Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_zener;

    Bubble->Liquid->get_pressurederivative_viscous_impl = apecss_liquid_pressurederivative_viscous_nonimpl;
  }
  else if (Bubble->Liquid->Type == APECSS_LIQUID_OLDROYDB)
  {
    Bubble->Liquid->get_pressure_bubblewall = apecss_liquid_pressure_bubblewall_oldroydb;
    Bubble->Liquid->get_pressurederivative_bubblewall_expl = apecss_liquid_pressurederivative_bubblewall_exploldroydb;
  }

  Bubble->Liquid->Kref = Bubble->Liquid->rhoref / (APECSS_POW(Bubble->Liquid->pref + Bubble->Liquid->B, 1.0 / Bubble->Liquid->Gamma) *
                                                   (1.0 - Bubble->Liquid->b * Bubble->Liquid->rhoref));

  if (Bubble->Liquid->Kref <= 0.0)
  {
    char str[APECSS_STRINGLENGTH_SPRINTF];
    sprintf(str, "K reference factor of the liquid is unphysical (Kref <= 0).");
    apecss_erroronscreen(-1, str);
  }

  return 0;
}

// -------------------------------------------------------------------
// PROPERTIES
// -------------------------------------------------------------------

APECSS_FLOAT apecss_liquid_density_fixed(APECSS_FLOAT p, struct APECSS_Liquid *Liquid) { return (Liquid->rhoref); }

APECSS_FLOAT apecss_liquid_density_nasg(APECSS_FLOAT p, struct APECSS_Liquid *Liquid)
{
  APECSS_FLOAT peff = APECSS_POW((p + Liquid->B), (1.0 / Liquid->Gamma));
  return (Liquid->Kref * peff / (1.0 + Liquid->b * Liquid->Kref * peff));
}

APECSS_FLOAT apecss_liquid_soundspeed_fixed(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid) { return (Liquid->cref); }

APECSS_FLOAT apecss_liquid_soundspeed_nasg(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid)
{
  return (APECSS_SQRT(Liquid->Gamma * (p + Liquid->B) / (rho * (1.0 - Liquid->b * rho))));
}

APECSS_FLOAT apecss_liquid_enthalpy_quasiacoustic(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid) { return (p / rho); }

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
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble));
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_kelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble->Interface) -
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble) -
          4.0 * APECSS_ONETHIRD * Bubble->Liquid->G * (1.0 - APECSS_POW3(Bubble->R0) / APECSS_POW3(Sol[1])));
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble->Interface) + 3.0 * Sol[2]);
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_oldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble->Interface) -
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble) + Sol[2] + Sol[3]);
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_expl(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble->Interface) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_explkelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble->Interface) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) -
          4.0 * Bubble->Liquid->G * APECSS_POW3(Bubble->R0) * Sol[0] / APECSS_POW4(Sol[1]));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble->Interface) +
          3.0 * apecss_viscoelastic_zenervarsigma_ode(Sol, t, Bubble));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_exploldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble->Interface) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) + apecss_viscoelastic_oldroydb1_ode(Sol, t, Bubble) +
          apecss_viscoelastic_oldroydb2_ode(Sol, t, Bubble));
}

// -------------------------------------------------------------------
// VISCOUS PRESSURE
// -------------------------------------------------------------------

APECSS_FLOAT apecss_liquid_pressure_viscous_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble) { return (4.0 * Bubble->Liquid->mu * U / R); }

APECSS_FLOAT apecss_liquid_pressure_viscous_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (4.0 * Bubble->Liquid->mu * U / R + 4.0 * Bubble->Interface->Viscosity * U / APECSS_POW2(R));
}

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_cleanexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (4.0 * Bubble->Liquid->mu * APECSS_POW2(U / R));
}

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_marmottantexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (4.0 * Bubble->Liquid->mu * APECSS_POW2(U / R) + 8.0 * Bubble->Interface->Viscosity * APECSS_POW2(U) / APECSS_POW3(R));
}

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_nonimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble) { return (0.0); }

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_cleanimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble) { return (4.0 * Bubble->Liquid->mu / R); }

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_marmottantimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble)
{
  return (4.0 * (Bubble->Liquid->mu * R + Bubble->Interface->Viscosity) / APECSS_POW2(R));
}

// -------------------------------------------------------------------
// PRESSURE AT INFINITY
// -------------------------------------------------------------------

APECSS_FLOAT apecss_liquid_pressure_infinity_noexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble) { return (Bubble->p0); }

APECSS_FLOAT apecss_liquid_pressure_infinity_sinexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->p0 - Bubble->Excitation->dp * APECSS_SIN(2.0 * APECSS_PI * Bubble->Excitation->f * t));
}

APECSS_FLOAT apecss_liquid_pressurederivative_infinity_noexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble) { return (0.0); }

APECSS_FLOAT apecss_liquid_pressurederivative_infinity_sinexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (-Bubble->Excitation->dp * 2.0 * APECSS_PI * Bubble->Excitation->f * APECSS_COS(2.0 * APECSS_PI * Bubble->Excitation->f * t));
}