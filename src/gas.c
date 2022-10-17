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

int apecss_gas_setdefaultoptions(struct APECSS_Bubble *Bubble)
{
  Bubble->Gas = (struct APECSS_Gas *) malloc(sizeof(struct APECSS_Gas));

  Bubble->Gas->EOS = APECSS_GAS_IG;
  Bubble->Gas->Gamma = 1.4;
  Bubble->Gas->B = 0.0;
  Bubble->Gas->b = 0.0;
  Bubble->Gas->h = 0.0;
  Bubble->Gas->dmol = -1.0;
  Bubble->Gas->mmol = -1.0;
  Bubble->Gas->pref = -1.0e10;  // Either set by user input or set to ambient pressure
  Bubble->Gas->rhoref = 1.2;

  Bubble->Gas->get_pressure = apecss_gas_pressure_ig;
  Bubble->Gas->get_pressurederivative = apecss_gas_pressurederivative_ig;

  return 0;
}

int apecss_gas_processoptions(struct APECSS_Bubble *Bubble)
{
  // Set the appropriate function pointers associated with the gas equation of state
  if (Bubble->Gas->EOS == APECSS_GAS_IG)
  {
    Bubble->Gas->get_pressure = apecss_gas_pressure_ig;
    Bubble->Gas->get_pressurederivative = apecss_gas_pressurederivative_ig;
  }
  else if (Bubble->Gas->EOS == APECSS_GAS_NASG)
  {
    Bubble->Gas->get_pressure = apecss_gas_pressure_nasg;
    Bubble->Gas->get_pressurederivative = apecss_gas_pressurederivative_nasg;
  }
  else if (Bubble->Gas->EOS == APECSS_GAS_HC)
  {
    Bubble->Gas->get_pressure = apecss_gas_pressure_hc;
    Bubble->Gas->get_pressurederivative = apecss_gas_pressurederivative_hc;
  }
  else
  {
    apecss_erroronscreen(1, "Defined gas model unknown.");
  }

  // Set the reference pressure
  if (Bubble->Gas->pref < -Bubble->Gas->B) Bubble->Gas->pref = Bubble->p0;

  return 0;
}

// -------------------------------------------------------------------
// PROPERTIES
// -------------------------------------------------------------------

APECSS_FLOAT apecss_gas_density_constmass(APECSS_FLOAT R, struct APECSS_Bubble *Bubble) { return (Bubble->rhoG0 * APECSS_POW3(Bubble->R0 / R)); }

APECSS_FLOAT apecss_gas_densityderivative_constmass(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (-3.0 * apecss_gas_density_constmass(R, Bubble) * U / R);
}

// -------------------------------------------------------------------
// PRESSURE
// -------------------------------------------------------------------

APECSS_FLOAT apecss_gas_pressure_ig(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble)
{
  return (Bubble->pG0 * APECSS_POW(apecss_gas_density_constmass(Sol[1], Bubble) / Bubble->rhoG0, Bubble->Gas->Gamma));
}

APECSS_FLOAT apecss_gas_pressure_hc(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble)
{
  return (Bubble->pG0 *
          APECSS_POW((APECSS_POW3(Bubble->R0) - APECSS_POW3(Bubble->Gas->h)) / (APECSS_POW3(Sol[1]) - APECSS_POW3(Bubble->Gas->h)), Bubble->Gas->Gamma));
}

APECSS_FLOAT apecss_gas_pressure_nasg(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT rhoG = apecss_gas_density_constmass(Sol[1], Bubble);
  return (-Bubble->Gas->B +
          (Bubble->pG0 + Bubble->Gas->B) *
              APECSS_POW((rhoG * (1.0 - Bubble->Gas->b * Bubble->rhoG0)) / (Bubble->rhoG0 * (1.0 - Bubble->Gas->b * rhoG)), Bubble->Gas->Gamma));
}

APECSS_FLOAT apecss_gas_pressurederivative_ig(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (apecss_gas_densityderivative_constmass(Sol[1], Sol[0], Bubble) * Bubble->Gas->Gamma * apecss_gas_pressure_ig(Sol, Bubble) /
          apecss_gas_density_constmass(Sol[1], Bubble));
}

APECSS_FLOAT apecss_gas_pressurederivative_hc(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (-3.0 * apecss_gas_pressure_hc(Sol, Bubble) * Bubble->Gas->Gamma * APECSS_POW2(Sol[1]) * Sol[0] / (APECSS_POW3(Sol[1]) - APECSS_POW3(Bubble->Gas->h)));
}

APECSS_FLOAT apecss_gas_pressurederivative_nasg(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT rhoG = apecss_gas_density_constmass(Sol[1], Bubble);
  return (apecss_gas_densityderivative_constmass(Sol[1], Sol[0], Bubble) * Bubble->Gas->Gamma * (apecss_gas_pressure_nasg(Sol, Bubble) + Bubble->Gas->B) /
          (rhoG * (1.0 - Bubble->Gas->b * rhoG)));
}