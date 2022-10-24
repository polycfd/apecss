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

int apecss_gas_setdefaultoptions(struct APECSS_Gas *Gas)
{
  Gas->EoS = APECSS_GAS_IG;
  Gas->Gamma = 1.4;
  Gas->B = 0.0;
  Gas->b = 0.0;
  Gas->dmol = -1.0;
  Gas->mmol = -1.0;
  Gas->pref = 1.0e5;
  Gas->rhoref = 1.2;

  Gas->get_pressure = apecss_gas_pressure_ig;
  Gas->get_pressurederivative = apecss_gas_pressurederivative_ig;

  return (0);
}

int apecss_gas_processoptions(struct APECSS_Gas *Gas)
{
  if (Gas->EoS == APECSS_GAS_IG)
  {
    // Set the function pointers for the pressure
    Gas->get_pressure = apecss_gas_pressure_ig;
    Gas->get_pressurederivative = apecss_gas_pressurederivative_ig;

    // Compute the reference coefficient for the gas
    Gas->Kref = Gas->rhoref / (APECSS_POW(Gas->pref, 1.0 / Gas->Gamma));
  }
  else if (Gas->EoS == APECSS_GAS_NASG)
  {
    // Set the function pointers for the pressure
    Gas->get_pressure = apecss_gas_pressure_nasg;
    Gas->get_pressurederivative = apecss_gas_pressurederivative_nasg;

    if (Gas->mmol > 0.0 && Gas->dmol > 0.0)
    {
      // Compute the co-volume based on the molecular properties
      Gas->b = APECSS_AVOGADRO * 4.0 * APECSS_PI * APECSS_POW3(Gas->dmol) / (6.0 * Gas->mmol);
    }

    // Compute the reference coefficient for the gas
    Gas->Kref = Gas->rhoref / (APECSS_POW(Gas->pref + Gas->B, 1.0 / Gas->Gamma) * (1.0 - Gas->b * Gas->rhoref));
  }
  else if (Gas->EoS == APECSS_GAS_HC)
  {
    // Set the function pointers for the pressure
    Gas->get_pressure = apecss_gas_pressure_hc;
    Gas->get_pressurederivative = apecss_gas_pressurederivative_hc;

    // Compute the reference coefficient for the gas
    Gas->Kref = Gas->rhoref / (APECSS_POW(Gas->pref, 1.0 / Gas->Gamma));
  }
  else
  {
    apecss_erroronscreen(1, "Defined equation of state for the gas is unknown.");
  }

  return (0);
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
  return (Bubble->pG0 * APECSS_POW(Bubble->R0 / Sol[1], 3.0 * Bubble->Gas->Gamma));
}

APECSS_FLOAT apecss_gas_pressure_hc(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble)
{
  return (Bubble->pG0 *
          APECSS_POW((APECSS_POW3(Bubble->R0) - APECSS_POW3(Bubble->r_hc)) / (APECSS_POW3(Sol[1]) - APECSS_POW3(Bubble->r_hc)), Bubble->Gas->Gamma));
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
  return (-3.0 * apecss_gas_pressure_ig(Sol, Bubble) * Bubble->Gas->Gamma * Sol[0] / Sol[1]);
}

APECSS_FLOAT apecss_gas_pressurederivative_hc(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (-3.0 * apecss_gas_pressure_hc(Sol, Bubble) * Bubble->Gas->Gamma * APECSS_POW2(Sol[1]) * Sol[0] / (APECSS_POW3(Sol[1]) - APECSS_POW3(Bubble->r_hc)));
}

APECSS_FLOAT apecss_gas_pressurederivative_nasg(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT rhoG = apecss_gas_density_constmass(Sol[1], Bubble);
  return (apecss_gas_densityderivative_constmass(Sol[1], Sol[0], Bubble) * Bubble->Gas->Gamma * (apecss_gas_pressure_nasg(Sol, Bubble) + Bubble->Gas->B) /
          (rhoG * (1.0 - Bubble->Gas->b * rhoG)));
}