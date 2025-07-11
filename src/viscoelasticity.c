// This source file is part of APECSS, an open-source software toolbox
// for the computation of pressure-driven bubble dynamics and acoustic
// emissions in spherical symmetry.
//
// Copyright (C) 2022-2025 The APECSS Developers
//
// The APECSS Developers are listed in the README.md file available in
// the GitHub repository at https://github.com/polycfd/apecss.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "apecss.h"

// -------------------------------------------------------------------
// ZENER
// -------------------------------------------------------------------
// Functions with additional ODEs defining the stresses in a Zener
// solid.
// -------------------------------------------------------------------
// Functions are hooked up to the function pointers Bubble->ode[2] and
// Bubble->ode[3] in apecss_bubble_processoptions().
// -------------------------------------------------------------------

APECSS_FLOAT apecss_viscoelastic_zenertaurr_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  /** ODE for the radial stress at the bubble wall, according to Eq. (7) of Hua & Johnsen, Physics of Fluids 25 (2013), 083101 **/
  APECSS_FLOAT U = Sol[0];
  APECSS_FLOAT R = Sol[1];
  APECSS_FLOAT taurr = Sol[3];
  APECSS_FLOAT S = 4.0 * (APECSS_ONETHIRD * Bubble->Liquid->G * (1.0 - APECSS_POW3(Bubble->R0) / APECSS_POW3(Bubble->R)) + Bubble->Liquid->mu * U / R);

  return ((-S - taurr) / (Bubble->Liquid->lambda + Bubble->dt));
}

APECSS_FLOAT apecss_viscoelastic_zenervarsigma_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  /** ODE for the auxilliary variable associated with the stress at the bubble wall, according to Eq. (6) of Hua & Johnsen, Physics of Fluids 25 (2013), 083101 **/
  APECSS_FLOAT U = Sol[0];
  APECSS_FLOAT R = Sol[1];
  APECSS_FLOAT varsigma = Sol[2];
  APECSS_FLOAT taurr = Sol[3];
  APECSS_FLOAT S = 4.0 * (APECSS_ONETHIRD * Bubble->Liquid->G * (1.0 - APECSS_POW3(Bubble->R0) / APECSS_POW3(Bubble->R)) + Bubble->Liquid->mu * U / R);

  return ((-APECSS_ONETHIRD * S - Bubble->Liquid->lambda * taurr * U / R - varsigma) / (Bubble->Liquid->lambda + Bubble->dt));
}

// -------------------------------------------------------------------
// UCM/OLDROYD-B
// -------------------------------------------------------------------
// Functions with additional ODEs defining the stresses in a upper-
// convected Maxwell or Oldroyd-B fluid.
// -------------------------------------------------------------------
// Functions are hooked up to the function pointers Bubble->ode[2] and
// Bubble->ode[3] in apecss_bubble_processoptions().
// -------------------------------------------------------------------

APECSS_FLOAT apecss_viscoelastic_oldroydb1_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  /** ODE for auxilliary variable 1 associated with the stress at the bubble wall, according to Eq. (6) of Jimenez-Fernandez & Crespo, Ultrasonics 43 (2005), 643–651 **/
  APECSS_FLOAT U = Sol[0];
  APECSS_FLOAT R = Sol[1];
  APECSS_FLOAT S1 = Sol[2];

  return ((-(4.0 * Bubble->Liquid->lambda * S1 + 2.0 * Bubble->Liquid->eta) * U / R - S1) / (Bubble->Liquid->lambda + Bubble->dt));
}

APECSS_FLOAT apecss_viscoelastic_oldroydb2_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  /** ODE for auxilliary variable 2 associated with the stress at the bubble wall, according to Eq. (7) of Jimenez-Fernandez & Crespo, Ultrasonics 43 (2005), 643–651 **/
  APECSS_FLOAT U = Sol[0];
  APECSS_FLOAT R = Sol[1];
  APECSS_FLOAT S2 = Sol[3];

  return ((-(Bubble->Liquid->lambda * S2 + 2.0 * Bubble->Liquid->eta) * U / R - S2) / (Bubble->Liquid->lambda + Bubble->dt));
}