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

APECSS_FLOAT apecss_rp_bubbleradius_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble) { return (Sol[0]); }

APECSS_FLOAT apecss_rp_rayleighplessetvelocity_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (((Bubble->Liquid->get_pressure_bubblewall(Sol, t, Bubble) - Bubble->Liquid->get_pressure_infinity(t, Bubble)) / Bubble->Liquid->rhoref -
           1.5 * APECSS_POW2(Sol[0])) /
          Sol[1]);
}

APECSS_FLOAT apecss_rp_rayleighplessetacousticrationvelocity_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT inv_c = 1.0 / Bubble->Liquid->cref;
  APECSS_FLOAT inv_rho = 1.0 / Bubble->Liquid->rhoref;
  return (((Bubble->Liquid->get_pressure_bubblewall(Sol, t, Bubble) - Bubble->Liquid->get_pressure_infinity(t, Bubble)) * inv_rho - 1.5 * APECSS_POW2(Sol[0])) /
              Sol[1] +
          Bubble->Gas->get_pressurederivative(Sol, t, Bubble) * inv_rho * inv_c);
}

APECSS_FLOAT apecss_rp_kellermiksisvelocity_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  /** Formulation according to Prosperetti & Lezzi 1986 **/
  APECSS_FLOAT inv_c = 1.0 / Bubble->Liquid->cref;
  APECSS_FLOAT inv_rho = 1.0 / Bubble->Liquid->rhoref;
  APECSS_FLOAT rhs =
      ((1.0 + Sol[0] * inv_c) * (Bubble->Liquid->get_pressure_bubblewall(Sol, t, Bubble) - Bubble->Liquid->get_pressure_infinity(t, Bubble)) * inv_rho -
       1.5 * (1.0 - (Sol[0] * APECSS_ONETHIRD * inv_c)) * APECSS_POW2(Sol[0])) /
          Sol[1] +
      Bubble->Liquid->get_pressurederivative_bubblewall_expl(Sol, t, Bubble) * inv_c * inv_rho;
  return (rhs / (1.0 - Sol[0] * inv_c + Bubble->Liquid->get_pressurederivative_viscous_impl(Sol[1], Bubble) * inv_c * inv_rho));
}

APECSS_FLOAT apecss_rp_gilmorevelocity_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT pL = Bubble->Liquid->get_pressure_bubblewall(Sol, t, Bubble);
  APECSS_FLOAT pInf = Bubble->Liquid->get_pressure_infinity(t, Bubble);
  APECSS_FLOAT rhoL = apecss_liquid_density_nasg(pL, Bubble->Liquid);
  APECSS_FLOAT rhoInf = apecss_liquid_density_nasg(pInf, Bubble->Liquid);
  APECSS_FLOAT H = apecss_liquid_enthalpy_nasg(pL, rhoL, Bubble->Liquid) - apecss_liquid_enthalpy_nasg(pInf, rhoInf, Bubble->Liquid);
  APECSS_FLOAT dot_Hexpl =
      Bubble->Liquid->get_pressurederivative_bubblewall_expl(Sol, t, Bubble) / rhoL - Bubble->Liquid->get_pressurederivative_infinity(t, Bubble) / rhoInf;
  APECSS_FLOAT inv_cL = 1.0 / apecss_liquid_soundspeed_nasg(pL, rhoL, Bubble->Liquid);
  APECSS_FLOAT GilmoreCoeffB = 1.0 + Bubble->Liquid->get_pressurederivative_viscous_impl(Sol[1], Bubble) * inv_cL / rhoL;
  return ((((1.0 + Sol[0] * inv_cL) * H - 1.5 * (1.0 - Sol[0] * APECSS_ONETHIRD * inv_cL) * APECSS_POW2(Sol[0])) / ((1.0 - Sol[0] * inv_cL) * Sol[1]) +
           dot_Hexpl * inv_cL) /
          GilmoreCoeffB);
}