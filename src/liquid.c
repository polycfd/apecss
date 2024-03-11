// This source file is part of APECSS, an open-source software toolbox
// for the computation of pressure-driven bubble dynamics and acoustic
// emissions in spherical symmetry.
//
// Copyright (C) 2022-2024 The APECSS Developers
//
// The APECSS Developers are listed in the README.md file available in
// the GitHub repository at https://github.com/polycfd/apecss.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "apecss.h"

// -------------------------------------------------------------------
// OPTIONS
// -------------------------------------------------------------------
// Functions initializing, processing and handling the options of the
// liquid phase.
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

int apecss_liquid_readoptions(struct APECSS_Liquid *Liquid, char *OptionsDir)
{
  int l = 0;
  int line = 0;
  FILE *OptionsFile;
  char str[APECSS_STRINGLENGTH_SPRINTF];
  char option[APECSS_STRINGLENGTH], option2[APECSS_STRINGLENGTH], option3[APECSS_STRINGLENGTH];
  int StatusFile = 1;
  int StatusSection = 1;

  if ((OptionsFile = fopen(OptionsDir, "r")) == (FILE *) NULL)
  {
    sprintf(str, "File %s cannot be opened for reading.\n", OptionsDir);
    apecss_erroronscreen(1, str);
  }

  while ((l = apecss_readoneoption(OptionsFile, option)) != EOF && StatusFile == 1)
  {
    line += l;

    if (strncasecmp(option, "liquid", 6) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else if (strncasecmp(option2, "liquidtype", 10) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "newtonian", 9) == 0)
            {
              Liquid->Type = APECSS_LIQUID_NEWTONIAN;
            }
            else if (strncasecmp(option3, "kelvinvoigt", 11) == 0)
            {
              Liquid->Type = APECSS_LIQUID_KELVINVOIGT;
            }
            else if (strncasecmp(option3, "zener", 5) == 0)
            {
              Liquid->Type = APECSS_LIQUID_ZENER;
            }
            else if (strncasecmp(option3, "oldroydb", 8) == 0)
            {
              Liquid->Type = APECSS_LIQUID_OLDROYDB;
            }
          }
        }
        else if (strncasecmp(option2, "eos", 3) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "tait", 4) == 0)
            {
              Liquid->EoS = APECSS_LIQUID_TAIT;
            }
            else if (strncasecmp(option3, "nasg", 4) == 0)
            {
              Liquid->EoS = APECSS_LIQUID_NASG;
            }
          }
        }
        else if (strncasecmp(option2, "polytropicexponent", 18) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->Gamma = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencepressure", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->pref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencedensity", 16) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->rhoref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencesoundspeed", 19) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->cref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "covolume", 8) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->b = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "taitpressureconst", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->B = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "viscosity", 9) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->mu = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "shearmodulus", 12) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->G = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "polymerviscosity", 16) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->eta = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "relaxationtime", 14) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Liquid->lambda = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of LIQUID is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    else if (strncasecmp(option, "bubble", 6) == 0 || strncasecmp(option, "gas", 3) == 0 || strncasecmp(option, "interface", 9) == 0 ||
             strncasecmp(option, "results", 7) == 0 || strncasecmp(option, "odesolver", 9) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else
        {
          // Nothing to be done here.
        }
      }
    }
    else
    {
      sprintf(str, "An unknown Section is given: %s, line %i", option, line);
      apecss_erroronscreen(1, str);
    }
  }

  fclose(OptionsFile);

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
// Functions defining the thermodynamic properties of the liquid phase
// -------------------------------------------------------------------
// The functions are chosen in apecss_liquid_processoptions() and
// associated with the function pointers:
// - Liquid->get_density()
// - Liquid->get_soundspeed()
// - Liquid->get_enthalpy()
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
// Functions defining the bubble wall pressure and its derivatives,
// dependent on the rheology of the liquid, the viscous stresses and
// surface tension contribution.
// -------------------------------------------------------------------
// The functions are chosen in apecss_liquid_processoptions() and
// associated with the function pointers:
// - Gas->get_pressure_bubblewall()
// - Gas->get_pressurederivative_bubblewall_expl()
// -------------------------------------------------------------------

APECSS_FLOAT apecss_liquid_pressure_bubblewall(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble) -
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble) - Bubble->Interface->get_pressure_viscous(Sol[1], Sol[0], Bubble->Interface));
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_kelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble) -
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble) - Bubble->Interface->get_pressure_viscous(Sol[1], Sol[0], Bubble->Interface) -
          4.0 * APECSS_ONETHIRD * Bubble->Liquid->G * (1.0 - APECSS_POW3(Bubble->R0) / APECSS_POW3(Sol[1])));
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble) + 3.0 * Sol[2] -
          Bubble->Interface->get_pressure_viscous(Sol[1], Sol[0], Bubble->Interface));
}

APECSS_FLOAT apecss_liquid_pressure_bubblewall_oldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressure(Sol, Bubble) - Bubble->Interface->get_pressure_surfacetension(Sol[1], Bubble) -
          Bubble->Liquid->get_pressure_viscous(Sol[1], Sol[0], Bubble) - Bubble->Interface->get_pressure_viscous(Sol[1], Sol[0], Bubble->Interface) + Sol[2] +
          Sol[3]);
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_expl(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) +
          Bubble->Interface->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble->Interface));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_explkelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) +
          Bubble->Interface->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble->Interface) -
          4.0 * Bubble->Liquid->G * APECSS_POW3(Bubble->R0) * Sol[0] / APECSS_POW4(Sol[1]));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble) +
          3.0 * apecss_viscoelastic_zenervarsigma_ode(Sol, t, Bubble) +
          Bubble->Interface->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble->Interface));
}

APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_exploldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Gas->get_pressurederivative(Sol, t, Bubble) + Bubble->Interface->get_pressurederivative_surfacetension(Sol[1], Sol[0], Bubble) +
          Bubble->Liquid->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble) +
          Bubble->Interface->get_pressurederivative_viscous_expl(Sol[1], Sol[0], Bubble->Interface) + apecss_viscoelastic_oldroydb1_ode(Sol, t, Bubble) +
          apecss_viscoelastic_oldroydb2_ode(Sol, t, Bubble));
}

// -------------------------------------------------------------------
// VISCOUS PRESSURE
// -------------------------------------------------------------------
// Functions defining the pressure contribution due to viscous
// stresses.
// -------------------------------------------------------------------
// The functions are chosen in apecss_liquid_processoptions() and
// associated with the function pointers:
// - Liquid->get_pressure_viscous()
// - Liquid->get_pressurederivative_viscous_expl()
// - Liquid->get_pressurederivative_viscous_impl()
// -------------------------------------------------------------------

APECSS_FLOAT apecss_liquid_pressure_viscous(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (2.0 * Bubble->dimensionality * Bubble->Liquid->mu * U / R);
}

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_expl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (2.0 * Bubble->dimensionality * Bubble->Liquid->mu * APECSS_POW2(U / R));
}

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_impl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble)
{
  return (2.0 * Bubble->dimensionality * Bubble->Liquid->mu / R);
}

APECSS_FLOAT apecss_liquid_pressurederivative_viscous_nonimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble) { return (0.0); }
