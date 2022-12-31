// This source file is part of APECSS, an open-source software toolbox
// for the computation of pressure-driven bubble dynamics and acoustic
// emissions in spherical symmetry.
//
// Copyright (C) 2022-2023 The APECSS Developers
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

int apecss_gas_readoptions(struct APECSS_Gas *Gas, char *OptionsDir)
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

    if (strncasecmp(option, "gas", 3) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
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
            if (strncasecmp(option3, "ig", 2) == 0)
            {
              Gas->EoS = APECSS_GAS_IG;
            }
            else if (strncasecmp(option3, "hc", 2) == 0)
            {
              Gas->EoS = APECSS_GAS_HC;
            }
            else if (strncasecmp(option3, "nasg", 4) == 0)
            {
              Gas->EoS = APECSS_GAS_NASG;
            }
          }
        }
        else if (strncasecmp(option2, "polytropicexponent", 18) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Gas->Gamma = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencepressure", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Gas->pref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencedensity", 16) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Gas->rhoref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "covolume", 8) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Gas->b = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "taitpressureconst", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Gas->B = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "molecularweight", 15) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Gas->mmol = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "moleculardiameter", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Gas->dmol = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of GAS is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    else if (strncasecmp(option, "bubble", 6) == 0 || strncasecmp(option, "liquid", 6) == 0 || strncasecmp(option, "interface", 9) == 0 ||
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