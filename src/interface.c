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
// OPTIONS
// -------------------------------------------------------------------
// Functions initializing, processing and handling the options of the
// gas-liquid interface.
// -------------------------------------------------------------------

int apecss_interface_setdefaultoptions(struct APECSS_Interface *Interface)
{
  Interface->sigma = 0.0;
  Interface->LipidCoatingModel = APECSS_LIPIDCOATING_NONE;
  Interface->sigma0 = 0.02;
  Interface->Elasticity = 0.5;
  Interface->Viscosity = 7.5e-9;
  Interface->GompertzC = 0.0;

  Interface->get_surfacetension = apecss_interface_surfacetension_clean;
  Interface->get_pressure_surfacetension = apecss_interface_surfacetensionpressure_clean;
  Interface->get_pressurederivative_surfacetension = apecss_interface_surfacetensionpressurederivative_clean;

  Interface->get_pressure_viscous = apecss_interface_pressure_viscous_clean;
  Interface->get_pressurederivative_viscous_expl = apecss_interface_pressurederivative_viscous_cleanexpl;
  Interface->get_pressurederivative_viscous_impl = apecss_interface_pressurederivative_viscous_cleanimpl;

  return (0);
}

int apecss_interface_readoptions(struct APECSS_Interface *Interface, char *OptionsDir)
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

    if (strncasecmp(option, "interface", 9) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else if (strncasecmp(option2, "surfacetensioncoeff", 19) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Interface->sigma = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "lipidcoatingmodel", 17) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "none", 4) == 0)
            {
              Interface->LipidCoatingModel = APECSS_LIPIDCOATING_NONE;
            }
            else if (strncasecmp(option3, "marmottant", 10) == 0)
            {
              Interface->LipidCoatingModel = APECSS_LIPIDCOATING_MARMOTTANT;
            }
            else if (strncasecmp(option3, "gompertz-marmottant", 19) == 0)
            {
              Interface->LipidCoatingModel = APECSS_LIPIDCOATING_MARMOTTANT + APECSS_LIPIDCOATING_GOMPERTZFUNCTION;
            }
          }
        }
        else if (strncasecmp(option2, "sigmainit", 9) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Interface->sigma0 = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "elasticity", 10) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Interface->Elasticity = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "dilatationalviscosity", 21) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Interface->Viscosity = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of INTERFACE is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    else if (strncasecmp(option, "bubble", 6) == 0 || strncasecmp(option, "gas", 3) == 0 || strncasecmp(option, "liquid", 6) == 0 ||
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

      Interface->GompertzC = 2.0 * Interface->Elasticity * APECSS_E * APECSS_SQRT(1.0 + Interface->sigma * 0.5 / Interface->Elasticity) / Interface->sigma;
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
// Functions defining the surface tension coefficient and its
// derivative, based on the chosen model.
// -------------------------------------------------------------------
// The function for surface tension is chosen in
// apecss_interface_processoptions() and associated with the function
// pointer Interface->get_surfacetension().
// -------------------------------------------------------------------

APECSS_FLOAT apecss_interface_surfacetension_clean(APECSS_FLOAT R, struct APECSS_Bubble *Bubble) { return (Bubble->Interface->sigma); }

APECSS_FLOAT apecss_interface_surfacetension_marmottant(APECSS_FLOAT R, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT sigma;
  if (R < Bubble->Rbuck)
    sigma = 0.0;
  else if (R > Bubble->Rrupt)
    sigma = Bubble->Interface->sigma;
  else
    sigma = Bubble->Interface->Elasticity * (APECSS_POW2(R) / APECSS_POW2(Bubble->Rbuck) - 1.0);

  return (sigma);
}

APECSS_FLOAT apecss_interface_surfacetension_gompertzmarmottant(APECSS_FLOAT R, struct APECSS_Bubble *Bubble)
{
  return (Bubble->Interface->sigma * APECSS_EXP(-Bubble->GompertzB * APECSS_EXP(Bubble->Interface->GompertzC * (1.0 - R / Bubble->Rbuck))));
}

APECSS_FLOAT apecss_interface_surfacetensionderivative_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT dot_sigma;

  if (R < Bubble->Rbuck)
    dot_sigma = 0.0;
  else if (R > Bubble->Rrupt)
    dot_sigma = 0.0;
  else
    dot_sigma = 2.0 * Bubble->Interface->Elasticity * R * U / APECSS_POW2(Bubble->Rbuck);

  return (dot_sigma);
}

APECSS_FLOAT apecss_interface_surfacetensionderivative_gompertzmarmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (apecss_interface_surfacetension_gompertzmarmottant(R, Bubble) * U * Bubble->GompertzB * Bubble->Interface->GompertzC *
          APECSS_EXP(Bubble->Interface->GompertzC * (1.0 - R / Bubble->Rbuck)) / Bubble->Rbuck);
}

// -------------------------------------------------------------------
// SURFACE TENSION PRESSURE
// -------------------------------------------------------------------
// Functions defining the pressure contribution due to surface tension
// -------------------------------------------------------------------
// The functions are chosen in apecss_interface_processoptions() and
// associated with the function pointers:
// - Interface->get_pressure_surfacetension()
// - Interface->get_pressurederivative_surfacetension()
// -------------------------------------------------------------------

APECSS_FLOAT apecss_interface_surfacetensionpressure_clean(APECSS_FLOAT R, struct APECSS_Bubble *Bubble)
{
  return (Bubble->dimensionality * Bubble->Interface->sigma / R);
}

APECSS_FLOAT apecss_interface_surfacetensionpressure_marmottant(APECSS_FLOAT R, struct APECSS_Bubble *Bubble)
{
  return (2.0 * apecss_interface_surfacetension_marmottant(R, Bubble) / R);
}

APECSS_FLOAT apecss_interface_surfacetensionpressure_gompertzmarmottant(APECSS_FLOAT R, struct APECSS_Bubble *Bubble)
{
  return (2.0 * apecss_interface_surfacetension_gompertzmarmottant(R, Bubble) / R);
}

APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (Bubble->dimensionality * Bubble->Interface->sigma * U / APECSS_POW2(R));
}

APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (2.0 * (apecss_interface_surfacetension_marmottant(R, Bubble) * U / APECSS_POW2(R) -
                 apecss_interface_surfacetensionderivative_marmottant(R, U, Bubble) / R));
}

APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_gompertzmarmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble)
{
  return (2.0 * (apecss_interface_surfacetension_gompertzmarmottant(R, Bubble) * U / APECSS_POW2(R) -
                 apecss_interface_surfacetensionderivative_gompertzmarmottant(R, U, Bubble) / R));
}

// -------------------------------------------------------------------
// VISCOUS PRESSURE
// -------------------------------------------------------------------
// Functions defining the pressure contribution due to viscous
// stresses of the interface coating.
// -------------------------------------------------------------------
// The functions are chosen in apecss_interface_processoptions() and
// associated with the function pointers:
// - Interface->get_pressure_viscous()
// - Interface->get_pressurederivative_viscous_expl()
// - Interface->get_pressurederivative_viscous_impl()
// -------------------------------------------------------------------

APECSS_FLOAT apecss_interface_pressure_viscous_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface) { return (0.0); }

APECSS_FLOAT apecss_interface_pressure_viscous_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface)
{
  return (4.0 * Interface->Viscosity * U / APECSS_POW2(R));
}

APECSS_FLOAT apecss_interface_pressurederivative_viscous_cleanexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface) { return (0.0); }

APECSS_FLOAT apecss_interface_pressurederivative_viscous_marmottantexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface)
{
  return (8.0 * Interface->Viscosity * APECSS_POW2(U) / APECSS_POW3(R));
}

APECSS_FLOAT apecss_interface_pressurederivative_viscous_cleanimpl(APECSS_FLOAT R, struct APECSS_Interface *Interface) { return (0.0); }

APECSS_FLOAT apecss_interface_pressurederivative_viscous_marmottantimpl(APECSS_FLOAT R, struct APECSS_Interface *Interface)
{
  return (4.0 * Interface->Viscosity / APECSS_POW2(R));
}