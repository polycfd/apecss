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

// -------------------------------------------------------------------
// APECSS standalone example for an ultrasound-driven bubble
// including the temperature model of Stricker et al. (2011).
// -------------------------------------------------------------------

#include <time.h>
#include "apecss.h"

// Declaration of the function containing the temperature ODE
APECSS_FLOAT gas_energy_strickersingleode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

// Globally-defined variables (indicated for clarity with leading and trailing underscores)
int _pos_energy_ode_;
APECSS_FLOAT _cv_, _k_, _T0_;

int main(int argc, char **args)
{
  char str[APECSS_STRINGLENGTH_SPRINTF];
  char OptionsDir[APECSS_STRINGLENGTH];

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Initialize the case-dependent simulation parameters
  double tEnd = 0.0;
  double fa = 0.0;
  double pa = 0.0;
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  apecss_infoscreen();

  /* Read commandline options */
  sprintf(OptionsDir, "./run.apecss");
  int i = 1;  // First argument is the call to the executable
  while (i < argc)
  {
    if (strcmp("-options", args[i]) == 0)
    {
      sprintf(OptionsDir, "%s", args[i + 1]);
      i += 2;
    }
    else if (strcmp("-tend", args[i]) == 0)
    {
      sscanf(args[i + 1], "%le", &tEnd);
      i += 2;
    }
    else if (strcmp("-freq", args[i]) == 0)
    {
      sscanf(args[i + 1], "%le", &fa);
      i += 2;
    }
    else if (strcmp("-amp", args[i]) == 0)
    {
      sscanf(args[i + 1], "%le", &pa);
      i += 2;
    }
    else
    {
      char str[APECSS_STRINGLENGTH_SPRINTF];
      sprintf(str, "Unknown command line options: %s", args[i]);
      apecss_erroronscreen(1, str);
      ++i;
    }
  }

  /* Allocate and initialize Bubble structure */
  struct APECSS_Bubble *Bubble = (struct APECSS_Bubble *) malloc(sizeof(struct APECSS_Bubble));
  apecss_bubble_initializestruct(Bubble);

  /* Set default options for the bubble and the fluids */
  apecss_options_setdefault(Bubble);

  /* Read the options file */
  apecss_options_readfile(Bubble, OptionsDir);

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Set the case-dependent simulation parameters
  Bubble->tStart = 0.0;
  Bubble->tEnd = (APECSS_FLOAT) tEnd;
  Bubble->dt = APECSS_MIN(1.0e-7, Bubble->tEnd - Bubble->tStart);  // Initial time-step
  Bubble->Excitation = (struct APECSS_Excitation *) malloc(sizeof(struct APECSS_Excitation));
  Bubble->Excitation->type = APECSS_EXCITATION_SIN;
  Bubble->Excitation->f = (APECSS_FLOAT) fa;
  Bubble->Excitation->dp = (APECSS_FLOAT) pa;

  // Set the parameters and functions for the gas energy model
  _cv_ = 720.0;  // Isochoric heat capacity of the gas
  _k_ = 0.025;  // Thermal conductivity
  _T0_ = 293.15;  // Ambient temperature

  // (Optional) The solutions of the additional ODEs may be stored with the RP solution
  if (Bubble->Results != NULL && Bubble->Results->RayleighPlesset != NULL)
  {
    Bubble->Results->RayleighPlesset->UserODEsName = malloc(sizeof(char *));
    Bubble->Results->RayleighPlesset->UserODEsSol = malloc(sizeof(APECSS_FLOAT *));

    Bubble->Results->RayleighPlesset->UserODEsName[0] = malloc(APECSS_STRINGLENGTH * sizeof(char));
    sprintf(Bubble->Results->RayleighPlesset->UserODEsName[0], "TG");

    // The solutions of the first nUserODEs additional ODEs defined by the user are written to file
    Bubble->Results->RayleighPlesset->nUserODEs = 1;
  }

  // Set the total number of additional user-defined ODEs
  Bubble->nUserODEs = 1;

  // Bubble.nODEs represents at this point the total count of ODEs that is used for allocating the corresponding arrays
  Bubble->nODEs += Bubble->nUserODEs;

  // Note that Bubble.nODEs is reset after initial allocation in apecss_options_process() and,
  // subsequently, counted upwards again as the functions of the ODEs are defined.

  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Process all options */
  apecss_options_process(Bubble);

  /* Initialize the bubble based on the selected options */
  apecss_bubble_initialize(Bubble);

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Add the temperature ODE to the set of solved ODEs
  _pos_energy_ode_ = Bubble->nODEs;  // Simplfies finding the corresponding solution
  Bubble->ode[_pos_energy_ode_] = gas_energy_strickersingleode;
  Bubble->ODEsSol[_pos_energy_ode_] = _T0_;
  Bubble->nODEs++;
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Solve the bubble dynamics */
  clock_t starttimebubble = clock();
  apecss_bubble_solve(Bubble);

  sprintf(str, "Solver concluded %i time-steps and %i sub-iterations in %.3f s.", Bubble->dtNumber, Bubble->nSubIter,
          (double) (clock() - starttimebubble) / CLOCKS_PER_SEC);
  apecss_writeonscreen(str);

  /* Write out all desired results */
  apecss_results_rayleighplesset_write(Bubble);
  apecss_results_emissionsspace_write(Bubble);
  apecss_results_emissionsnodespecific_write(Bubble);
  apecss_results_emissionsnodeminmax_write(Bubble);

  /* Make sure all allocated memory is freed */
  apecss_bubble_freestruct(Bubble);

  return (0);
}

APECSS_FLOAT gas_energy_strickersingleode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT U = Sol[0];
  APECSS_FLOAT R = Sol[1];
  APECSS_FLOAT TG = Sol[_pos_energy_ode_];

  // see Stricker, Prosperetti & Lohse, J. Acoust. Soc. Am. 130 (2011), 3243
  APECSS_FLOAT lth = APECSS_MIN(APECSS_SQRT(R * _k_ / (APECSS_ABS(U) * apecss_gas_density_constmass(R, Bubble) * Bubble->Gas->Gamma * _cv_ + APECSS_SMALL)),
                                R / APECSS_PI);  // Thermal lengthscale
  APECSS_FLOAT Qcond = 4.0 * APECSS_PI * APECSS_POW2(R) * _k_ * (_T0_ - TG) / lth;  // Heat due to conduction between liquid and gas
  APECSS_FLOAT Wvol = Bubble->Gas->get_pressure(Sol, Bubble) * 4.0 * APECSS_PI * APECSS_POW2(R) * U;  // Volume work

  return ((Qcond - Wvol) / (_cv_ * Bubble->rhoG0 * 4.0 * APECSS_ONETHIRD * APECSS_PI * APECSS_POW3(Bubble->R0)));
}