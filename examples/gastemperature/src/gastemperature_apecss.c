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
// int _pos_energy_ode_;
// APECSS_FLOAT _cv_, _k_, _T0_;

struct AdditionalGasProperties
{
  APECSS_FLOAT cv, k, T0;
};

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
  sprintf(OptionsDir, "./run.apecss");  // This is the default
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

  /* Set default options and read the options for the bubble */
  apecss_bubble_setdefaultoptions(Bubble);
  apecss_bubble_readoptions(Bubble, OptionsDir);

  /* Allocate the structures for the fluid properties and ODE solver parameters */
  struct APECSS_Gas *Gas = (struct APECSS_Gas *) malloc(sizeof(struct APECSS_Gas));
  struct APECSS_Liquid *Liquid = (struct APECSS_Liquid *) malloc(sizeof(struct APECSS_Liquid));
  struct APECSS_Interface *Interface = (struct APECSS_Interface *) malloc(sizeof(struct APECSS_Interface));
  struct APECSS_NumericsODE *NumericsODE = (struct APECSS_NumericsODE *) malloc(sizeof(struct APECSS_NumericsODE));

  /* Set the default options for the fluid properties and solver parameters  */
  apecss_gas_setdefaultoptions(Gas);
  apecss_liquid_setdefaultoptions(Liquid);
  apecss_interface_setdefaultoptions(Interface);
  apecss_odesolver_setdefaultoptions(NumericsODE);

  /* Read the options file for the fluid properties and solver parameters  */
  apecss_gas_readoptions(Gas, OptionsDir);
  apecss_liquid_readoptions(Liquid, OptionsDir);
  apecss_interface_readoptions(Interface, OptionsDir);
  apecss_odesolver_readoptions(NumericsODE, OptionsDir);

  /* Associate the bubble with the relevant fluid properties and solver parameters */
  Bubble->Gas = Gas;
  Bubble->Liquid = Liquid;
  Bubble->Interface = Interface;
  Bubble->NumericsODE = NumericsODE;

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Set the case-dependent simulation parameters
  Bubble->tStart = 0.0;
  Bubble->tEnd = (APECSS_FLOAT) tEnd;
  Bubble->dt = APECSS_MIN(1.0e-7, Bubble->tEnd - Bubble->tStart);  // Initial time-step

  struct APECSS_Excitation *Excitation = (struct APECSS_Excitation *) malloc(sizeof(struct APECSS_Excitation));
  Excitation->type = APECSS_EXCITATION_SIN;
  Excitation->f = (APECSS_FLOAT) fa;
  Excitation->dp = (APECSS_FLOAT) pa;

  Bubble->Excitation = Excitation;

  // Allocate and set the structure holding the additional properties of the gas
  struct AdditionalGasProperties *gas_data = (struct AdditionalGasProperties *) malloc(sizeof(struct AdditionalGasProperties));
  gas_data->cv = 720.0;  // Isochoric heat capacity of the gas
  gas_data->k = 0.025;  // Thermal conductivity
  gas_data->T0 = 293.15;  // Ambient temperature
  Gas->user_data = gas_data;  // Hook case-dependent data structure to the void data pointer

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
  apecss_gas_processoptions(Gas);
  apecss_liquid_processoptions(Liquid);
  apecss_interface_processoptions(Interface);
  apecss_odesolver_processoptions(NumericsODE);
  apecss_bubble_processoptions(Bubble);

  /* Initialize the bubble based on the selected options */
  apecss_bubble_initialize(Bubble);

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Add the temperature ODE to the set of solved ODEs
  int pos_energy_ode = Bubble->nODEs;  // Simplfies finding the corresponding solution
  Bubble->ode[pos_energy_ode] = gas_energy_strickersingleode;
  Bubble->ODEsSol[pos_energy_ode] = gas_data->T0;
  Bubble->user_data = &pos_energy_ode;  // Hook the position of the additional ODE to the void data pointer
  Bubble->nODEs++;
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Solve the bubble dynamics */
  clock_t starttimebubble = clock();
  apecss_bubble_solver_initialize(Bubble);
  apecss_bubble_solver_run(Bubble->tEnd, Bubble);
  apecss_bubble_solver_finalize(Bubble);

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

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Freeing the structure holding the additional properties of the gas
  free(gas_data);
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  free(Bubble);
  free(Gas);
  free(Liquid);
  free(Interface);
  free(NumericsODE);
  free(Excitation);

  return (0);
}

APECSS_FLOAT gas_energy_strickersingleode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  int *pos_energy_ode = Bubble->user_data;
  struct AdditionalGasProperties *gas_data = Bubble->Gas->user_data;

  APECSS_FLOAT U = Sol[0];
  APECSS_FLOAT R = Sol[1];
  APECSS_FLOAT TG = Sol[*pos_energy_ode];

  // see Stricker, Prosperetti & Lohse, J. Acoust. Soc. Am. 130 (2011), 3243
  APECSS_FLOAT lth =
      APECSS_MIN(APECSS_SQRT(R * gas_data->k / (APECSS_ABS(U) * apecss_gas_density_constmass(R, Bubble) * Bubble->Gas->Gamma * gas_data->cv + APECSS_SMALL)),
                 R / APECSS_PI);  // Thermal lengthscale
  APECSS_FLOAT Qcond = 4.0 * APECSS_PI * APECSS_POW2(R) * gas_data->k * (gas_data->T0 - TG) / lth;  // Heat due to conduction between liquid and gas
  APECSS_FLOAT Wvol = Bubble->Gas->get_pressure(Sol, Bubble) * 4.0 * APECSS_PI * APECSS_POW2(R) * U;  // Volume work

  return ((Qcond - Wvol) / (gas_data->cv * Bubble->rhoG0 * 4.0 * APECSS_ONETHIRD * APECSS_PI * APECSS_POW3(Bubble->R0)));
}