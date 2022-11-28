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
// APECSS standalone example for an ultrasound-driven bubble.
// -------------------------------------------------------------------

#include <time.h>
#include "apecss.h"

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
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Process all options */
  apecss_gas_processoptions(Gas);
  apecss_liquid_processoptions(Liquid);
  apecss_interface_processoptions(Interface);
  apecss_odesolver_processoptions(NumericsODE);
  apecss_bubble_processoptions(Bubble);

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Set the function pointers for the progress bar
  if (Bubble->Emissions != NULL)
  {
    // Optional! Displays a progress bar during the simulation, here only
    // if emissions are calculated.
    Bubble->progress_initial = apecss_bubble_solver_progress_initialscreen;
    Bubble->progress_update = apecss_bubble_solver_progress_updatescreen;
    Bubble->progress_final = apecss_bubble_solver_progress_finalscreen;
  }
  else
  {
    Bubble->progress_initial = apecss_bubble_solver_progress_initialnone;
    Bubble->progress_update = apecss_bubble_solver_progress_updatenone;
    Bubble->progress_final = apecss_bubble_solver_progress_finalnone;
  }
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Initialize the bubble based on the selected options */
  apecss_bubble_initialize(Bubble);

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

  free(Bubble);
  free(Gas);
  free(Liquid);
  free(Interface);
  free(NumericsODE);
  free(Excitation);

  return (0);
}