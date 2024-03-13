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

// -------------------------------------------------------------------
// APECSS standalone example for an acoustic emitter.
// -------------------------------------------------------------------

#include <time.h>
#include "apecss.h"

APECSS_FLOAT emitter_liquid_pressure_emitterwall(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
int planaremitter_emissions_updatelinkedlist(struct APECSS_Bubble *Bubble);

APECSS_FLOAT fa = 0.0;  // Acoustic frequency
APECSS_FLOAT dpa = 0.0;  // Acoustic pressure amplitude

int main(int argc, char **args)
{
  char str[APECSS_STRINGLENGTH_SPRINTF];
  char OptionsDir[APECSS_STRINGLENGTH];

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Initialize the case-dependent simulation parameters
  APECSS_FLOAT tend = 0.0;  // End time of the simulation
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
    else if (strcmp("-fa", args[i]) == 0)
    {
      sscanf(args[i + 1], "%le", &fa);
      i += 2;
    }
    else if (strcmp("-dpa", args[i]) == 0)
    {
      sscanf(args[i + 1], "%le", &dpa);
      i += 2;
    }
    else if (strcmp("-tend", args[i]) == 0)
    {
      sscanf(args[i + 1], "%le", &tend);
      i += 2;
    }
    else if (strcmp("-h", args[i]) == 0)
    {
      apecss_helpscreen();
    }
    else
    {
      char str[APECSS_STRINGLENGTH_SPRINTF];
      sprintf(str, "Unknown command line options: %s", args[i]);
      apecss_erroronscreen(1, str);
      ++i;
    }
  }

  // Desired time step (200 time steps per period)
  APECSS_FLOAT dt_desired = 1.0 / (200.0 * fa);

  /* Allocate and initialize Bubble structure */
  struct APECSS_Bubble *Bubble = (struct APECSS_Bubble *) malloc(sizeof(struct APECSS_Bubble));
  apecss_bubble_initializestruct(Bubble);
  apecss_bubble_setdefaultoptions(Bubble);
  apecss_bubble_readoptions(Bubble, OptionsDir);

  struct APECSS_Liquid *Liquid = (struct APECSS_Liquid *) malloc(sizeof(struct APECSS_Liquid));
  apecss_liquid_setdefaultoptions(Liquid);
  apecss_liquid_readoptions(Liquid, OptionsDir);
  Bubble->Liquid = Liquid;

  apecss_bubble_processoptions(Bubble);
  apecss_liquid_processoptions(Liquid);

  Bubble->Liquid->get_pressure_bubblewall = emitter_liquid_pressure_emitterwall;

  // To only simulate a part of a very long wave train emitted by the planar emitter
  if (0 == Bubble->dimensionality) Bubble->emissions_update = planaremitter_emissions_updatelinkedlist;

  Bubble->t = Bubble->tStart;
  Bubble->dt = dt_desired;
  Bubble->R = Bubble->R0;
  Bubble->U = 0.0;

  Bubble->dtNumber = 0;
  Bubble->emissions_initialize(Bubble);

  clock_t starttimebubble = clock();

  // Solve the acoustic emissions
  while (Bubble->t < tend - APECSS_SMALL)
  {
    APECSS_FLOAT next_event_time = apecss_results_emissionstime_checktime(tend, Bubble);
    Bubble->dt = APECSS_MIN(dt_desired, next_event_time - Bubble->t);

    ++(Bubble->dtNumber);
    Bubble->t += Bubble->dt;

    APECSS_FLOAT p = Liquid->get_pressure_bubblewall(Bubble->ODEsSol, Bubble->t, Bubble);
    APECSS_FLOAT rho = Liquid->get_density(p, Liquid);
    APECSS_FLOAT c = Liquid->get_soundspeed(p, rho, Liquid);

    APECSS_FLOAT omega = 2.0 * APECSS_PI * fa;
    Bubble->U = dpa * APECSS_COS(omega * Bubble->t - 0.5 * APECSS_PI) / (rho * c);
    Bubble->R = Bubble->R0 - dpa * APECSS_SIN(omega * Bubble->t - 0.5 * APECSS_PI) / ((omega + APECSS_SMALL) * rho * c);

    // Acoustic emissions (if applicable)
    Bubble->emissions_update(Bubble);

    // Store results (if applicable)
    Bubble->results_emissionsspace_store(Bubble);

    // Write one-off results (if applicable)
    Bubble->results_emissionstime_write(Bubble);
  }

  apecss_bubble_solver_finalize(Bubble);

  sprintf(str, "Solver concluded %i time-steps in %.3f s.", Bubble->dtNumber, (double) (clock() - starttimebubble) / CLOCKS_PER_SEC);
  apecss_writeonscreen(str);

  /* Write out all desired results */
  apecss_results_emissionsspace_write(Bubble, APECSS_RESULTS_WRITE);
  apecss_results_emissionsnodespecific_write(Bubble);
  apecss_results_emissionsnodeminmax_write(Bubble);

  /* Make sure all allocated memory is freed */
  apecss_bubble_freestruct(Bubble);

  free(Bubble);
  free(Liquid);

  return (0);
}

APECSS_FLOAT emitter_liquid_pressure_emitterwall(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->p0 + dpa * APECSS_COS(2.0 * APECSS_PI * fa * Bubble->t - 0.5 * APECSS_PI));
}

int planaremitter_emissions_updatelinkedlist(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Emissions->nNodes) Bubble->Emissions->advance(Bubble);

  // Only the waves emitted during a short time interval (10 excitation periods) are tracked.
  if (Bubble->t * fa < 10.0) apecss_emissions_addnode(Bubble);
  if (Bubble->Emissions->LastNode->r > Bubble->Emissions->CutOffDistance) apecss_emissions_removenode(Bubble);

  return (0);
}
