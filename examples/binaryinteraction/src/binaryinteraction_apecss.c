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
// APECSS standalone example of acoustically-interacting microbubbles,
// based on Jiang et al., Ultrasonics Sonochemistry 34 (2017), 90-97.
// -------------------------------------------------------------------

#include <time.h>
#include "apecss.h"

// Declaration of additional case-dependent functions
APECSS_FLOAT interaction_bubble_pressure_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

// Declaration of the structure holding the interaction variables of each bubble
struct Interaction
{
  APECSS_FLOAT dp_neighbor;
};

int main(int argc, char **args)
{
  char OptionsDir[APECSS_STRINGLENGTH];

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Set the case-dependent simulation parameters
  const int nBubbles = 2;  // Number of bubbles
  APECSS_FLOAT bubble_bubble_dist = 200.0e-6;  // Bubble-bubble distance

  // Interbubble time-step, defining the frequency with which the neighbor influence is updated
  APECSS_FLOAT dt_interbubble = 1.0e-8;

  // Initialize the simulation parameters given by the execution command
  double tEnd = 0.0;
  double fa = 0.0;
  double pa = 0.0;
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  apecss_infoscreen();

  /* Read commandline options */
  sprintf(OptionsDir, "./run.apecss");  // This is the default
  int j = 1;  // First argument is the call to the executable
  while (j < argc)
  {
    if (strcmp("-options", args[j]) == 0)
    {
      sprintf(OptionsDir, "%s", args[j + 1]);
      j += 2;
    }
    else if (strcmp("-tend", args[j]) == 0)
    {
      sscanf(args[j + 1], "%le", &tEnd);
      j += 2;
    }
    else if (strcmp("-freq", args[j]) == 0)
    {
      sscanf(args[j + 1], "%le", &fa);
      j += 2;
    }
    else if (strcmp("-amp", args[j]) == 0)
    {
      sscanf(args[j + 1], "%le", &pa);
      j += 2;
    }
    else if (strcmp("-h", args[j]) == 0)
    {
      apecss_helpscreen();
    }
    else
    {
      char str[APECSS_STRINGLENGTH_SPRINTF];
      sprintf(str, "Unknown command line options: %s", args[j]);
      apecss_erroronscreen(1, str);
      ++j;
    }
  }

  /* Allocate and initialize Bubble structure */
  struct APECSS_Bubble *Bubbles[nBubbles];
  for (register int i = 0; i < nBubbles; i++) Bubbles[i] = (struct APECSS_Bubble *) malloc(nBubbles * sizeof(struct APECSS_Bubble));
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_initializestruct(Bubbles[i]);

  /* Set default options and read the options for the bubble */
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_setdefaultoptions(Bubbles[i]);
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_readoptions(Bubbles[i], OptionsDir);

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Allocate and set bubble-associated variables for the interaction
  for (register int i = 0; i < nBubbles; i++)
  {
    struct Interaction *interaction_data = (struct Interaction *) malloc(sizeof(struct Interaction));
    interaction_data->dp_neighbor = 0.0;  // Pressure excerted by the neighbor bubbles

    // Hook interaction-structure to the void data pointer
    Bubbles[i]->user_data = interaction_data;
  }
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
  for (register int i = 0; i < nBubbles; i++)
  {
    Bubbles[i]->Gas = Gas;
    Bubbles[i]->Liquid = Liquid;
    Bubbles[i]->Interface = Interface;
    Bubbles[i]->NumericsODE = NumericsODE;
  }

  /* Allocate and set the excitation parameters */
  struct APECSS_Excitation *Excitation = (struct APECSS_Excitation *) malloc(sizeof(struct APECSS_Excitation));
  Excitation->type = APECSS_EXCITATION_SIN;
  Excitation->f = (APECSS_FLOAT) fa;
  Excitation->dp = (APECSS_FLOAT) pa;
  for (register int i = 0; i < nBubbles; i++) Bubbles[i]->Excitation = Excitation;

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Create individual folders for the results of each bubble
  for (register int i = 0; i < nBubbles; i++)
  {
    if (Bubbles[i]->Results != NULL)
    {
      sprintf(Bubbles[i]->Results->dir, "./Bubble_%i/", i);
      struct stat st = {0};
      if (stat(Bubbles[i]->Results->dir, &st) == -1) mkdir(Bubbles[i]->Results->dir, 0700);
    }
  }
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Process all options */
  apecss_gas_processoptions(Gas);
  apecss_liquid_processoptions(Liquid);
  apecss_interface_processoptions(Interface);
  apecss_odesolver_processoptions(NumericsODE);
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_processoptions(Bubbles[i]);

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Use the revised pressure at infinity, including neighbor contributions
  for (register int i = 0; i < nBubbles; i++) Bubbles[i]->get_pressure_infinity = interaction_bubble_pressure_infinity;

  // Initialisize interaction structure
  for (register int i = 0; i < nBubbles; i++) Bubbles[i]->Interaction = (struct APECSS_Interaction *) malloc(sizeof(struct APECSS_Interaction));

  // Update interaction structure
  for (register int i = 0; i < nBubbles; i++) Bubbles[i]->Interaction->nBubbles = nBubbles;

  // Define the size of each bubble
  for (register int i = 0; i < nBubbles; i++)
  {
    if (0 == i)
      Bubbles[i]->R0 = 5.0e-6;
    else if (1 == i)
      Bubbles[i]->R0 = 10.0e-6;

    Bubbles[i]->r_hc = Bubbles[i]->R0 / 8.54;
  }

  // Define center location for each bubble
  for (register int i = 0; i < nBubbles; i++)
  {
    if (0 == i)
    {
      Bubbles[i]->Interaction->location[0] = 0.0;
      Bubbles[i]->Interaction->location[1] = 0.0;
      Bubbles[i]->Interaction->location[2] = 0.0;
    }
    else if (1 == i)
    {
      Bubbles[i]->Interaction->location[0] = bubble_bubble_dist;
      Bubbles[i]->Interaction->location[1] = 0.0;
      Bubbles[i]->Interaction->location[2] = 0.0;
    }
  }
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  clock_t starttimebubble = clock();
  APECSS_FLOAT tSim = 0.0;  // Simulation time of the coupled system

  for (register int i = 0; i < nBubbles; i++)
  {
    Bubbles[i]->tStart = tSim;
    Bubbles[i]->tEnd = (APECSS_FLOAT) tEnd;
    Bubbles[i]->dt = APECSS_MIN(1.0e-11, dt_interbubble);  // Initial time-step
  }

  /* Initialize the bubble based on the selected options */
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_initialize(Bubbles[i]);

  /* Initialize */
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_solver_initialize(Bubbles[i]);

  /* Solve the bubble dynamics */
  while (tSim < (APECSS_FLOAT) tEnd)  // Interaction loop, corresponding to the time-intervals at which interactions are considered
  {
    APECSS_FLOAT dtSim = APECSS_MIN(dt_interbubble, (APECSS_FLOAT) tEnd - tSim);
    tSim += dtSim;

    for (register int i = 0; i < nBubbles; i++) apecss_bubble_solver_run(tSim, Bubbles[i]);

    for (register int i = 0; i < nBubbles; i++) Bubbles[i]->Interaction->dp_neighbor = 0.0;

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Update the contribution of the neighbor bubble
    apecss_instantaneous_interactions(Bubbles);

    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }

  /* Finalize the simulation*/
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_solver_finalize(Bubbles[i]);

  char str[APECSS_STRINGLENGTH_SPRINTF];
  for (register int i = 0; i < nBubbles; i++)
  {
    sprintf(str, "Bubble %i: Solver concluded %i time-steps and %i sub-iterations in %.3f s.", i, Bubbles[i]->dtNumber, Bubbles[i]->nSubIter,
            (double) (clock() - starttimebubble) / CLOCKS_PER_SEC);
    apecss_writeonscreen(str);
  }

  for (register int i = 0; i < nBubbles; i++) apecss_results_rayleighplesset_write(Bubbles[i], APECSS_RESULTS_WRITE);
  for (register int i = 0; i < nBubbles; i++) apecss_results_emissionsspace_write(Bubbles[i], APECSS_RESULTS_WRITE);

  /* Make sure all allocated memory is freed */
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_freestruct(Bubbles[i]);

  for (register int i = 0; i < nBubbles; i++) free(Bubbles[i]);
  free(Gas);
  free(Liquid);
  free(Interface);
  free(NumericsODE);
  free(Excitation);

  return (0);
}

APECSS_FLOAT interaction_bubble_pressure_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  struct Interaction *temp_struct = Bubble->user_data;
  return (Bubble->p0 - Bubble->Excitation->dp * APECSS_SIN(2.0 * APECSS_PI * Bubble->Excitation->f * t) + temp_struct->dp_neighbor);
}