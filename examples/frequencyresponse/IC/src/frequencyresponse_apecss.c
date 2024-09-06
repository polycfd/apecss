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
APECSS_FLOAT interaction_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_rp_kellermiksisvelocity_ode_simplified(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
int apecss_new_bubble_solver_run(APECSS_FLOAT tend, APECSS_FLOAT tEnd, struct APECSS_Bubble *Bubble);

APECSS_FLOAT APECSS_TAN(APECSS_FLOAT x) { return APECSS_SIN(x) / APECSS_COS(x); }

APECSS_FLOAT maxR;

int main(int argc, char **args)
{
  char OptionsDir[APECSS_STRINGLENGTH];

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Set the case-dependent simulation parameters
  int nBubbles = 3;  // Number of bubbles
  APECSS_FLOAT bubble_bubble_dist = 5.0e-6;  // Bubble-bubble distance

  // Interbubble time-step, defining the frequency with which the neighbor influence is updated
  APECSS_FLOAT dt_interbubble = 1.0e-8;

  // Initialize the simulation parameters given by the execution command
  double tEnd = 0.0;
  double fa = 0.0;
  double pa = 0.0;
  double dist = 5.0e-6;
  double dt_inter = 1.0e-8;
  int ode = 0;
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
    else if (strcmp("-nb", args[j]) == 0)
    {
      sscanf(args[j + 1], "%d", &nBubbles);
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
    else if (strcmp("-dist", args[j]) == 0)
    {
      sscanf(args[j + 1], "%le", &dist);
      j += 2;
    }
    else if (strcmp("-dt_inter", args[j]) == 0)
    {
      sscanf(args[j + 1], "%le", &dt_inter);
      j += 2;
    }
    else if (strcmp("-ode", args[j]) == 0)
    {
      sscanf(args[j + 1], "%d", &ode);
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

  /* Properly update parameters with commandline options */
  bubble_bubble_dist = (APECSS_FLOAT) dist;
  dt_interbubble = (APECSS_FLOAT) dt_inter;
  tEnd = 30 / fa;

  /* Check if the number of bubbles is right */
  if ((nBubbles != 3) && (nBubbles != 4))
  {
    nBubbles = 3;
  }

  /* Allocate and initialize Bubble structure */
  struct APECSS_Bubble *Bubbles[nBubbles];
  for (register int i = 0; i < nBubbles; i++) Bubbles[i] = (struct APECSS_Bubble *) malloc(nBubbles * sizeof(struct APECSS_Bubble));
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_initializestruct(Bubbles[i]);

  /* Set default options and read the options for the bubble */
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_setdefaultoptions(Bubbles[i]);
  for (register int i = 0; i < nBubbles; i++) apecss_bubble_readoptions(Bubbles[i], OptionsDir);

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
  for (register int i = 0; i < nBubbles; i++) Bubbles[i]->get_pressurederivative_infinity = interaction_bubble_pressurederivative_infinity;

  // Update the Keller-Miksis ODE
  if (ode == 1)
  {
    for (register int i = 0; i < nBubbles; i++) Bubbles[i]->ode[0] = apecss_rp_kellermiksisvelocity_ode_simplified;
  }

  // Initialize interaction structure
  for (register int i = 0; i < nBubbles; i++) Bubbles[i]->Interaction = (struct APECSS_Interaction *) malloc(sizeof(struct APECSS_Interaction));

  // Update interaction structure
  for (register int i = 0; i < nBubbles; i++)
  {
    Bubbles[i]->Interaction->nBubbles = nBubbles;
    Bubbles[i]->Interaction->dp_neighbor = 0.0;
    Bubbles[i]->Interaction->last_t_1 = 0.0;
    Bubbles[i]->Interaction->last_t_2 = 0.0;
    Bubbles[i]->Interaction->last_p_1 = 0.0;
    Bubbles[i]->Interaction->last_p_2 = 0.0;
  }

  // Define the size of each bubble
  for (register int i = 0; i < nBubbles; i++)
  {
    if (0 == i)
      Bubbles[i]->R0 = 1.0e-6;
    else if (1 == i)
      Bubbles[i]->R0 = 0.8e-6;
    else if (2 == i)
      Bubbles[i]->R0 = 0.5e-6;
    else if (3 == i)
      Bubbles[i]->R0 = 1.5e-6;
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
    else if (2 == i)
    {
      Bubbles[i]->Interaction->location[0] = 0.5 * bubble_bubble_dist;
      Bubbles[i]->Interaction->location[1] = -APECSS_SIN(APECSS_PI / 3) * bubble_bubble_dist;
      Bubbles[i]->Interaction->location[2] = 0.0;
    }
    else if (3 == i)
    {
      Bubbles[i]->Interaction->location[0] = 0.5 * bubble_bubble_dist;
      Bubbles[i]->Interaction->location[1] = -APECSS_TAN(APECSS_PI / 6) * 0.5 * bubble_bubble_dist;
      Bubbles[i]->Interaction->location[2] = APECSS_SQRT(APECSS_POW2(bubble_bubble_dist) - APECSS_POW2(0.5 * bubble_bubble_dist) -
                                                         APECSS_POW2(APECSS_TAN(APECSS_PI / 6) * 0.5 * bubble_bubble_dist));
    }
  }

  // Create array to gather maximum radius value for each bubble
  APECSS_FLOAT max_bubble[nBubbles];
  for (register int i = 0; i < nBubbles; i++)
  {
    max_bubble[i] = Bubbles[i]->R0;
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

    for (register int i = 0; i < nBubbles; i++)
    {
      maxR = Bubbles[i]->R0;
      apecss_new_bubble_solver_run(tSim, (APECSS_FLOAT) tEnd, Bubbles[i]);
      if (maxR > max_bubble[i])
      {
        max_bubble[i] = maxR;
      }
    }

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Update the contribution of the neighbor bubbles
    apecss_interactions_instantaneous(Bubbles);
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }

  /* Retrieve results */
  FILE *file_max_radius;
  file_max_radius = fopen("max_radius.txt", "a");
  fprintf(file_max_radius, "nb %d f(Hz) %e p(Pa) %e dist %e ODE %d R0(m);Rmax(m)", nBubbles, fa, pa, dist, ode);
  for (register int i = 0; i < nBubbles; i++)
  {
    fprintf(file_max_radius, " %e;%e", Bubbles[i]->R0, max_bubble[i]);
  }
  fprintf(file_max_radius, "\n");
  fclose(file_max_radius);

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
  APECSS_FLOAT x = Bubble->Interaction->location[0];
  return (Bubble->p0 - Bubble->Excitation->dp * APECSS_SIN(2.0 * APECSS_PI * Bubble->Excitation->f * (t - x / Bubble->Liquid->cref)) +
          Bubble->Interaction->dp_neighbor);
}

APECSS_FLOAT interaction_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  // Approximate numerical computation of p_infinity derivative
  APECSS_FLOAT x = Bubble->Interaction->location[0];
  APECSS_FLOAT derivative =
      -Bubble->Excitation->dp * 2.0 * APECSS_PI * Bubble->Excitation->f * APECSS_COS(2.0 * APECSS_PI * Bubble->Excitation->f * (t - x / Bubble->Liquid->cref));
  return (derivative);
}

APECSS_FLOAT apecss_rp_kellermiksisvelocity_ode_simplified(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  /** Simplified model for comparison with Haghi **/
  APECSS_FLOAT inv_c = 1.0 / Bubble->Liquid->cref;
  APECSS_FLOAT inv_rho = 1.0 / Bubble->Liquid->rhoref;

  APECSS_FLOAT rhs = ((Bubble->Liquid->get_pressure_bubblewall(Sol, t, Bubble) - Bubble->get_pressure_infinity(t, Bubble)) * inv_rho -
                      1.5 * (1.0 - (Sol[0] * APECSS_ONETHIRD * inv_c)) * APECSS_POW2(Sol[0])) /
                     Sol[1];

  return (rhs / (1.0 - Sol[0] * inv_c));
}

int apecss_new_bubble_solver_run(APECSS_FLOAT tend, APECSS_FLOAT tEnd, struct APECSS_Bubble *Bubble)
{
  while (Bubble->t < tend - 0.01 * Bubble->NumericsODE->dtMin)
  {
    // Store the previous solution for sub-iterations
    for (register int i = 0; i < Bubble->nODEs; i++) Bubble->ODEsSolOld[i] = Bubble->ODEsSol[i];

    // Check what comes next: the end of the solver run or, if applicable, a time instance to output the emissions
    APECSS_FLOAT next_event_time = Bubble->results_emissionstime_check(tend, Bubble);

    // Set the time-step for the ODEs
    apecss_odesolver_settimestep(Bubble->NumericsODE, Bubble->err, next_event_time - Bubble->t, &(*Bubble).dt);

    // Solve the ODEs
    Bubble->err = apecss_odesolver(Bubble);

    // Perform sub-iterations on the control of the time-step when err > tol
    register int subiter = 0;
    while ((Bubble->err > Bubble->NumericsODE->tol) && (subiter < Bubble->NumericsODE->maxSubIter))
    {
      ++subiter;
      // Rewind the solution
      for (register int i = 0; i < Bubble->nODEs; i++) Bubble->ODEsSol[i] = Bubble->ODEsSolOld[i];
      // Set the time-step for the ODEs
      apecss_odesolver_settimestep(Bubble->NumericsODE, Bubble->err, next_event_time - Bubble->t, &(*Bubble).dt);
      // Solve the ODEs again
      Bubble->err = apecss_odesolver(Bubble);
    }
    Bubble->nSubIter += subiter;

    // Set new values
    ++(Bubble->dtNumber);
    Bubble->t += Bubble->dt;
    Bubble->U = Bubble->ODEsSol[0];
    Bubble->R = Bubble->ODEsSol[1];

    // Detect max value of radius
    if (Bubble->t > tEnd - 10 / Bubble->Excitation->f) maxR = APECSS_MAX(maxR, Bubble->R);

    // Update allocation of the results (if necessary)
    Bubble->results_emissionsnodeminmax_identify(Bubble);
    Bubble->results_emissionsnode_alloc(Bubble);

    // Acoustic emissions (if applicable)
    Bubble->emissions_update(Bubble);

    // Store results (if applicable)
    Bubble->results_rayleighplesset_store(Bubble);
    Bubble->results_emissionsspace_store(Bubble);

    // Write one-off results (if applicable)
    Bubble->results_emissionstime_write(Bubble);

    // Update the last-step solution of the RK54 scheme
    for (register int i = 0; i < Bubble->nODEs; i++) Bubble->kLast[i] = Bubble->k7[i];

    // Update progress screen in the terminal (if applicable)
    Bubble->progress_update(&(*Bubble).progress, Bubble->t - Bubble->tStart, Bubble->tEnd - Bubble->tStart);
  }

  return (0);
}