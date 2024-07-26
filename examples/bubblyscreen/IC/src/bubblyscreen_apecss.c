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
// APECSS standalone example of acoustically-interacting microbubbles
// in a finite bubbly screen based on Fan, Y., Li, H., & Fuster,
// D. (2021). Time-delayed interactions on acoustically driven bubbly
// screens. The Journal of the Acoustical Society of America, 150(6),
// 4219‑4231. https://doi.org/10.1121/10.0008905
// -------------------------------------------------------------------

#include <time.h>
#include "apecss.h"

// Declaration of additional case-dependent functions
APECSS_FLOAT interaction_bubble_pressure_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT interaction_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
int apecss_bubble_new_solver_run(APECSS_FLOAT tend, APECSS_FLOAT tEnd, struct APECSS_Bubble *Bubble);

APECSS_FLOAT maxR;
APECSS_FLOAT minR;

int main(int argc, char **args)
{
  char OptionsDir[APECSS_STRINGLENGTH];

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Set the case-dependent simulation parameters
  const int nBubbles_x = 51;
  const int nBubbles = nBubbles_x * nBubbles_x;  // Number of bubbles
  APECSS_FLOAT bubble_bubble_dist = 400.0e-6;  // Bubble-bubble distance

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
    Bubbles[i]->R0 = 1.0e-6;
    // Bubbles[i]->r_hc = Bubbles[i]->R0 / 8.54;
  }

  // Define center location for each bubble
  int space = (int) (0.5 * (nBubbles_x - 1));
  int row = 0;
  int col = 0;
  for (register int i = 0; i < nBubbles; i++)
  {
    Bubbles[i]->Interaction->location[0] = (col - space) * bubble_bubble_dist;
    Bubbles[i]->Interaction->location[1] = (row - space) * bubble_bubble_dist;
    Bubbles[i]->Interaction->location[2] = 0.0;
    if ((col < 2 * space) && (col >= 0))
    {
      col += 1;
    }
    else
    {
      row += 1;
      col = 0;
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

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Lists to gather maximum and minimum radius achieved by each bubble in steady state
  APECSS_FLOAT max_radii[nBubbles];
  APECSS_FLOAT min_radii[nBubbles];

  for (register int i = 0; i < nBubbles; i++)
  {
    max_radii[i] = Bubbles[i]->R;
    min_radii[i] = Bubbles[i]->R;
  }

  // File to retrieve the maximum and minimum radius evolution achieved in steady state
  FILE *file_extremum;
  file_extremum = fopen("bubblyscreen_extremum.txt", "w");

  APECSS_FLOAT p0 = Bubbles[0]->p0;
  APECSS_FLOAT poly = Gas->Gamma;
  APECSS_FLOAT rho = Liquid->rhoref;
  APECSS_FLOAT w0 = APECSS_SQRT((3 * poly * p0) / (APECSS_POW2(Bubbles[0]->R0) * rho));

  fprintf(file_extremum, "w0(rad/s) %e f(Hz) %e p(Pa) %e D/R0 %e\n", w0, fa, pa, bubble_bubble_dist / Bubbles[0]->R0);
  fprintf(file_extremum, "label x(m) y(m) z(m) R0(m) Rmin(m) Rmax(m)\n");

  // File to retrieve the radius evolution during computation
  FILE *file_radii;
  file_radii = fopen("bubblyscreen_radii.txt", "w");
  fprintf(file_radii, "w0(rad/s) %e f(Hz) %e p(Pa) %e D/R0 %e\n", w0, fa, pa, bubble_bubble_dist / Bubbles[0]->R0);
  fprintf(file_radii, "t(s) R(m)\n");
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Solve the bubble dynamics */
  while (tSim < (APECSS_FLOAT) tEnd)  // Interaction loop, corresponding to the time-intervals at which interactions are considered
  {
    APECSS_FLOAT dtSim = APECSS_MIN(dt_interbubble, (APECSS_FLOAT) tEnd - tSim);
    tSim += dtSim;

    for (register int i = 0; i < nBubbles; i++)
    {
      maxR = 0.0;
      minR = Bubbles[i]->R0;
      apecss_bubble_new_solver_run(tSim, tEnd, Bubbles[i]);

      if (maxR > max_radii[i])
      {
        max_radii[i] = maxR;
      }
      if (minR < min_radii[i])
      {
        min_radii[i] = minR;
      }
    }

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Update the contribution of the neighbor bubble
    apecss_interactions_instantaneous(Bubbles);
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    fprintf(file_radii, "%e", tSim);
    for (register int i = 0; i < nBubbles; i++)
    {
      fprintf(file_radii, " %e", Bubbles[i]->R);
    }
    fprintf(file_radii, "\n");

    // for (register int i = 0; i < nBubbles; i++)
    // {
    //   Bubbles[i]->Interaction->last_t_2 = Bubbles[i]->Interaction->last_t_1;
    //   Bubbles[i]->Interaction->last_p_2 = Bubbles[i]->Interaction->last_p_1;

    //   Bubbles[i]->Interaction->last_t_1 = tSim;
    //   Bubbles[i]->Interaction->last_p_1 = Bubbles[i]->Interaction->dp_neighbor;
    // }
  }

  /* Complete results file */
  for (register int i = 0; i < nBubbles; i++)
  {
    fprintf(file_extremum, "%d %e %e %e %e %e %e\n", i, Bubbles[i]->Interaction->location[0], Bubbles[i]->Interaction->location[1],
            Bubbles[i]->Interaction->location[2], Bubbles[i]->R0, min_radii[i], max_radii[i]);
  }
  fclose(file_extremum);

  fclose(file_radii);

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
  return (Bubble->p0 - Bubble->Excitation->dp * APECSS_SIN(2.0 * APECSS_PI * Bubble->Excitation->f * t) + Bubble->Interaction->dp_neighbor);
}

APECSS_FLOAT interaction_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  // Approximate numerical computation of p_infinity derivative
  // APECSS_FLOAT delta_t = Bubble->Interaction->last_t_1 - Bubble->Interaction->last_t_2;
  APECSS_FLOAT derivative = -Bubble->Excitation->dp * 2.0 * APECSS_PI * Bubble->Excitation->f * APECSS_COS(2.0 * APECSS_PI * Bubble->Excitation->f * t);
  return (derivative);
  // if (delta_t > Bubble->dt)
  // {
  //   return (derivative + ((Bubble->Interaction->last_p_1 - Bubble->Interaction->last_p_2) / (Bubble->Interaction->last_t_1 - Bubble->Interaction->last_t_2)));
  // }
  // else
  // {
  //   return (derivative);
  // }
}

int apecss_bubble_new_solver_run(APECSS_FLOAT tend, APECSS_FLOAT tEnd, struct APECSS_Bubble *Bubble)
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

    // Retrieve data to compute radius oscillation amplitude during
    if (Bubble->t > tEnd - 10.0 / Bubble->Excitation->f)
    {
      maxR = APECSS_MAX(maxR, Bubble->R);
    }
    if (Bubble->t > tEnd - 10.0 / Bubble->Excitation->f)
    {
      minR = APECSS_MIN(minR, Bubble->R);
    }

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