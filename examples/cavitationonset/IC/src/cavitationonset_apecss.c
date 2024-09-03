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
// APECSS standalone example of cavitation onset with acoustically-in-
// teracting microbubbles, based on Ida (2009), Physics of Fluids 21
// (11), 113302, DOI : 10.1063/1.3265547
// -------------------------------------------------------------------

#include <time.h>
#include "apecss.h"

// Declaration of additional case-dependent functions
APECSS_FLOAT interaction_bubble_pressure_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT interaction_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

// Declaration of the structure holding the interaction variables of each bubble
struct Interaction
{
  APECSS_FLOAT dp_neighbor;
};

int main(int argc, char **args)
{
  char OptionsDir[APECSS_STRINGLENGTH];

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Interbubble time-step, defining the frequency with which the neighbor influence is updated
  APECSS_FLOAT dt_interbubble = 1.0e-8;

  // Initialize the simulation parameters given by the execution command
  int nBubbles = 2;
  double tEnd = 0.0;
  double fa = 0.0;
  double pa = 0.0;
  int cluster_distrib = 0;
  double cluster_size = 0;
  int inttype = 0;
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
    else if (strcmp("-nbb", args[j]) == 0)
    {
      sscanf(args[j + 1], "%d", &nBubbles);
      j += 2;
    }
    else if (strcmp("-cldistrib", args[j]) == 0)
    {
      sscanf(args[j + 1], "%d", &cluster_distrib);
      j += 2;
    }
    else if (strcmp("-clsize", args[j]) == 0)
    {
      sscanf(args[j + 1], "%le", &cluster_size);
      j += 2;
    }
    else if (strcmp("-inttype", args[j]) == 0)
    {
      sscanf(args[j + 1], "%d", &inttype);
      j += 2;
    }
    else if (strcmp("-dt_inter", args[j]) == 0)
    {
      sscanf(args[j + 1], "%le", &dt_interbubble);
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
  dt_interbubble = (APECSS_FLOAT) dt_interbubble;

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
  if (cluster_distrib == 0)
  {
    // Two bubbles of different initial size interacting
    Bubbles[0]->R0 = 2.0e-06;
    Bubbles[1]->R0 = 20.0e-06;
  }
  else
  {
    // Monodispersed multibubble distributions
    for (register int i = 0; i < nBubbles; i++)
    {
      Bubbles[i]->R0 = 20.0e-06;
    }
  }

  // Define center location for each bubble
  if (cluster_distrib == 0)
  {
    // Two bubbles of different initial size interacting
    Bubbles[0]->Interaction->location[0] = 0.0;
    Bubbles[0]->Interaction->location[1] = 0.0;
    Bubbles[0]->Interaction->location[2] = 0.0;

    Bubbles[1]->Interaction->location[0] = (APECSS_FLOAT) cluster_size * (Bubbles[0]->R0 + Bubbles[1]->R0);
    Bubbles[1]->Interaction->location[1] = 0.0;
    Bubbles[1]->Interaction->location[2] = 0.0;
  }
  else
  {
    // Monodispersed multibubble distributions
    APECSS_FLOAT D = 400.0e-06;
    if (nBubbles == 1)
    {
      // Single bubble
      Bubbles[0]->Interaction->location[0] = 0.0;
      Bubbles[0]->Interaction->location[1] = 0.0;
      Bubbles[0]->Interaction->location[2] = 0.0;
    }
    else if (nBubbles == 2)
    {
      // Two bubbles
      Bubbles[0]->Interaction->location[0] = 0.0;
      Bubbles[0]->Interaction->location[1] = 0.0;
      Bubbles[0]->Interaction->location[2] = 0.0;

      Bubbles[1]->Interaction->location[0] = D;
      Bubbles[1]->Interaction->location[1] = 0.0;
      Bubbles[1]->Interaction->location[2] = 0.0;
    }
    else if (nBubbles == 3)
    {
      // Regular triangle
      Bubbles[0]->Interaction->location[0] = 0.0;
      Bubbles[0]->Interaction->location[1] = 0.0;
      Bubbles[0]->Interaction->location[2] = 0.0;

      Bubbles[1]->Interaction->location[0] = D;
      Bubbles[1]->Interaction->location[1] = 0.0;
      Bubbles[1]->Interaction->location[2] = 0.0;

      Bubbles[2]->Interaction->location[0] = 0.5 * D;
      Bubbles[2]->Interaction->location[1] = D * APECSS_SIN(APECSS_ONETHIRD * APECSS_PI);
      Bubbles[2]->Interaction->location[2] = 0.0;
    }
    else if (nBubbles == 4)
    {
      // Regular tetragon
      Bubbles[0]->Interaction->location[0] = 0.0;
      Bubbles[0]->Interaction->location[1] = 0.0;
      Bubbles[0]->Interaction->location[2] = 0.0;

      Bubbles[1]->Interaction->location[0] = D;
      Bubbles[1]->Interaction->location[1] = 0.0;
      Bubbles[1]->Interaction->location[2] = 0.0;

      Bubbles[2]->Interaction->location[0] = 0.0;
      Bubbles[2]->Interaction->location[1] = D;
      Bubbles[2]->Interaction->location[2] = 0.0;

      Bubbles[3]->Interaction->location[0] = D;
      Bubbles[3]->Interaction->location[1] = D;
      Bubbles[3]->Interaction->location[2] = 0.0;
    }
    else if (nBubbles == 8)
    {
      // Regular hexaedron
      Bubbles[0]->Interaction->location[0] = 0.0;
      Bubbles[0]->Interaction->location[1] = 0.0;
      Bubbles[0]->Interaction->location[2] = 0.0;

      Bubbles[1]->Interaction->location[0] = D;
      Bubbles[1]->Interaction->location[1] = 0.0;
      Bubbles[1]->Interaction->location[2] = 0.0;

      Bubbles[2]->Interaction->location[0] = 0.0;
      Bubbles[2]->Interaction->location[1] = D;
      Bubbles[2]->Interaction->location[2] = 0.0;

      Bubbles[3]->Interaction->location[0] = D;
      Bubbles[3]->Interaction->location[1] = D;
      Bubbles[3]->Interaction->location[2] = 0.0;

      Bubbles[4]->Interaction->location[0] = 0.0;
      Bubbles[4]->Interaction->location[1] = 0.0;
      Bubbles[4]->Interaction->location[2] = D;

      Bubbles[5]->Interaction->location[0] = D;
      Bubbles[5]->Interaction->location[1] = 0.0;
      Bubbles[5]->Interaction->location[2] = D;

      Bubbles[6]->Interaction->location[0] = 0.0;
      Bubbles[6]->Interaction->location[1] = D;
      Bubbles[6]->Interaction->location[2] = D;

      Bubbles[7]->Interaction->location[0] = D;
      Bubbles[7]->Interaction->location[1] = D;
      Bubbles[7]->Interaction->location[2] = D;
    }
    else
    {
      char str[APECSS_STRINGLENGTH_SPRINTF];
      sprintf(str, "Wrong number of bubbles as input for this specific cluster distribution (1, 2, 3, 4 or 8)");
      apecss_erroronscreen(1, str);
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
  // File to retrieve all valuable information for cavitation onset test case
  FILE *file_ida2009;
  file_ida2009 = fopen("Ida2009_results.txt", "w");
  fprintf(file_ida2009, "%d Bubbles p0(pa) %e png(Pa) %e D_multiplier(-) %e cl_distrib %d Interaction-type %d\n", nBubbles, Liquid->pref, pa, cluster_size,
          cluster_distrib, inttype);
  fprintf(file_ida2009, "Initial_radii(m)");
  for (register int i = 0; i < nBubbles; i++) fprintf(file_ida2009, " %e", Bubbles[i]->R0);
  fprintf(file_ida2009, "\n");
  fprintf(file_ida2009, "#Time(s) R(m) Pt(Pa)\n");
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Solve the bubble dynamics */
  while (tSim < (APECSS_FLOAT) tEnd)  // Interaction loop, corresponding to the time-intervals at which interactions are considered
  {
    APECSS_FLOAT dtSim = APECSS_MIN(dt_interbubble, (APECSS_FLOAT) tEnd - tSim);
    tSim += dtSim;

    // printf("%e", tSim);
    // printf(" %e", Bubbles[0]->ode[0](Bubbles[0]->ODEsSol, Bubbles[0]->t, Bubbles[0]));
    for (register int i = 0; i < nBubbles; i++) apecss_bubble_solver_run(tSim, Bubbles[i]);

    // for (register int i = 0; i < nBubbles; i++) Bubbles[i]->Interaction->dp_neighbor = 0.0;

    // printf(" %e", Bubbles[0]->R);
    // printf(" %e", Bubbles[0]->ode[0](Bubbles[0]->ODEsSol, Bubbles[0]->t, Bubbles[0]));
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Update the contribution of the neighbor bubbles
    if (inttype == 1)
    {
      apecss_interactions_instantaneous(Bubbles);
    }
    else if (inttype == 2)
    {
      apecss_interactions_quasi_acoustic(Bubbles);
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // for (register int i = 0; i < nBubbles; i++)
    // {
    //   printf(" %e", Bubbles[i]->get_pressure_infinity(tSim, Bubbles[i]));
    // }
    // printf(" %e", Bubbles[0]->get_pressure_infinity(Bubbles[0]->t, Bubbles[0]));
    // printf(" %e", Bubbles[0]->ode[0](Bubbles[0]->ODEsSol, Bubbles[0]->t, Bubbles[0]));
    // printf("\n");

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Retrieve data
    fprintf(file_ida2009, "%e", tSim);
    for (register int i = 0; i < nBubbles; i++)
    {
      fprintf(file_ida2009, " %e", Bubbles[i]->R);
    }
    for (register int i = 0; i < nBubbles; i++)
    {
      fprintf(file_ida2009, " %e", Bubbles[i]->get_pressure_infinity(tSim, Bubbles[i]));
    }
    for (register int i = 0; i < nBubbles; i++)
    {
      fprintf(file_ida2009, " %e", Bubbles[i]->ode[0](Bubbles[i]->ODEsSol, Bubbles[i]->t, Bubbles[i]));
    }
    for (register int i = 0; i < nBubbles; i++)
    {
      fprintf(file_ida2009, " %e", Bubbles[i]->Interaction->dp_neighbor);
    }
    for (register int i = 0; i < nBubbles; i++)
    {
      fprintf(file_ida2009, " %e", Bubbles[i]->U);
    }
    fprintf(file_ida2009, "\n");
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // int index = 0;
    // // APECSS_FLOAT derivative = (Bubbles[index]->Interaction->last_p_1 - Bubbles[index]->Interaction->last_p_2) /
    // // (Bubbles[index]->Interaction->last_t_1 - Bubbles[index]->Interaction->last_t_2);
    // printf("%e Bubble %d R %e U %e A %e t_1 %e t_2 %e inv_t %e p_1 %e p_2 %e diff_p %e derivative %e pinfty %e\n", tSim, index, Bubbles[index]->R,
    //        Bubbles[index]->U, Bubbles[index]->ode[0](Bubbles[index]->ODEsSol, Bubbles[index]->t, Bubbles[index]), Bubbles[index]->Interaction->last_t_1,
    //        Bubbles[index]->Interaction->last_t_2, 1 / (Bubbles[index]->Interaction->last_t_1 - Bubbles[index]->Interaction->last_t_2),
    //        Bubbles[index]->Interaction->last_p_1, Bubbles[index]->Interaction->last_p_2,
    //        (Bubbles[index]->Interaction->last_p_1 - Bubbles[index]->Interaction->last_p_2),
    //        Bubbles[index]->get_pressurederivative_infinity(tSim, Bubbles[index]), Bubbles[index]->get_pressure_infinity(tSim, Bubbles[index]));

    for (register int i = 0; i < nBubbles; i++)
    {
      Bubbles[i]->Interaction->last_t_2 = Bubbles[i]->Interaction->last_t_1;
      Bubbles[i]->Interaction->last_p_2 = Bubbles[i]->Interaction->last_p_1;

      Bubbles[i]->Interaction->last_t_1 = tSim;
      Bubbles[i]->Interaction->last_p_1 = Bubbles[i]->Interaction->dp_neighbor;
    }
  }

  fclose(file_ida2009);

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
  APECSS_FLOAT T = 10.0e-06;
  if (t < T)
  {
    return (Bubble->p0 + Bubble->Interaction->dp_neighbor);
  }
  else if (t > 2 * T)
  {
    return (Bubble->Excitation->dp + Bubble->Interaction->dp_neighbor);
  }
  else
  {
    APECSS_FLOAT W = 0.5 * (1 - APECSS_COS((t + T) * APECSS_PI / T));
    return (Bubble->p0 + W * (Bubble->Excitation->dp - Bubble->p0) + Bubble->Interaction->dp_neighbor);
  }
}

APECSS_FLOAT interaction_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  // Approximate numerical computation of p_infinity derivative
  APECSS_FLOAT T = 10.0e-06;
  APECSS_FLOAT derivative = 0.0;
  if ((t >= T) && (t <= 2 * T))
  {
    APECSS_FLOAT inv_T = 1 / T;
    derivative = 0.5 * APECSS_PI * inv_T * APECSS_SIN(APECSS_PI * (t + T) * inv_T) * (Bubble->Excitation->dp - Bubble->p0);
  }
  return (derivative);

  APECSS_FLOAT delta_t = Bubble->Interaction->last_t_1 - Bubble->Interaction->last_t_2;
  if (delta_t > Bubble->dt)
  {
    APECSS_FLOAT inv_delta_t = 1 / delta_t;
    return (derivative + ((Bubble->Interaction->last_p_1 - Bubble->Interaction->last_p_2) * inv_delta_t));
  }
  else
  {
    return (derivative);
  }
}