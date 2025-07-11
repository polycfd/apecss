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

// -------------------------------------------------------------------
// APECSS standalone example of acoustically-interacting microbubbles
// in a spherical cluster with the goal to study cavitation onset
// -------------------------------------------------------------------

#include <time.h>
#include <mpi.h>
#include "apecss.h"

APECSS_FLOAT rand_range(double min, double max)
{
  double random = (drand48());
  double range = (max - min) * random;
  double number = min + range;

  return (APECSS_FLOAT) number;
}

struct APECSS_Parallel_Cluster
{
  int rank, size;
  int nBubbles_local, nBubbles_global;

  int *bubblerank;  // max id of bubbles for each rank
  APECSS_FLOAT *bubbleglobal_R, *bubbleglobal_x, *bubbleglobal_y, *bubbleglobal_z;  // Radius and x,y,z location of all bubbles
  APECSS_FLOAT *sumGU_rank;  // Contributions from local bubbles to all bubbles
};

// Declaration of additional case-dependent functions
APECSS_FLOAT parallel_interactions_bubble_pressure_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT parallel_interactions_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

int parallel_interactions_quasi_acoustic(struct APECSS_Bubble *Bubble[], struct APECSS_Parallel_Cluster *RankInfo);
int parallel_interactions_proper_cutoffdistance(struct APECSS_Bubble *Bubbles[], struct APECSS_Parallel_Cluster *RankInfo);

int main(int argc, char **args)
{
  /* Initialize MPI */
  int mpi_rank, mpi_size;
  MPI_Init(&argc, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  char OptionsDir[APECSS_STRINGLENGTH];

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Set the case-dependent simulation parameters
  const int nBubbles = 250;  // Number of bubbles
  APECSS_FLOAT bubble_bubble_dist = 20.0e-6;  // Bubble-bubble minimal distance
  APECSS_FLOAT cluster_radius = 232e-6;  // Spherical cluster radius

  // Interbubble time-step, defining the frequency with which the neighbor influence is updated
  APECSS_FLOAT dt_interbubble = 1.0e-9;

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

  /* Allocate structure for parallel data */
  struct APECSS_Parallel_Cluster *RankInfo = (struct APECSS_Parallel_Cluster *) malloc(sizeof(struct APECSS_Parallel_Cluster));
  RankInfo->rank = mpi_rank;
  RankInfo->size = mpi_size;

  /* Determine the number of bubbles per rank */
  int max_per_rank = ceil((double) nBubbles / (double) mpi_size);
  RankInfo->nBubbles_local = APECSS_MAX(0, APECSS_MIN(max_per_rank, nBubbles - mpi_rank * max_per_rank));

  /* Share the parallel distribution of bubbles with all ranks */
  RankInfo->bubblerank = malloc((RankInfo->size + 1) * sizeof(int));

  RankInfo->nBubbles_global = 0;
  RankInfo->bubblerank[0] = 0;
  for (int r = 0; r < RankInfo->size; r++)
  {
    int temp = RankInfo->nBubbles_local;
    MPI_Bcast(&temp, 1, MPI_INT, r, MPI_COMM_WORLD);
    RankInfo->bubblerank[r + 1] = RankInfo->bubblerank[r] + temp;
    RankInfo->nBubbles_global += temp;
  }

  /* Allocate and initialize Bubble structure */
  struct APECSS_Bubble *Bubbles[RankInfo->nBubbles_local];
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i] = (struct APECSS_Bubble *) malloc(sizeof(struct APECSS_Bubble));
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_initializestruct(Bubbles[i]);

  /* Set default options and read the options for the bubble */
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_setdefaultoptions(Bubbles[i]);
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_readoptions(Bubbles[i], OptionsDir);

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
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
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
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->Excitation = Excitation;

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Create individual folders for the results of each bubble
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    if (Bubbles[i]->Results != NULL)
    {
      sprintf(Bubbles[i]->Results->dir, "./Bubble_%i/", RankInfo->bubblerank[RankInfo->rank] + i);
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
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_processoptions(Bubbles[i]);

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Use the revised pressure at infinity, including neighbor contributions
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->get_pressure_infinity = parallel_interactions_bubble_pressure_infinity;
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
    Bubbles[i]->get_pressurederivative_infinity = parallel_interactions_bubble_pressurederivative_infinity;

  // Allocate interaction structure
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    Bubbles[i]->Interaction = (struct APECSS_Interaction *) malloc(sizeof(struct APECSS_Interaction));
    Bubbles[i]->Interaction->dp_neighbor = 0.0;
    Bubbles[i]->Interaction->last_t_1 = 0.0;
    Bubbles[i]->Interaction->last_t_2 = 0.0;
    Bubbles[i]->Interaction->last_p_1 = 0.0;
    Bubbles[i]->Interaction->last_p_2 = 0.0;
  }

  // Update interaction structure
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->Interaction->nBubbles = nBubbles;

  // Define the size of each bubble
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    Bubbles[i]->R0 = 2.0e-6;
    Bubbles[i]->R = Bubbles[i]->R0;
  }

  // Define center location for each bubble
  APECSS_FLOAT Bubble_Center[nBubbles][3];

  for (register int i = 0; i < nBubbles; i++)
  {
    if (i == 0)
    {
      Bubble_Center[i][0] = 0.0;
      Bubble_Center[i][1] = 0.0;
      Bubble_Center[i][2] = 0.0;
    }
    else
    {
      APECSS_FLOAT x = rand_range((double) -cluster_radius, (double) cluster_radius);
      APECSS_FLOAT y = rand_range((double) -cluster_radius, (double) cluster_radius);
      APECSS_FLOAT z = rand_range((double) -cluster_radius, (double) cluster_radius);

      APECSS_FLOAT radius = APECSS_SQRT(APECSS_POW2(x) + APECSS_POW2(y) + APECSS_POW2(z));

      while (radius > cluster_radius)
      {
        x = rand_range((double) -cluster_radius, (double) cluster_radius);
        y = rand_range((double) -cluster_radius, (double) cluster_radius);
        z = rand_range((double) -cluster_radius, (double) cluster_radius);

        radius = APECSS_SQRT(APECSS_POW2(x) + APECSS_POW2(y) + APECSS_POW2(z));
      }

      Bubble_Center[i][0] = x;
      Bubble_Center[i][1] = y;
      Bubble_Center[i][2] = z;

      for (register int k = 0; k < i; k++)
      {
        APECSS_FLOAT bubbledist = APECSS_SQRT(APECSS_POW2(Bubble_Center[i][0] - Bubble_Center[k][0]) + APECSS_POW2(Bubble_Center[i][1] - Bubble_Center[k][1]) +
                                              APECSS_POW2(Bubble_Center[i][2] - Bubble_Center[k][2]));
        if (bubbledist < bubble_bubble_dist)
        {
          i--;
          break;
        }
      }
    }
  }

  for (register int n = 0; n < RankInfo->nBubbles_local; n++)
  {
    Bubbles[n]->Interaction->location[0] = Bubble_Center[RankInfo->bubblerank[RankInfo->rank] + n][0];
    Bubbles[n]->Interaction->location[1] = Bubble_Center[RankInfo->bubblerank[RankInfo->rank] + n][1];
    Bubbles[n]->Interaction->location[2] = Bubble_Center[RankInfo->bubblerank[RankInfo->rank] + n][2];
  }

  // Share the location of each bubble with all ranks

  RankInfo->bubbleglobal_R = malloc(RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));
  RankInfo->bubbleglobal_x = malloc(RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));
  RankInfo->bubbleglobal_y = malloc(RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));
  RankInfo->bubbleglobal_z = malloc(RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));

  for (int r = 0; r < RankInfo->size; r++)
  {
    int temp = RankInfo->nBubbles_local;
    MPI_Bcast(&temp, 1, MPI_INT, r, MPI_COMM_WORLD);

    APECSS_FLOAT temp_array[4];

    for (int i = 0; i < temp; i++)
    {
      if (r == RankInfo->rank)
      {
        temp_array[0] = Bubbles[i]->Interaction->location[0];
        temp_array[1] = Bubbles[i]->Interaction->location[1];
        temp_array[2] = Bubbles[i]->Interaction->location[2];
        temp_array[3] = Bubbles[i]->R;
      }

      MPI_Bcast(temp_array, 4, MPI_DOUBLE, r, MPI_COMM_WORLD);

      RankInfo->bubbleglobal_x[RankInfo->bubblerank[r] + i] = temp_array[0];
      RankInfo->bubbleglobal_y[RankInfo->bubblerank[r] + i] = temp_array[1];
      RankInfo->bubbleglobal_z[RankInfo->bubblerank[r] + i] = temp_array[2];
      RankInfo->bubbleglobal_R[RankInfo->bubblerank[r] + i] = temp_array[3];
    }
  }

  // Update cut off distance for each bubble
  parallel_interactions_proper_cutoffdistance(Bubbles, RankInfo);
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  clock_t starttimebubble = clock();
  APECSS_FLOAT tSim = 0.0;  // Simulation time of the coupled system

  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    Bubbles[i]->tStart = tSim;
    Bubbles[i]->tEnd = (APECSS_FLOAT) tEnd;
    Bubbles[i]->dt = APECSS_MIN(1.0e-11, dt_interbubble);  // Initial time-step
  }

  /* Initialize the bubble based on the selected options */
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_initialize(Bubbles[i]);

  /* Initialize */
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_solver_initialize(Bubbles[i]);

  // Allocate the pressure contribution array
  RankInfo->sumGU_rank = malloc(2 * RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // Files to retrieve all valuable information for cavitation onset test case
  // Bubbles locations
  FILE *file_loc = fopen("bubble_loc.txt", "w");
  fprintf(file_loc, "#number x(m) y(m) z(m)\n");
  for (register int i = 0; i < RankInfo->nBubbles_global; i++)
  {
    fprintf(file_loc, "%d %e %e %e\n", i, RankInfo->bubbleglobal_x[i], RankInfo->bubbleglobal_y[i], RankInfo->bubbleglobal_z[i]);
  }
  fclose(file_loc);

  // Radius and pressure evolution for each bubble during computation
  FILE *file_tension;
  char file_name[APECSS_STRINGLENGTH_SPRINTF_LONG];
  sprintf(file_name, "tension_results_%d.txt", RankInfo->rank);
  file_tension = fopen(file_name, "w");

  fprintf(file_tension, "%d Bubbles p0(pa) %e p1(Pa) %e\n", nBubbles, Liquid->pref, pa);
  fprintf(file_tension, "Initial_radii(m)");
  for (register int i = 0; i < RankInfo->nBubbles_global; i++)
  {
    fprintf(file_tension, " %e", RankInfo->bubbleglobal_R[i]);
  }
  fprintf(file_tension, "\n");
  fprintf(file_tension, "#Time(s) R(m) Pt(Pa)\n");
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /* Solve the bubble dynamics */
  while (tSim < (APECSS_FLOAT) tEnd)  // Interaction loop, corresponding to the time-intervals at which interactions are considered
  {
    APECSS_FLOAT dtSim = APECSS_MIN(dt_interbubble, (APECSS_FLOAT) tEnd - tSim);
    tSim += dtSim;

    for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_solver_run(tSim, Bubbles[i]);

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Update the contribution of the neighbor bubble
    parallel_interactions_quasi_acoustic(Bubbles, RankInfo);
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Retrieve data
    fprintf(file_tension, "%e", tSim);
    for (register int i = 0; i < RankInfo->nBubbles_local; i++)
    {
      fprintf(file_tension, " %e", Bubbles[i]->R);
    }
    for (register int i = 0; i < RankInfo->nBubbles_local; i++)
    {
      fprintf(file_tension, " %e", Bubbles[i]->get_pressure_infinity(Bubbles[i]->t, Bubbles[i]));
    }
    fprintf(file_tension, "\n");
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    for (register int i = 0; i < RankInfo->nBubbles_local; i++)
    {
      Bubbles[i]->Interaction->last_t_2 = Bubbles[i]->Interaction->last_t_1;
      Bubbles[i]->Interaction->last_p_2 = Bubbles[i]->Interaction->last_p_1;

      Bubbles[i]->Interaction->last_t_1 = tSim;
      Bubbles[i]->Interaction->last_p_1 = Bubbles[i]->Interaction->dp_neighbor;
    }
  }

  fclose(file_tension);

  /* Finalize the simulation*/
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_solver_finalize(Bubbles[i]);

  char str[APECSS_STRINGLENGTH_SPRINTF];
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    sprintf(str, "Bubble %i: Solver concluded %i time-steps and %i sub-iterations in %.3f s.", RankInfo->bubblerank[RankInfo->rank] + i, Bubbles[i]->dtNumber,
            Bubbles[i]->nSubIter, (double) (clock() - starttimebubble) / CLOCKS_PER_SEC);
    apecss_writeonscreen(str);
  }

  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_results_rayleighplesset_write(Bubbles[i], APECSS_RESULTS_WRITE);
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_results_emissionsspace_write(Bubbles[i], APECSS_RESULTS_WRITE);

  /* Make sure all allocated memory is freed */
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_freestruct(Bubbles[i]);

  for (register int i = 0; i < RankInfo->nBubbles_local; i++) free(Bubbles[i]);
  free(Gas);
  free(Liquid);
  free(Interface);
  free(NumericsODE);
  free(Excitation);

  free(RankInfo->bubblerank);
  RankInfo->bubblerank = NULL;
  free(RankInfo->bubbleglobal_R);
  RankInfo->bubbleglobal_R = NULL;
  free(RankInfo->bubbleglobal_x);
  RankInfo->bubbleglobal_x = NULL;
  free(RankInfo->bubbleglobal_y);
  RankInfo->bubbleglobal_y = NULL;
  free(RankInfo->bubbleglobal_z);
  RankInfo->bubbleglobal_z = NULL;
  free(RankInfo->sumGU_rank);
  RankInfo->sumGU_rank = NULL;
  free(RankInfo);

  MPI_Finalize();

  return (0);
}

APECSS_FLOAT parallel_interactions_bubble_pressure_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT tau = 1.75e-6;
  if (t < tau)
  {
    return (Bubble->p0 - (Bubble->p0 - Bubble->Excitation->dp) * APECSS_POW2(APECSS_SIN(APECSS_PI * t / tau)) + Bubble->Interaction->dp_neighbor);
  }
  else
  {
    return (Bubble->p0 + Bubble->Interaction->dp_neighbor);
  }
}

APECSS_FLOAT parallel_interactions_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  // Approximate numerical computation of p_infinity derivative
  APECSS_FLOAT tau = 1.75e-6;
  APECSS_FLOAT derivative = 0.0;
  if (t < tau)
  {
    APECSS_FLOAT inv_tau = 1 / tau;
    derivative = -2.0 * APECSS_PI * inv_tau * (Bubble->p0 - Bubble->Excitation->dp) * APECSS_COS(APECSS_PI * t * inv_tau) * APECSS_SIN(APECSS_PI * t * inv_tau);
  }
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

int parallel_interactions_quasi_acoustic(struct APECSS_Bubble *Bubbles[], struct APECSS_Parallel_Cluster *RankInfo)
{
  // All bubbles are supposed to be in the same liquid
  struct APECSS_Liquid *Liquid = Bubbles[0]->Liquid;

  // Update bubble radii info
  APECSS_FLOAT *tempR = malloc(RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));
  for (register int j = 0; j < RankInfo->nBubbles_global; j++) tempR[j] = 0.0;
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) tempR[RankInfo->bubblerank[RankInfo->rank] + i] = Bubbles[i]->R;
  MPI_Allreduce(tempR, RankInfo->bubbleglobal_R, RankInfo->nBubbles_global, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  free(tempR);

  // Reset pressure contributions of the neighours
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->Interaction->dp_neighbor = 0.0;

  // Locally compute the contribution to each bubble
  for (register int j = 0; j < RankInfo->nBubbles_global; j++)
  {
    APECSS_FLOAT x_j = RankInfo->bubbleglobal_x[j];
    APECSS_FLOAT y_j = RankInfo->bubbleglobal_y[j];
    APECSS_FLOAT z_j = RankInfo->bubbleglobal_z[j];

    RankInfo->sumGU_rank[j] = 0.0;
    RankInfo->sumGU_rank[j + RankInfo->nBubbles_global] = 0.0;

    double R_j = RankInfo->bubbleglobal_R[j];

    for (register int i = 0; i < RankInfo->nBubbles_local; i++)
    {
      if (j != RankInfo->bubblerank[RankInfo->rank] + i)
      {
        APECSS_FLOAT sumU_bubble = 0.0;
        APECSS_FLOAT sumG_bubble = 0.0;
        int nodecount_bubble = 0;

        APECSS_FLOAT x_i = Bubbles[i]->Interaction->location[0];
        APECSS_FLOAT y_i = Bubbles[i]->Interaction->location[1];
        APECSS_FLOAT z_i = Bubbles[i]->Interaction->location[2];

        APECSS_FLOAT interbubble_dist = APECSS_SQRT(APECSS_POW2(x_i - x_j) + APECSS_POW2(y_i - y_j) + APECSS_POW2(z_i - z_j));

        // The loop over the emission nodes begins with the most far away from the emitting bubble
        struct APECSS_EmissionNode *Current = Bubbles[i]->Emissions->LastNode;
        while (Current != NULL)
        {
          if (Current->r < interbubble_dist - R_j)
          {
            // The following emission nodes have not yet reached the bubble of interest
            Current = NULL;
          }
          else if (Current->r > interbubble_dist + R_j)
          {
            // The bubble of interest is still not reached
            Current = Current->backward;
          }
          else
          {
            // The current emission node is located inside the bubble of interest
            APECSS_FLOAT gr = Current->g / Current->r;
            sumU_bubble += (Current->f / APECSS_POW2(Current->r)) + gr / Liquid->cref;
            sumG_bubble += gr;
            nodecount_bubble++;
            Current = Current->backward;
          }
        }

        if (nodecount_bubble)
        {
          APECSS_FLOAT inv_nodecount = 1.0 / (APECSS_FLOAT) nodecount_bubble;
          RankInfo->sumGU_rank[j] += sumG_bubble * inv_nodecount;
          RankInfo->sumGU_rank[j + RankInfo->nBubbles_global] += sumU_bubble * inv_nodecount;
        }
      }
    }
  }

  APECSS_FLOAT *sumGU_all = malloc(2 * RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));
  MPI_Allreduce(RankInfo->sumGU_rank, sumGU_all, 2 * RankInfo->nBubbles_global, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Compute the total pressure contributions of the neighbours for all local bubbles
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    Bubbles[i]->Interaction->dp_neighbor =
        Liquid->rhoref * (sumGU_all[RankInfo->bubblerank[RankInfo->rank] + i] -
                          0.5 * APECSS_POW2(sumGU_all[RankInfo->nBubbles_global + RankInfo->bubblerank[RankInfo->rank] + i]));
  }

  free(sumGU_all);

  return (0);
}

int parallel_interactions_proper_cutoffdistance(struct APECSS_Bubble *Bubbles[], struct APECSS_Parallel_Cluster *RankInfo)
{
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    APECSS_FLOAT x_i = Bubbles[i]->Interaction->location[0];
    APECSS_FLOAT y_i = Bubbles[i]->Interaction->location[1];
    APECSS_FLOAT z_i = Bubbles[i]->Interaction->location[2];
    APECSS_FLOAT dist = 0.0;

    for (register int j = 0; j < RankInfo->nBubbles_global; j++)
    {
      APECSS_FLOAT x_j = RankInfo->bubbleglobal_x[j];
      APECSS_FLOAT y_j = RankInfo->bubbleglobal_y[j];
      APECSS_FLOAT z_j = RankInfo->bubbleglobal_z[j];

      APECSS_FLOAT dist_ij = APECSS_SQRT(APECSS_POW2(x_i - x_j) + APECSS_POW2(y_i - y_j) + APECSS_POW2(z_i - z_j));

      if (dist_ij > dist)
      {
        dist = dist_ij;
      }
    }
    if (Bubbles[i]->Emissions != NULL) Bubbles[i]->Emissions->CutOffDistance = 1.25 * dist;
  }
  return (0);
}