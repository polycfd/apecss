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
// demonstrating a simple parallelization of APECSS.
// -------------------------------------------------------------------

#include <time.h>
#include </usr/include/mpich-x86_64/mpi.h>
#include "apecss.h"

struct APECSS_Parallel_Cluster
{
  int rank, size;
  int nBubbles_local, nBubbles_global;

  int *bubblerank;  // max id of bubbles for each rank
  APECSS_FLOAT *bubbleglobal_x, *bubbleglobal_y, *bubbleglobal_z;  // x,y,z location of all bubbles
  APECSS_FLOAT *sumGU_rank;  // Contributions from local bubbles to all bubbles
};

// Declaration of additional case-dependent functions
APECSS_FLOAT parallel_interactions_bubble_pressure_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
int parallel_interactions_quasi_acoustic(struct APECSS_Bubble *Bubble[], struct APECSS_Parallel_Cluster *RankInfo);

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
  const int nBubbles_x = 5;
  const int nBubbles = nBubbles_x * nBubbles_x;  // Number of bubbles
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

  /* Allocate structure for parallel data */
  struct APECSS_Parallel_Cluster *RankInfo = (struct APECSS_Parallel_Cluster *) malloc(sizeof(struct APECSS_Parallel_Cluster));
  RankInfo->rank = mpi_rank;
  RankInfo->size = mpi_size;

  /* Determine the number of bubbles per rank */
  int max_per_rank = ceil((double) nBubbles / (double) mpi_size);
  RankInfo->nBubbles_local = APECSS_MIN(max_per_rank, nBubbles - mpi_rank * max_per_rank);

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

  // Allocate interaction structure
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->Interaction = (struct APECSS_Interaction *) malloc(sizeof(struct APECSS_Interaction));

  // Update interaction structure
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->Interaction->nBubbles = nBubbles;  // Not used?

  // Define the size of each bubble
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    if ((i + mpi_rank * max_per_rank) % 2)
      Bubbles[i]->R0 = 5.0e-6;
    else
      Bubbles[i]->R0 = 10.0e-6;

    Bubbles[i]->r_hc = Bubbles[i]->R0 / 8.54;
  }

  // Define center location for each bubble

  int ri = 0, rj = 0;
  if (mpi_rank)
  {
    ri = mpi_rank * max_per_rank % nBubbles_x;
    rj = mpi_rank * max_per_rank / nBubbles_x;
  }

  for (register int n = 0; n < RankInfo->nBubbles_local; n++)
  {
    Bubbles[n]->Interaction->location[0] = ((APECSS_FLOAT) ri) * bubble_bubble_dist;
    Bubbles[n]->Interaction->location[1] = ((APECSS_FLOAT) rj) * bubble_bubble_dist;
    Bubbles[n]->Interaction->location[2] = 0.0;

    if (ri < nBubbles_x - 1)
    {
      ri++;
    }
    else
    {
      ri = 0;
      rj++;
    }
  }

  // Share the location of each bubble with all ranks

  RankInfo->bubbleglobal_x = malloc(RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));
  RankInfo->bubbleglobal_y = malloc(RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));
  RankInfo->bubbleglobal_z = malloc(RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));

  for (int r = 0; r < RankInfo->size; r++)
  {
    int temp = RankInfo->nBubbles_local;
    MPI_Bcast(&temp, 1, MPI_INT, r, MPI_COMM_WORLD);

    APECSS_FLOAT temp_array[3];

    for (int i = 0; i < temp; i++)
    {
      if (r == RankInfo->rank)
      {
        temp_array[0] = Bubbles[i]->Interaction->location[0];
        temp_array[1] = Bubbles[i]->Interaction->location[1];
        temp_array[2] = Bubbles[i]->Interaction->location[2];
      }

      MPI_Bcast(temp_array, 3, MPI_DOUBLE, r, MPI_COMM_WORLD);

      RankInfo->bubbleglobal_x[RankInfo->bubblerank[r] + i] = temp_array[0];
      RankInfo->bubbleglobal_y[RankInfo->bubblerank[r] + i] = temp_array[1];
      RankInfo->bubbleglobal_z[RankInfo->bubblerank[r] + i] = temp_array[2];
    }
  }

  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // apecss_erroronscreen(-1, "OUT");

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

  /* Solve the bubble dynamics */
  while (tSim < (APECSS_FLOAT) tEnd)  // Interaction loop, corresponding to the time-intervals at which interactions are considered
  {
    APECSS_FLOAT dtSim = APECSS_MIN(dt_interbubble, (APECSS_FLOAT) tEnd - tSim);
    tSim += dtSim;
    printf("%e\n", tSim);

    for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_solver_run(tSim, Bubbles[i]);

    for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->Interaction->dp_neighbor = 0.0;

    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // Update the contribution of the neighbor bubble
    parallel_interactions_quasi_acoustic(Bubbles, RankInfo);
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }

  /* Finalize the simulation*/
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) apecss_bubble_solver_finalize(Bubbles[i]);

  char str[APECSS_STRINGLENGTH_SPRINTF];
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    sprintf(str, "Bubble %i: Solver concluded %i time-steps and %i sub-iterations in %.3f s.", i, Bubbles[i]->dtNumber, Bubbles[i]->nSubIter,
            (double) (clock() - starttimebubble) / CLOCKS_PER_SEC);
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

  free(RankInfo->sumGU_rank);
  RankInfo->sumGU_rank = NULL;
  free(RankInfo);

  MPI_Finalize();

  return (0);
}

APECSS_FLOAT parallel_interactions_bubble_pressure_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->p0 - Bubble->Excitation->dp * APECSS_SIN(2.0 * APECSS_PI * Bubble->Excitation->f * t) + Bubble->Interaction->dp_neighbor);
}

int parallel_interactions_quasi_acoustic(struct APECSS_Bubble *Bubbles[], struct APECSS_Parallel_Cluster *RankInfo)
{
  // All bubbles are supposed to be in the same liquid
  struct APECSS_Liquid *Liquid = Bubbles[0]->Liquid;

  // Reset pressure contributions of the neighours
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->Interaction->dp_neighbor = 0;

  // Locally compute the contribution to each bubble
  for (register int i = 0; i < RankInfo->nBubbles_global; i++)
  {
    APECSS_FLOAT x_i = RankInfo->bubbleglobal_x[i];
    APECSS_FLOAT y_i = RankInfo->bubbleglobal_y[i];
    APECSS_FLOAT z_i = RankInfo->bubbleglobal_z[i];

    RankInfo->sumGU_rank[i] = 0.0;
    RankInfo->sumGU_rank[i + RankInfo->nBubbles_global] = 0.0;

    for (register int j = 0; j < RankInfo->nBubbles_local; j++)
    {
      if (j != i)
      {
        APECSS_FLOAT sumU_bubble = 0.0;
        APECSS_FLOAT sumG_bubble = 0.0;
        int nodecount_bubble = 0;

        APECSS_FLOAT x_j = Bubbles[j]->Interaction->location[0];
        APECSS_FLOAT y_j = Bubbles[j]->Interaction->location[1];
        APECSS_FLOAT z_j = Bubbles[j]->Interaction->location[2];

        APECSS_FLOAT interbubble_dist = APECSS_SQRT(APECSS_POW2(x_i - x_j) + APECSS_POW2(y_i - y_j) + APECSS_POW2(z_i - z_j));

        double r_b = Bubbles[i]->R;

        // The loop over the emission nodes begins with the most far away from the emitting bubble
        struct APECSS_EmissionNode *Current = Bubbles[j]->Emissions->LastNode;
        while (Current != NULL)
        {
          if (Current->r < interbubble_dist - r_b)
            // The following emission nodes have not yet reached the bubble of interest
            Current = NULL;
          else if (Current->r > interbubble_dist + r_b)
            // The bubble of interest is still not reached
            Current = Current->backward;
          else
          {
            // The current emission node is located inside the bubble of interest
            APECSS_FLOAT gr = Current->g / Current->r;
            sumU_bubble += (Current->f / APECSS_POW2(Current->r)) + (gr) / Liquid->cref;
            sumG_bubble += gr;
            nodecount_bubble++;
            Current = Current->backward;
          }
        }

        if (nodecount_bubble)
        {
          APECSS_FLOAT inv_nodecount = 1.0 / (APECSS_FLOAT) nodecount_bubble;
          RankInfo->sumGU_rank[i] += sumG_bubble * inv_nodecount;
          RankInfo->sumGU_rank[i + RankInfo->nBubbles_global] += sumU_bubble * inv_nodecount;
        }
      }
    }
  }

  APECSS_FLOAT *sumGU_all = malloc(2 * RankInfo->nBubbles_global * sizeof(APECSS_FLOAT));
  for (register int i = 0; i < 2 * RankInfo->nBubbles_local; i++) sumGU_all[i] = 0.0;

  // Share the contributions from each rank to all bubbles
  for (int r = 0; r < RankInfo->size; r++)
  {
    APECSS_FLOAT *temp = RankInfo->sumGU_rank;
    MPI_Bcast(temp, 2 * RankInfo->nBubbles_global, MPI_INT, r, MPI_COMM_WORLD);
    for (register int i = 0; i < 2 * RankInfo->nBubbles_global; i++) sumGU_all[i] += temp[i];
  }

  // Compute the total pressure contributions of the neighbours for all local bubbles
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
    Bubbles[i]->Interaction->dp_neighbor = Liquid->rhoref * (sumGU_all[i] - 0.5 * APECSS_POW2(sumGU_all[RankInfo->nBubbles_global]));

  free(sumGU_all);

  return (0);
}