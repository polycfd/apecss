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
#include <mpi.h>
#include "apecss.h"

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
int apecss_bubble_new_solver_run(APECSS_FLOAT tend, APECSS_FLOAT tEnd, struct APECSS_Bubble *Bubble);

APECSS_FLOAT maxR;
APECSS_FLOAT minR;

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
  for (register int i = 0; i < RankInfo->nBubbles_local; i++) Bubbles[i]->Interaction->nBubbles = nBubbles;  // Not used?

  // Define the size of each bubble
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    Bubbles[i]->R0 = 1.0e-6;
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
    Bubbles[n]->Interaction->location[0] = ((APECSS_FLOAT) ri) * bubble_bubble_dist - (0.5 * (nBubbles_x - 1)) * bubble_bubble_dist;
    Bubbles[n]->Interaction->location[1] = ((APECSS_FLOAT) rj) * bubble_bubble_dist - (0.5 * (nBubbles_x - 1)) * bubble_bubble_dist;
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
  // Lists to gather maximum and minimum radius achieved by each bubble in steady state
  APECSS_FLOAT max_radii[RankInfo->nBubbles_local];
  APECSS_FLOAT min_radii[RankInfo->nBubbles_local];

  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    max_radii[i] = Bubbles[i]->R;
    min_radii[i] = Bubbles[i]->R;
  }

  // File to retrieve the maximum and minimum radius evolution achieved in steady state
  FILE *file_extremum;
  char file_name[APECSS_STRINGLENGTH_SPRINTF_LONG];
  sprintf(file_name, "bubblyscreen_extremum_%d.txt", RankInfo->rank);
  file_extremum = fopen(file_name, "w");

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

    // See progress during computation
    if (RankInfo->rank == 0)
    {
      printf("D:%le, f:%le, time:%e\n", bubble_bubble_dist, fa, tSim);
    }

    for (register int i = 0; i < RankInfo->nBubbles_local; i++)
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
    parallel_interactions_quasi_acoustic(Bubbles, RankInfo);
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if (RankInfo->rank == 0)
    {
      fprintf(file_radii, "%e", tSim);
      for (register int i = 0; i < RankInfo->nBubbles_global; i++)
      {
        fprintf(file_radii, " %e", RankInfo->bubbleglobal_R[i]);
      }
      fprintf(file_radii, "\n");
    }

    for (register int i = 0; i < RankInfo->nBubbles_local; i++)
    {
      Bubbles[i]->Interaction->last_t_2 = Bubbles[i]->Interaction->last_t_1;
      Bubbles[i]->Interaction->last_p_2 = Bubbles[i]->Interaction->last_p_1;

      Bubbles[i]->Interaction->last_t_1 = tSim;
      Bubbles[i]->Interaction->last_p_1 = Bubbles[i]->Interaction->dp_neighbor;
    }
  }

  /* Complete results file */
  for (register int i = 0; i < RankInfo->nBubbles_local; i++)
  {
    fprintf(file_extremum, "%d %e %e %e %e %e %e\n", RankInfo->bubblerank[RankInfo->rank] + i, Bubbles[i]->Interaction->location[0],
            Bubbles[i]->Interaction->location[1], Bubbles[i]->Interaction->location[2], Bubbles[i]->R0, min_radii[i], max_radii[i]);
  }
  fclose(file_extremum);

  fclose(file_radii);

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
  return (Bubble->p0 - Bubble->Excitation->dp * APECSS_SIN(2.0 * APECSS_PI * Bubble->Excitation->f * t) + Bubble->Interaction->dp_neighbor);
}

APECSS_FLOAT parallel_interactions_bubble_pressurederivative_infinity(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT derivative = -2.0 * APECSS_PI * Bubble->Excitation->f * Bubble->Excitation->dp * APECSS_COS(2.0 * APECSS_PI * Bubble->Excitation->f * t);
  APECSS_FLOAT delta_t = Bubble->Interaction->last_t_1 - Bubble->Interaction->last_t_2;
  if (delta_t > Bubble->dt)
  {
    return (derivative + (Bubble->Interaction->last_p_1 - Bubble->Interaction->last_p_2) / delta_t);
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