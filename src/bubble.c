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

#include "apecss.h"

// -------------------------------------------------------------------
// SET OPTIONS
// -------------------------------------------------------------------

int apecss_bubble_setdefaultoptions(struct APECSS_Bubble *Bubble)
{
  // Governing RP model
  Bubble->RPModel = APECSS_BUBBLEMODEL_RP;

  // Ambient conditions
  Bubble->p0 = 1.0e5;
  Bubble->T0 = 293.15;

  // Initial gas pressure in the bubble
  Bubble->pG0 = -1.0e10;

  // ODEs
  Bubble->dt = 1.0e-10;
  Bubble->nODEs = 2;  // ODEs for the bubble radius and bubble wall velocity
  Bubble->nUserODEs = 0;  // Number of ODEs defined by the user

  Bubble->NumericsODE = (struct APECSS_NumericsODE *) malloc(sizeof(struct APECSS_NumericsODE));
  Bubble->NumericsODE->RKtype = APECSS_RK54_7M;
  Bubble->NumericsODE->tol = 1.0e-10;
  Bubble->NumericsODE->maxSubIter = 20;
  Bubble->NumericsODE->dtMin = 1.0e-13;
  Bubble->NumericsODE->dtMax = 1.0e-6;
  Bubble->NumericsODE->minScale = 0.5;
  Bubble->NumericsODE->maxScale = 2.0;
  Bubble->NumericsODE->control_coeff_alpha = 0.9;
  Bubble->NumericsODE->control_coeff_q = 0.1;

  // Initialize pointers to optional structures
  Bubble->Excitation = NULL;
  Bubble->Emissions = NULL;
  Bubble->Results = NULL;

  // Initialize pointers to functions
  Bubble->emissions_initialize = apecss_emissions_initializenone;
  Bubble->emissions_update = apecss_emissions_updatenone;
  Bubble->emissions_free = apecss_emissions_freenone;
  Bubble->progress_initial = apecss_bubble_solverprogress_initialnone;
  Bubble->progress_update = apecss_bubble_solverprogress_updatenone;
  Bubble->progress_final = apecss_bubble_solverprogress_finalnone;
  Bubble->results_rayleighplesset_store = apecss_results_rayleighplesset_storenone;
  Bubble->results_emissionstime_write = apecss_results_emissionstime_writenone;
  Bubble->results_emissionstime_check = apecss_results_emissionstime_checknone;
  Bubble->results_emissionsspace_store = apecss_results_emissionsspace_storenone;
  Bubble->results_emissionsnode_alloc = apecss_results_emissionsnode_allocnone;
  Bubble->results_emissionsnode_store = apecss_results_emissionsnode_storenone;
  Bubble->results_emissionsnodeminmax_identify = apecss_results_emissionsnodeminmax_identifynone;

  return 0;
}

int apecss_bubble_processoptions(struct APECSS_Bubble *Bubble)
{
  // ODEs to describe viscoelasticity
  if (Bubble->Liquid->Type == APECSS_LIQUID_ZENER || Bubble->Liquid->Type == APECSS_LIQUID_OLDROYDB) Bubble->nODEs += 2;

  // Allocate the solution array and the array for the function pointers of the ODEs
  Bubble->ode = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->ODEsSol = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  for (register int i = 0; i < Bubble->nODEs; i++) Bubble->ODEsSol[i] = 0.0;

  Bubble->k2 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k3 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k4 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k5 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k6 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k7 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->kLast = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));

  // ---------------------------------------
  // RP model

  if (Bubble->RPModel & APECSS_BUBBLEMODEL_GILMORE)
    Bubble->ode[0] = apecss_rp_gilmorevelocity_ode;
  else if (Bubble->RPModel & APECSS_BUBBLEMODEL_KELLERMIKSIS)
    Bubble->ode[0] = apecss_rp_kellermiksisvelocity_ode;
  else if (Bubble->RPModel & APECSS_BUBBLEMODEL_RP_ACOUSTICRADIATION)
    Bubble->ode[0] = apecss_rp_rayleighplessetacousticrationvelocity_ode;
  else
    Bubble->ode[0] = apecss_rp_rayleighplessetvelocity_ode;

  Bubble->ode[1] = apecss_rp_bubbleradius_ode;

  // Bubble->nODEs has to reflect the number of actually defined ODEs, e.g. so that
  // additional user-defined ODEs can easily be added.
  Bubble->nODEs = 2;

  // ---------------------------------------
  // Viscoelasticity

  if (Bubble->Liquid->Type == APECSS_LIQUID_ZENER)
  {
    Bubble->ode[2] = apecss_viscoelastic_zenervarsigma_ode;
    Bubble->ode[3] = apecss_viscoelastic_zenertaurr_ode;
    Bubble->nODEs = 4;
  }
  else if (Bubble->Liquid->Type == APECSS_LIQUID_OLDROYDB)
  {
    Bubble->ode[2] = apecss_viscoelastic_oldroydb1_ode;
    Bubble->ode[3] = apecss_viscoelastic_oldroydb2_ode;
    Bubble->nODEs = 4;
  }

  // ---------------------------------------
  // ODE solver

  apecss_odesolver_rungekuttacoeffs(Bubble->NumericsODE, Bubble->nODEs);

  // ---------------------------------------
  // Emissions

  if (Bubble->Emissions != NULL)
  {
    Bubble->emissions_initialize = apecss_emissions_initializelinkedlist;
    Bubble->emissions_update = apecss_emissions_updatelinkedlist;
    Bubble->emissions_free = apecss_emissions_freelinkedlist;

    if (Bubble->Emissions->Type == APECSS_EMISSION_FINITE_TIME_INCOMPRESSIBLE)
    {
      Bubble->Emissions->advance = apecss_emissions_advance_finitetimeincompressible;
      Bubble->Emissions->get_advectingvelocity = apecss_emissions_getadvectingvelocity_returnzero;
    }
    else if (Bubble->Emissions->Type == APECSS_EMISSION_QUASIACOUSTIC)
    {
      Bubble->Emissions->advance = apecss_emissions_advance_quasiacoustic;
      Bubble->Emissions->get_advectingvelocity = apecss_emissions_getadvectingvelocity_returnzero;
    }
    else if (Bubble->Emissions->Type == APECSS_EMISSION_KIRKWOODBETHE)
    {
      Bubble->Emissions->advance = apecss_emissions_advance_kirkwoodbethe;
      Bubble->Emissions->get_advectingvelocity = apecss_emissions_getadvectingvelocity_returnvelocity;
    }
  }
  else
  {
    Bubble->emissions_initialize = apecss_emissions_initializenone;
    Bubble->emissions_update = apecss_emissions_updatenone;
    Bubble->emissions_free = apecss_emissions_freenone;
  }

  // ---------------------------------------
  // Results

  if (Bubble->Results != NULL)
  {
    if (Bubble->Results->RayleighPlesset != NULL && Bubble->Results->RayleighPlesset->freq > 0)
      Bubble->results_rayleighplesset_store = apecss_results_rayleighplesset_storeall;
    else
      Bubble->results_rayleighplesset_store = apecss_results_rayleighplesset_storenone;

    if (Bubble->Emissions != NULL && Bubble->Results->Emissions != NULL)
    {
      if (Bubble->Results->Emissions->nTimeInstances)
      {
        Bubble->results_emissionstime_write = apecss_results_emissionstime_writeall;
        Bubble->results_emissionstime_check = apecss_results_emissionstime_checktime;
      }
      else
      {
        Bubble->results_emissionstime_write = apecss_results_emissionstime_writenone;
        Bubble->results_emissionstime_check = apecss_results_emissionstime_checknone;
      }

      if (Bubble->Results->Emissions->nSpaceLocations)
        Bubble->results_emissionsspace_store = apecss_results_emissionsspace_storeall;
      else
        Bubble->results_emissionsspace_store = apecss_results_emissionsspace_storenone;

      Bubble->results_emissionsnode_alloc = apecss_results_emissionsnode_allocnone;  // default
      Bubble->results_emissionsnode_store = apecss_results_emissionsnode_storenone;  // default

      if (Bubble->Results->Emissions->nNodes)
      {
        Bubble->results_emissionsnode_store = apecss_results_emissionsnode_storeall;
        Bubble->Results->Emissions->store_nodespecific = apecss_results_emissionsnodespecific_storeall;

        Bubble->results_emissionsnode_alloc = apecss_results_emissionsnode_allocall;
        Bubble->Results->Emissions->allocation_nodespecific = apecss_results_emissionsnodespecific_allocall;
      }
      else
      {
        Bubble->Results->Emissions->store_nodespecific = apecss_results_emissionsnodespecific_storenone;
        Bubble->Results->Emissions->allocation_nodespecific = apecss_results_emissionsnodespecific_allocnone;
      }

      if (Bubble->Results->Emissions->MinMaxPeriod)
      {
        Bubble->results_emissionsnode_store = apecss_results_emissionsnode_storeall;
        Bubble->Results->Emissions->store_nodeminmax = apecss_results_emissionsnodeminmax_storeall;

        Bubble->results_emissionsnode_alloc = apecss_results_emissionsnode_allocall;
        Bubble->Results->Emissions->allocation_nodeminmax = apecss_results_emissionsnodeminmax_allocall;

        Bubble->results_emissionsnodeminmax_identify = apecss_results_emissionsnodeminmax_identifyall;
      }
      else
      {
        Bubble->results_emissionsnodeminmax_identify = apecss_results_emissionsnodeminmax_identifynone;
        Bubble->Results->Emissions->allocation_nodeminmax = apecss_results_emissionsnodeminmax_allocnone;
        Bubble->Results->Emissions->store_nodeminmax = apecss_results_emissionsnodeminmax_storenone;
      }
    }
    else
    {
      Bubble->results_emissionstime_write = apecss_results_emissionstime_writenone;
      Bubble->results_emissionstime_check = apecss_results_emissionstime_checknone;
      Bubble->results_emissionsspace_store = apecss_results_emissionsspace_storenone;
      Bubble->results_emissionsnode_alloc = apecss_results_emissionsnode_allocnone;
      Bubble->results_emissionsnode_store = apecss_results_emissionsnode_storenone;
      Bubble->results_emissionsnodeminmax_identify = apecss_results_emissionsnodeminmax_identifynone;
    }
  }
  else
  {
    Bubble->results_emissionstime_write = apecss_results_emissionstime_writenone;
    Bubble->results_emissionstime_check = apecss_results_emissionstime_checknone;
    Bubble->results_rayleighplesset_store = apecss_results_rayleighplesset_storenone;
    Bubble->results_emissionsspace_store = apecss_results_emissionsspace_storenone;
    Bubble->results_emissionsnode_alloc = apecss_results_emissionsnode_allocnone;
    Bubble->results_emissionsnode_store = apecss_results_emissionsnode_storenone;
    Bubble->results_emissionsnodeminmax_identify = apecss_results_emissionsnodeminmax_identifynone;
  }

  return 0;
}

int apecss_bubble_initialize(struct APECSS_Bubble *Bubble)
{
  for (register int i = 0; i < Bubble->nODEs; i++) Bubble->ODEsSol[i] = 0.0;

  // ---------------------------------------
  // Gas

  if (Bubble->Gas->EOS & APECSS_GAS_NASG)
  {
    if (Bubble->Gas->mmol > 0.0 && Bubble->Gas->dmol > 0.0)
    {
      Bubble->Gas->b = APECSS_AVOGADRO * 16.0 * APECSS_PI * APECSS_POW3(0.5 * Bubble->Gas->dmol) / (3.0 * Bubble->Gas->mmol);
    }
    else if (Bubble->Gas->b < 0.0)
    {
      char str[APECSS_STRINGLENGTH_SPRINTF];
      sprintf(str, "Co-volume of the gas cannot be computed. Insufficient data.");
      apecss_erroronscreen(-1, str);
    }

    Bubble->Gas->Kref =
        Bubble->Gas->rhoref / (APECSS_POW(Bubble->Gas->pref + Bubble->Gas->B, 1.0 / Bubble->Gas->Gamma) * (1.0 - Bubble->Gas->b * Bubble->Gas->rhoref));
  }
  else if (Bubble->Gas->EOS & APECSS_GAS_HC)
  {
    if (Bubble->Gas->mmol > 0.0 && Bubble->Gas->dmol > 0.0)
    {
      APECSS_FLOAT mgref = Bubble->Gas->rhoref * 4.0 * APECSS_PI * APECSS_POW3(Bubble->R0) / 3.0;
      Bubble->Gas->h = APECSS_POW(APECSS_AVOGADRO * mgref * 4.0 * APECSS_POW3(0.5 * Bubble->Gas->dmol) / (Bubble->Gas->mmol), APECSS_ONETHIRD);
    }
    else if (Bubble->Gas->h < 0.0)
    {
      char str[APECSS_STRINGLENGTH_SPRINTF];
      sprintf(str, "Hardcore radius of the gas cannot be computed. Insufficient data.");
      apecss_erroronscreen(-1, str);
    }

    Bubble->Gas->Kref = Bubble->Gas->rhoref / (APECSS_POW(Bubble->Gas->pref, 1.0 / Bubble->Gas->Gamma));
  }
  else
  {
    Bubble->Gas->b = 0.0;
    Bubble->Gas->B = 0.0;
    Bubble->Gas->h = 0.0;
    Bubble->Gas->Kref = Bubble->Gas->rhoref / (APECSS_POW(Bubble->Gas->pref, 1.0 / Bubble->Gas->Gamma));
  }

  if (Bubble->Gas->Kref <= 0.0)
  {
    char str[APECSS_STRINGLENGTH_SPRINTF];
    sprintf(str, "K reference factor of the gas is unphysical (Kref <= 0).");
    apecss_erroronscreen(-1, str);
  }

  if (Bubble->Interface->LipidCoatingModel & APECSS_LIPIDCOATING_MARMOTTANT)
  {
    if (Bubble->pG0 < -Bubble->Gas->B) Bubble->pG0 = Bubble->p0 + 2.0 * Bubble->Interface->sigma0 / Bubble->R0;
    Bubble->Interface->Rbuck = Bubble->R0 / APECSS_SQRT(1.0 + Bubble->Interface->sigma0 / Bubble->Interface->Elasticity);
    Bubble->Interface->Rrupt = Bubble->Interface->Rbuck * APECSS_SQRT(1.0 + Bubble->Interface->sigma / Bubble->Interface->Elasticity);

    if (Bubble->Interface->LipidCoatingModel & APECSS_LIPIDCOATING_GOMPERTZFUNCTION)
    {
      Bubble->Interface->GompertzC = 2.0 * Bubble->Interface->Elasticity * APECSS_E *
                                     APECSS_SQRT(1.0 + Bubble->Interface->sigma * 0.5 / Bubble->Interface->Elasticity) / Bubble->Interface->sigma;
      Bubble->Interface->GompertzB = -APECSS_LOG(Bubble->Interface->sigma0 / Bubble->Interface->sigma) /
                                     APECSS_EXP(Bubble->Interface->GompertzC * (1.0 - Bubble->R0 / Bubble->Interface->Rbuck));
    }
  }
  else
  {
    if (Bubble->pG0 < -Bubble->Gas->B) Bubble->pG0 = Bubble->p0 + 2.0 * Bubble->Interface->sigma / Bubble->R0;
  }

  // ---------------------------------------
  // Initial gas pressure assumed to result from polytropic expansion/compression from the reference state
  if (Bubble->Gas->EOS & APECSS_GAS_NASG)
  {
    APECSS_FLOAT peff = APECSS_POW((Bubble->pG0 + Bubble->Gas->B), (1.0 / Bubble->Gas->Gamma));
    Bubble->rhoG0 = Bubble->Gas->Kref * peff / (1.0 + Bubble->Gas->b * Bubble->Gas->Kref * peff);
  }
  else
  {
    Bubble->rhoG0 = Bubble->Gas->rhoref * APECSS_POW((Bubble->pG0 / Bubble->Gas->pref), (1.0 / Bubble->Gas->Gamma));
  }

  // ---------------------------------------
  // Solution

  for (register int i = 0; i < Bubble->nODEs; i++)
  {
    Bubble->k2[i] = 0.0;
    Bubble->k3[i] = 0.0;
    Bubble->k4[i] = 0.0;
    Bubble->k5[i] = 0.0;
    Bubble->k6[i] = 0.0;
    Bubble->k7[i] = 0.0;
    Bubble->kLast[i] = 0.0;
  }

  // ---------------------------------------
  // Bubble

  Bubble->t = Bubble->tStart;
  Bubble->R = Bubble->R0;
  Bubble->U = Bubble->U0;

  // Set starting values for the bubble radius and bubble wall velocity
  Bubble->ODEsSol[0] = Bubble->U;
  Bubble->ODEsSol[1] = Bubble->R;

  // ---------------------------------------
  // Emissions

  if (Bubble->Emissions != NULL) Bubble->Emissions->CutOffDistance = APECSS_MAX(Bubble->Emissions->CutOffDistance, Bubble->R0);

  // ---------------------------------------
  // Results

  if (Bubble->Results != NULL && Bubble->Results->RayleighPlesset != NULL) apecss_results_rayleighplesset_initialize(Bubble);

  return 0;
}

int apecss_bubble_freearrays(struct APECSS_Bubble *Bubble)
{
  free(Bubble->ODEsSol);
  free(Bubble->ode);

  free(Bubble->k2);
  free(Bubble->k3);
  free(Bubble->k4);
  free(Bubble->k5);
  free(Bubble->k6);
  free(Bubble->k7);
  free(Bubble->kLast);

  if (Bubble->Results != NULL && Bubble->Results->Emissions != NULL)
  {
    if (Bubble->Results->Emissions->nTimeInstances)
    {
      free(Bubble->Results->Emissions->TimeInstances);
      Bubble->Results->Emissions->TimeInstances = NULL;
    }

    if (Bubble->Results->Emissions->nSpaceLocations)
    {
      free(Bubble->Results->Emissions->SpaceLocation);
      Bubble->Results->Emissions->nSpaceLocations = 0;
    }

    if (Bubble->Results->Emissions->nNodes)
    {
      free(Bubble->Results->Emissions->Node);
      Bubble->Results->Emissions->nNodes = 0;
    }
  }

  return 0;
}

// -------------------------------------------------------------------
// SOLVE
// -------------------------------------------------------------------

int apecss_bubble_solve(struct APECSS_Bubble *Bubble)
{
  Bubble->dtNumber = 0;
  Bubble->nSubIter = 0;

  // Set start time
  Bubble->t = Bubble->tStart;

  // ---------------------------------------
  // Getting ready

  // Emissions
  Bubble->emissions_initialize(Bubble);

  // Store the initial results
  Bubble->results_rayleighplesset_store(Bubble);

  // Old solution required for sub-iterations
  APECSS_FLOAT *OldSol = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));

  // Initialize the error variable
  APECSS_FLOAT err = Bubble->NumericsODE->tol;

  // Initialise the last-step solution of the RK54 scheme
  for (register int i = 0; i < Bubble->nODEs; i++) Bubble->kLast[i] = Bubble->ode[i](Bubble->ODEsSol, Bubble->tStart, Bubble);

  // ---------------------------------------
  // Time loop

  Bubble->progress_initial();
  int prog = 0;

  while (Bubble->t < Bubble->tEnd - APECSS_EPS)
  {
    // Store the previous solution for sub-iterations
    for (register int i = 0; i < Bubble->nODEs; i++) OldSol[i] = Bubble->ODEsSol[i];

    // Check what comes next: the end of the simulation or, if applicable, a time instance to output the emissions
    APECSS_FLOAT next_event_time = Bubble->results_emissionstime_check(Bubble);

    // Set the time-step for the ODEs
    apecss_odesolver_settimestep(Bubble->NumericsODE, err, next_event_time - Bubble->t, &(*Bubble).dt);

    // Solve the ODEs
    err = apecss_odesolver(Bubble);

    // Perform sub-iterations on the control of dt when err > tol
    int attempts = 0;
    while ((err > Bubble->NumericsODE->tol) && (attempts < Bubble->NumericsODE->maxSubIter))
    {
      for (register int i = 0; i < Bubble->nODEs; i++) Bubble->ODEsSol[i] = OldSol[i];
      apecss_odesolver_settimestep(Bubble->NumericsODE, err, next_event_time - Bubble->t, &(*Bubble).dt);
      err = apecss_odesolver(Bubble);
      ++attempts;
    }
    Bubble->nSubIter += attempts;

    // Set new values
    ++(Bubble->dtNumber);
    Bubble->t += Bubble->dt;
    Bubble->U = Bubble->ODEsSol[0];
    Bubble->R = Bubble->ODEsSol[1];

    // Emissions
    Bubble->results_emissionsnodeminmax_identify(Bubble);
    Bubble->results_emissionsnode_alloc(Bubble);
    Bubble->emissions_update(Bubble);

    // Store results
    Bubble->results_rayleighplesset_store(Bubble);
    Bubble->results_emissionsspace_store(Bubble);
    Bubble->results_emissionstime_write(Bubble);

    // Update the last-step solution of the RK54 scheme
    for (register int i = 0; i < Bubble->nODEs; i++) Bubble->kLast[i] = Bubble->k7[i];

    // Update progress screen in the terminal
    Bubble->progress_update(&prog, Bubble->t, Bubble->tEnd - Bubble->tStart);
  }

  Bubble->progress_final();

  // ---------------------------------------
  // Clean up

  Bubble->emissions_free(Bubble);
  free(OldSol);

  return 0;
}

// -------------------------------------------------------------------
// PROGRESS SCREEN
// -------------------------------------------------------------------

int apecss_bubble_solverprogress_initialnone() { return 0; }

int apecss_bubble_solverprogress_initialscreen()
{
  fprintf(stderr, "| APECSS | Progress %%: ");
  return 0;
}

int apecss_bubble_solverprogress_updatenone(int *prog, APECSS_FLOAT t, APECSS_FLOAT totaltime) { return 0; }

int apecss_bubble_solverprogress_updatescreen(int *prog, APECSS_FLOAT t, APECSS_FLOAT totaltime)
{
  if (t > (APECSS_FLOAT) *prog * 0.02 * totaltime)
  {
    if (!(*prog % 5))
      fprintf(stderr, "%i", *prog * 2);
    else
      fprintf(stderr, ".");

    (*prog)++;
  }

  return 0;
}

int apecss_bubble_solverprogress_finalnone() { return 0; }

int apecss_bubble_solverprogress_finalscreen()
{
  fprintf(stderr, "100\n");
  return 0;
}