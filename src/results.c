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

#include "apecss.h"

// -------------------------------------------------------------------
// INITIALIZE
// -------------------------------------------------------------------
// Functions initializing the structures that hold the results.
// -------------------------------------------------------------------

int apecss_results_initializestruct(struct APECSS_Bubble *Bubble)
{
  // Initialize the results structure attached to the bubble
  Bubble->Results = (struct APECSS_Results *) malloc(sizeof(struct APECSS_Results));
  sprintf(Bubble->Results->dir, "./");
  Bubble->Results->digits = 6;
  Bubble->Results->RayleighPlesset = NULL;
  Bubble->Results->Emissions = NULL;

  return (0);
}

int apecss_results_rayleighplesset_initializestruct(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Results == NULL) apecss_results_initializestruct(Bubble);

  // Initialize the structure for the Rayleigh-Plesset results
  Bubble->Results->RayleighPlesset = (struct APECSS_ResultsBubble *) malloc(sizeof(struct APECSS_ResultsBubble));
  Bubble->Results->RayleighPlesset->freq = 1;
  Bubble->Results->RayleighPlesset->n = 0;
  Bubble->Results->RayleighPlesset->nAllocated = 0;
  Bubble->Results->RayleighPlesset->nUserODEs = 0;

  return (0);
}

int apecss_results_emissions_initializestruct(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Results == NULL) apecss_results_initializestruct(Bubble);

  // Initialize the structure for the results of the acoustic emissions
  Bubble->Results->Emissions = (struct APECSS_ResultsEmissions *) malloc(sizeof(struct APECSS_ResultsEmissions));
  Bubble->Results->Emissions->nTimeInstances = 0;
  Bubble->Results->Emissions->nextTimeInstance = 0;
  Bubble->Results->Emissions->TimeInstances = NULL;
  Bubble->Results->Emissions->freqSpaceLocations = 1;
  Bubble->Results->Emissions->nSpaceLocations = 0;
  Bubble->Results->Emissions->SpaceLocation = NULL;
  Bubble->Results->Emissions->nNodes = 0;
  Bubble->Results->Emissions->Node = NULL;
  Bubble->Results->Emissions->MinMaxPeriod = 0;
  Bubble->Results->Emissions->Node_Rmin = NULL;
  Bubble->Results->Emissions->Node_Umin = NULL;
  Bubble->Results->Emissions->Node_pLmax = NULL;

  return (0);
}

// -------------------------------------------------------------------
// BUBBLE
// -------------------------------------------------------------------
// Functions for the handling of the results associated with the
// bubble dynamics, e.g. bubble radius and bubble wall velocity.
// -------------------------------------------------------------------

int apecss_results_rayleighplesset_initialize(struct APECSS_Bubble *Bubble)
{
  Bubble->Results->RayleighPlesset->t = NULL;
  Bubble->Results->RayleighPlesset->dt = NULL;
  Bubble->Results->RayleighPlesset->R = NULL;
  Bubble->Results->RayleighPlesset->U = NULL;
  Bubble->Results->RayleighPlesset->pG = NULL;
  Bubble->Results->RayleighPlesset->pL = NULL;
  Bubble->Results->RayleighPlesset->pinf = NULL;
  Bubble->Results->RayleighPlesset->cL = NULL;

  for (register int userode = 0; userode < Bubble->Results->RayleighPlesset->nUserODEs; userode++)
    Bubble->Results->RayleighPlesset->UserODEsSol[userode] = NULL;

  Bubble->Results->RayleighPlesset->n = 0;
  Bubble->Results->RayleighPlesset->nAllocated = 0;

  return (0);
}

int apecss_results_rayleighplesset_storenone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_results_rayleighplesset_storeall(struct APECSS_Bubble *Bubble)
{
  if (!(Bubble->dtNumber % Bubble->Results->RayleighPlesset->freq))  // Check if it's time to write out the RP results
  {
    if (Bubble->Results->RayleighPlesset->n == Bubble->Results->RayleighPlesset->nAllocated)  // Extend the length of the result arrays
    {
      // Allocate a temporary array to hold the results
      APECSS_FLOAT *temp;
      temp = malloc(Bubble->Results->RayleighPlesset->n * sizeof(APECSS_FLOAT));

      // Define the new length of the result arrays
      Bubble->Results->RayleighPlesset->nAllocated =
          Bubble->Results->RayleighPlesset->n + APECSS_CEIL(APECSS_DATA_ALLOC_INCREMENT / Bubble->Results->RayleighPlesset->freq);

      // Extend the length of all arrays
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->t[i];
      free(Bubble->Results->RayleighPlesset->t);
      Bubble->Results->RayleighPlesset->t = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->t[i] = temp[i];

      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->dt[i];
      free(Bubble->Results->RayleighPlesset->dt);
      Bubble->Results->RayleighPlesset->dt = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->dt[i] = temp[i];

      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->R[i];
      free(Bubble->Results->RayleighPlesset->R);
      Bubble->Results->RayleighPlesset->R = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->R[i] = temp[i];

      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->U[i];
      free(Bubble->Results->RayleighPlesset->U);
      Bubble->Results->RayleighPlesset->U = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->U[i] = temp[i];

      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->pG[i];
      free(Bubble->Results->RayleighPlesset->pG);
      Bubble->Results->RayleighPlesset->pG = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->pG[i] = temp[i];

      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->pL[i];
      free(Bubble->Results->RayleighPlesset->pL);
      Bubble->Results->RayleighPlesset->pL = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->pL[i] = temp[i];

      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->pinf[i];
      free(Bubble->Results->RayleighPlesset->pinf);
      Bubble->Results->RayleighPlesset->pinf = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->pinf[i] = temp[i];

      if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE)
      {
        for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->cL[i];
        free(Bubble->Results->RayleighPlesset->cL);
        Bubble->Results->RayleighPlesset->cL = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
        for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->cL[i] = temp[i];
      }

      for (register int userode = 0; userode < Bubble->Results->RayleighPlesset->nUserODEs; userode++)
      {
        for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) temp[i] = Bubble->Results->RayleighPlesset->UserODEsSol[userode][i];
        free(Bubble->Results->RayleighPlesset->UserODEsSol[userode]);
        Bubble->Results->RayleighPlesset->UserODEsSol[userode] = malloc(Bubble->Results->RayleighPlesset->nAllocated * sizeof(APECSS_FLOAT));
        for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++) Bubble->Results->RayleighPlesset->UserODEsSol[userode][i] = temp[i];
      }

      free(temp);
    }
    else if (Bubble->Results->RayleighPlesset->n > Bubble->Results->RayleighPlesset->nAllocated)
    {
      char str[APECSS_STRINGLENGTH_SPRINTF];
      sprintf(str, "Bubble results array allocation invalid (Result number: %i, Allocated %i)", Bubble->Results->RayleighPlesset->n,
              Bubble->Results->RayleighPlesset->nAllocated);
      apecss_erroronscreen(-1, str);
    }

    Bubble->Results->RayleighPlesset->t[Bubble->Results->RayleighPlesset->n] = Bubble->t;
    Bubble->Results->RayleighPlesset->dt[Bubble->Results->RayleighPlesset->n] = Bubble->dt;
    Bubble->Results->RayleighPlesset->R[Bubble->Results->RayleighPlesset->n] = Bubble->R;
    Bubble->Results->RayleighPlesset->U[Bubble->Results->RayleighPlesset->n] = Bubble->U;
    Bubble->Results->RayleighPlesset->pG[Bubble->Results->RayleighPlesset->n] = Bubble->Gas->get_pressure(Bubble->ODEsSol, Bubble);

    APECSS_FLOAT pL = Bubble->Liquid->get_pressure_bubblewall(Bubble->ODEsSol, Bubble->t, Bubble);
    Bubble->Results->RayleighPlesset->pL[Bubble->Results->RayleighPlesset->n] = pL;

    Bubble->Results->RayleighPlesset->pinf[Bubble->Results->RayleighPlesset->n] = Bubble->get_pressure_infinity(Bubble->t, Bubble);

    if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE)
      Bubble->Results->RayleighPlesset->cL[Bubble->Results->RayleighPlesset->n] =
          Bubble->Liquid->get_soundspeed(pL, Bubble->Liquid->get_density(pL, Bubble->Liquid), Bubble->Liquid);

    for (register int userode = 0; userode < Bubble->Results->RayleighPlesset->nUserODEs; userode++)
    {
      Bubble->Results->RayleighPlesset->UserODEsSol[userode][Bubble->Results->RayleighPlesset->n] =
          Bubble->ODEsSol[(Bubble->nODEs - Bubble->nUserODEs) + userode];
    }

    ++(Bubble->Results->RayleighPlesset->n);
  }

  return (0);
}

int apecss_results_rayleighplesset_write(struct APECSS_Bubble *Bubble, int write)
{
  if (Bubble->Results->RayleighPlesset != NULL)
  {
    if (Bubble->Results->RayleighPlesset->n && write)
    {
      char path[APECSS_STRINGLENGTH_SPRINTF_LONG], rptype[APECSS_STRINGLENGTH];

      if (Bubble->RPModel == APECSS_BUBBLEMODEL_RP)
        sprintf(rptype, "RP");
      else if (Bubble->RPModel == APECSS_BUBBLEMODEL_RP_ACOUSTICRADIATION)
        sprintf(rptype, "RPAR");
      else if (Bubble->RPModel == APECSS_BUBBLEMODEL_KELLERMIKSIS)
        sprintf(rptype, "KellerMiksis");
      else if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE)
        sprintf(rptype, "Gilmore");

      if (Bubble->Excitation != NULL)
      {
#if defined(APECSS_PRECISION_LONGDOUBLE)
        sprintf(path, "%s/%s_R%.3Le_fa%.3Le_pa%.3Le.txt", Bubble->Results->dir, rptype, Bubble->R0, Bubble->Excitation->f, Bubble->Excitation->dp);
#else
        sprintf(path, "%s/%s_R%.3e_fa%.3e_pa%.3e.txt", Bubble->Results->dir, rptype, Bubble->R0, Bubble->Excitation->f, Bubble->Excitation->dp);
#endif
      }
      else
      {
#if defined(APECSS_PRECISION_LONGDOUBLE)
        sprintf(path, "%s/%s_R%.3Le.txt", Bubble->Results->dir, rptype, Bubble->R0);
#else
        sprintf(path, "%s/%s_R%.3e.txt", Bubble->Results->dir, rptype, Bubble->R0);
#endif
      }

      FILE *file_ptr;
      if (write == APECSS_RESULTS_WRITE)  // Create new file and write all data in one go.
      {
        file_ptr = fopen(path, "w");

        fprintf(file_ptr, "# timeStep time dt R U pG pL pinf");
        if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE) fprintf(file_ptr, " cL");
        for (register int userode = 0; userode < Bubble->Results->RayleighPlesset->nUserODEs; userode++)
          fprintf(file_ptr, " %s", Bubble->Results->RayleighPlesset->UserODEsName[userode]);
        fprintf(file_ptr, " \n");
      }
      else  // Append file.
      {
        if (!(file_ptr = fopen(path, "r")))  // File does not exist yet -> write header.
        {
          file_ptr = fopen(path, "a");

          fprintf(file_ptr, "# timeStep time dt R U pG pL pinf");
          if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE) fprintf(file_ptr, " cL");
          for (register int userode = 0; userode < Bubble->Results->RayleighPlesset->nUserODEs; userode++)
            fprintf(file_ptr, " %s", Bubble->Results->RayleighPlesset->UserODEsName[userode]);
          fprintf(file_ptr, " \n");
        }
        else  // File exists already, append.
        {
          fclose(file_ptr);
          file_ptr = fopen(path, "a");
        }
      }

#if defined(APECSS_PRECISION_LONGDOUBLE)
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++)
      {
        fprintf(file_ptr, "%d %.*Le %.*Le %.*Le %.*Le %.*Le %.*Le %.*Le", i * Bubble->Results->RayleighPlesset->freq, Bubble->Results->digits,
                Bubble->Results->RayleighPlesset->t[i], Bubble->Results->digits, Bubble->Results->RayleighPlesset->dt[i], Bubble->Results->digits,
                Bubble->Results->RayleighPlesset->R[i], Bubble->Results->digits, Bubble->Results->RayleighPlesset->U[i], Bubble->Results->digits,
                Bubble->Results->RayleighPlesset->pG[i], Bubble->Results->digits, Bubble->Results->RayleighPlesset->pL[i], Bubble->Results->digits,
                Bubble->Results->RayleighPlesset->pinf[i]);

        if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE) fprintf(file_ptr, " %.*Le", Bubble->Results->digits, Bubble->Results->RayleighPlesset->cL[i]);

        for (register int userode = 0; userode < Bubble->Results->RayleighPlesset->nUserODEs; userode++)
          fprintf(file_ptr, " %.*Le", Bubble->Results->digits, Bubble->Results->RayleighPlesset->UserODEsSol[userode][i]);

        fprintf(file_ptr, " \n");
      }
#else
      for (register int i = 0; i < Bubble->Results->RayleighPlesset->n; i++)
      {
        fprintf(file_ptr, "%d %.*e %.*e %.*e %.*e %.*e %.*e %.*e", i * Bubble->Results->RayleighPlesset->freq, Bubble->Results->digits,
                Bubble->Results->RayleighPlesset->t[i], Bubble->Results->digits, Bubble->Results->RayleighPlesset->dt[i], Bubble->Results->digits,
                Bubble->Results->RayleighPlesset->R[i], Bubble->Results->digits, Bubble->Results->RayleighPlesset->U[i], Bubble->Results->digits,
                Bubble->Results->RayleighPlesset->pG[i], Bubble->Results->digits, Bubble->Results->RayleighPlesset->pL[i], Bubble->Results->digits,
                Bubble->Results->RayleighPlesset->pinf[i]);

        if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE) fprintf(file_ptr, " %.*e", Bubble->Results->digits, Bubble->Results->RayleighPlesset->cL[i]);

        for (register int userode = 0; userode < Bubble->Results->RayleighPlesset->nUserODEs; userode++)
          fprintf(file_ptr, " %.*e", Bubble->Results->digits, Bubble->Results->RayleighPlesset->UserODEsSol[userode][i]);

        fprintf(file_ptr, " \n");
      }
#endif

      fclose(file_ptr);
    }

    apecss_results_rayleighplesset_free(Bubble);
  }

  return (0);
}

int apecss_results_rayleighplesset_free(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Results->RayleighPlesset->nAllocated)
  {
    free(Bubble->Results->RayleighPlesset->t);
    Bubble->Results->RayleighPlesset->t = NULL;
    free(Bubble->Results->RayleighPlesset->dt);
    Bubble->Results->RayleighPlesset->dt = NULL;
    free(Bubble->Results->RayleighPlesset->R);
    Bubble->Results->RayleighPlesset->R = NULL;
    free(Bubble->Results->RayleighPlesset->U);
    Bubble->Results->RayleighPlesset->U = NULL;
    free(Bubble->Results->RayleighPlesset->pG);
    Bubble->Results->RayleighPlesset->pG = NULL;
    free(Bubble->Results->RayleighPlesset->pL);
    Bubble->Results->RayleighPlesset->pL = NULL;
    free(Bubble->Results->RayleighPlesset->pinf);
    Bubble->Results->RayleighPlesset->pinf = NULL;

    if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE)
    {
      free(Bubble->Results->RayleighPlesset->cL);
      Bubble->Results->RayleighPlesset->cL = NULL;
    }

    for (register int userode = 0; userode < Bubble->Results->RayleighPlesset->nUserODEs; userode++)
    {
      free(Bubble->Results->RayleighPlesset->UserODEsSol[userode]);
      Bubble->Results->RayleighPlesset->UserODEsSol[userode] = NULL;
    }

    Bubble->Results->RayleighPlesset->n = 0;
    Bubble->Results->RayleighPlesset->nAllocated = 0;
  }

  return (0);
}

// -------------------------------------------------------------------
// EMISSIONS  (time instance)
// -------------------------------------------------------------------
// Functions handling the results associated with the acoustics
// emissions recorded in space at user-defined instances of time.
// -------------------------------------------------------------------

int apecss_results_emissionstime_writenone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_results_emissionstime_writeall(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Results->Emissions->nextTimeInstance < Bubble->Results->Emissions->nTimeInstances &&
      APECSS_ABS(Bubble->t - Bubble->Results->Emissions->TimeInstances[Bubble->Results->Emissions->nextTimeInstance]) < APECSS_SMALL)
  {
    char path[APECSS_STRINGLENGTH_SPRINTF_LONG];
    sprintf(path, "%s/EmissionsTime_%.*e.txt", Bubble->Results->dir, Bubble->Results->digits, (double) Bubble->t);

    FILE *file_ptr;
    file_ptr = fopen(path, "w");

    fprintf(file_ptr, "# real-id r p u c pinf \n");

    struct APECSS_EmissionNode *Node = Bubble->Emissions->FirstNode;

    while (Node != NULL)
    {
      fprintf(file_ptr, "%i", Node->id);
      fprintf(file_ptr, " ");

#if defined(APECSS_PRECISION_LONGDOUBLE)
      fprintf(file_ptr, "%.*Le %.*Le %.*Le %.*Le %.*Le \n", Bubble->Results->digits, Node->r, Bubble->Results->digits, Node->p, Bubble->Results->digits,
              Node->u, Bubble->Results->digits, Bubble->Liquid->get_soundspeed(Node->p, Bubble->Liquid->get_density(Node->p, Bubble->Liquid), Bubble->Liquid),
              Bubble->Results->digits, Bubble->get_pressure_infinity(Bubble->t, Bubble));
#else
      fprintf(file_ptr, "%.*e %.*e %.*e %.*e %.*e \n", Bubble->Results->digits, Node->r, Bubble->Results->digits, Node->p, Bubble->Results->digits, Node->u,
              Bubble->Results->digits, Bubble->Liquid->get_soundspeed(Node->p, Bubble->Liquid->get_density(Node->p, Bubble->Liquid), Bubble->Liquid),
              Bubble->Results->digits, Bubble->get_pressure_infinity(Bubble->t, Bubble));
#endif

      Node = Node->forward;
    }

    fclose(file_ptr);

    Bubble->Results->Emissions->nextTimeInstance++;
  }

  return (0);
}

APECSS_FLOAT apecss_results_emissionstime_checknone(APECSS_FLOAT tend, struct APECSS_Bubble *Bubble) { return (tend); }

APECSS_FLOAT apecss_results_emissionstime_checktime(APECSS_FLOAT tend, struct APECSS_Bubble *Bubble)
{
  // The last event is obviously the end of the simulation
  APECSS_FLOAT nexttime = tend;

  // If there are any more time instances to be written out, they need to be considered
  if (Bubble->Results->Emissions->nextTimeInstance < Bubble->Results->Emissions->nTimeInstances)
    nexttime = APECSS_MIN(nexttime, Bubble->Results->Emissions->TimeInstances[Bubble->Results->Emissions->nextTimeInstance]);

  return (nexttime);
}

// -------------------------------------------------------------------
// EMISSIONS  (space locations)
// -------------------------------------------------------------------
// Functions handling the results associated with the acoustics
// emissions recorded in time at user-defined locations in space.
// -------------------------------------------------------------------

int apecss_results_emissionsspace_storenone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_results_emissionsspace_storeall(struct APECSS_Bubble *Bubble)
{
  for (register int l = 0; l < Bubble->Results->Emissions->nSpaceLocations; l++)  // Loop over all radial locations
  {
    if (!(Bubble->dtNumber % Bubble->Results->Emissions->freqSpaceLocations))  // Check if it's time to write out the RP results
    {
      if (Bubble->Results->Emissions->SpaceLocation[l].n == Bubble->Results->Emissions->SpaceLocation[l].nAllocated)  // Extend the length of the result arrays
      {
        // Allocate a temporary array to hold the results
        APECSS_FLOAT *temp;
        temp = malloc(Bubble->Results->Emissions->SpaceLocation[l].n * sizeof(APECSS_FLOAT));

        // Define the new length of the result arrays
        Bubble->Results->Emissions->SpaceLocation[l].nAllocated =
            Bubble->Results->Emissions->SpaceLocation[l].n + APECSS_CEIL(APECSS_DATA_ALLOC_INCREMENT / Bubble->Results->Emissions->freqSpaceLocations);

        // Extend the length of all arrays
        for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) temp[i] = Bubble->Results->Emissions->SpaceLocation[l].t[i];
        free(Bubble->Results->Emissions->SpaceLocation[l].t);
        Bubble->Results->Emissions->SpaceLocation[l].t = malloc(Bubble->Results->Emissions->SpaceLocation[l].nAllocated * sizeof(APECSS_FLOAT));
        for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) Bubble->Results->Emissions->SpaceLocation[l].t[i] = temp[i];

        for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) temp[i] = Bubble->Results->Emissions->SpaceLocation[l].p[i];
        free(Bubble->Results->Emissions->SpaceLocation[l].p);
        Bubble->Results->Emissions->SpaceLocation[l].p = malloc(Bubble->Results->Emissions->SpaceLocation[l].nAllocated * sizeof(APECSS_FLOAT));
        for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) Bubble->Results->Emissions->SpaceLocation[l].p[i] = temp[i];

        for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) temp[i] = Bubble->Results->Emissions->SpaceLocation[l].u[i];
        free(Bubble->Results->Emissions->SpaceLocation[l].u);
        Bubble->Results->Emissions->SpaceLocation[l].u = malloc(Bubble->Results->Emissions->SpaceLocation[l].nAllocated * sizeof(APECSS_FLOAT));
        for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) Bubble->Results->Emissions->SpaceLocation[l].u[i] = temp[i];

        if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE)
        {
          for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) temp[i] = Bubble->Results->Emissions->SpaceLocation[l].c[i];
          free(Bubble->Results->Emissions->SpaceLocation[l].c);
          Bubble->Results->Emissions->SpaceLocation[l].c = malloc(Bubble->Results->Emissions->SpaceLocation[l].nAllocated * sizeof(APECSS_FLOAT));
          for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) Bubble->Results->Emissions->SpaceLocation[l].c[i] = temp[i];
        }

        for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) temp[i] = Bubble->Results->Emissions->SpaceLocation[l].pInf[i];
        free(Bubble->Results->Emissions->SpaceLocation[l].pInf);
        Bubble->Results->Emissions->SpaceLocation[l].pInf = malloc(Bubble->Results->Emissions->SpaceLocation[l].nAllocated * sizeof(APECSS_FLOAT));
        for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++) Bubble->Results->Emissions->SpaceLocation[l].pInf[i] = temp[i];

        free(temp);
      }
      else if (Bubble->Results->Emissions->SpaceLocation[l].n > Bubble->Results->Emissions->SpaceLocation[l].nAllocated)
      {
        char str[APECSS_STRINGLENGTH_SPRINTF];
        sprintf(str, "Results array allocation invalid (Result number: %i, Allocated %i)", Bubble->Results->Emissions->SpaceLocation[l].n,
                Bubble->Results->Emissions->SpaceLocation[l].nAllocated);
        apecss_erroronscreen(-1, str);
      }

      APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
      APECSS_FLOAT r = Bubble->Results->Emissions->SpaceLocation[l].RadialLocation;

      // Find the relevant nodes and interpolate the quantities to be written out
      if (r > Bubble->R)
      {
        struct APECSS_EmissionNode *Current = Bubble->Emissions->LastNode;

        if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE)
        {
          if (r < Current->r)
          {
            do
            {
              Current = Current->backward;
            } while (r < Current->r);

            APECSS_FLOAT p = (Current->forward->p * (r - Current->r) + Current->p * (Current->forward->r - r)) / (Current->forward->r - Current->r);

            Bubble->Results->Emissions->SpaceLocation[l].p[Bubble->Results->Emissions->SpaceLocation[l].n] = p;
            Bubble->Results->Emissions->SpaceLocation[l].u[Bubble->Results->Emissions->SpaceLocation[l].n] =
                (Current->forward->u * (r - Current->r) + Current->u * (Current->forward->r - r)) / (Current->forward->r - Current->r);
            Bubble->Results->Emissions->SpaceLocation[l].c[Bubble->Results->Emissions->SpaceLocation[l].n] =
                Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);
          }
          else
          {
            Bubble->Results->Emissions->SpaceLocation[l].p[Bubble->Results->Emissions->SpaceLocation[l].n] = 0.0;
            Bubble->Results->Emissions->SpaceLocation[l].u[Bubble->Results->Emissions->SpaceLocation[l].n] = 0.0;
            Bubble->Results->Emissions->SpaceLocation[l].c[Bubble->Results->Emissions->SpaceLocation[l].n] = 0.0;
          }
        }
        else if (Bubble->Emissions->Type == APECSS_EMISSION_QUASIACOUSTIC || Bubble->Emissions->Type == APECSS_EMISSION_FINITE_SPEED_INCOMPRESSIBLE)
        {
          if (r < Current->r)
          {
            do
            {
              Current = Current->backward;
            } while (r < Current->r);

            Bubble->Results->Emissions->SpaceLocation[l].p[Bubble->Results->Emissions->SpaceLocation[l].n] =
                (Current->forward->p * (r - Current->r) + Current->p * (Current->forward->r - r)) / (Current->forward->r - Current->r);
            Bubble->Results->Emissions->SpaceLocation[l].u[Bubble->Results->Emissions->SpaceLocation[l].n] =
                (Current->forward->u * (r - Current->r) + Current->u * (Current->forward->r - r)) / (Current->forward->r - Current->r);
          }
          else
          {
            Bubble->Results->Emissions->SpaceLocation[l].p[Bubble->Results->Emissions->SpaceLocation[l].n] = 0.0;
            Bubble->Results->Emissions->SpaceLocation[l].u[Bubble->Results->Emissions->SpaceLocation[l].n] = 0.0;
          }
        }
        else if (Bubble->Emissions->Type == APECSS_EMISSION_INCOMPRESSIBLE)
        {
          APECSS_FLOAT dotU = Bubble->ode[0](Bubble->ODEsSol, Bubble->t, Bubble);
          Bubble->Results->Emissions->SpaceLocation[l].p[Bubble->Results->Emissions->SpaceLocation[l].n] =
              pinf + Bubble->Liquid->rhoref * ((2.0 * Bubble->R * APECSS_POW2(Bubble->U) + dotU * APECSS_POW2(Bubble->R)) / r -
                                               APECSS_POW4(Bubble->R) * APECSS_POW2(Bubble->U) / (2.0 * APECSS_POW4(r)));
          Bubble->Results->Emissions->SpaceLocation[l].u[Bubble->Results->Emissions->SpaceLocation[l].n] = APECSS_POW2(Bubble->R) * Bubble->U / APECSS_POW2(r);
        }
        else
        {
          apecss_erroronscreen(-1, "Unknown emission type specified in space location.");
        }
      }
      else
      {
        Bubble->Results->Emissions->SpaceLocation[l].p[Bubble->Results->Emissions->SpaceLocation[l].n] = 0.0;
        Bubble->Results->Emissions->SpaceLocation[l].u[Bubble->Results->Emissions->SpaceLocation[l].n] = 0.0;

        if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE)
          Bubble->Results->Emissions->SpaceLocation[l].c[Bubble->Results->Emissions->SpaceLocation[l].n] = 0.0;
      }

      Bubble->Results->Emissions->SpaceLocation[l].t[Bubble->Results->Emissions->SpaceLocation[l].n] = Bubble->t;
      Bubble->Results->Emissions->SpaceLocation[l].pInf[Bubble->Results->Emissions->SpaceLocation[l].n] = pinf;

      ++(Bubble->Results->Emissions->SpaceLocation[l].n);
    }
  }

  return (0);
}

int apecss_results_emissionsspace_write(struct APECSS_Bubble *Bubble, int write)
{
  if (Bubble->Results->Emissions != NULL)
  {
    if (write)
    {
      for (register int l = 0; l < Bubble->Results->Emissions->nSpaceLocations; l++)
      {
        if (Bubble->Results->Emissions->SpaceLocation[l].n)
        {
          char path[APECSS_STRINGLENGTH_SPRINTF_LONG];

#if defined(APECSS_PRECISION_LONGDOUBLE)
          sprintf(path, "%s/EmissionsSpace_%.3Le.txt", Bubble->Results->dir, Bubble->Results->Emissions->SpaceLocation[l].RadialLocation);
#else
          sprintf(path, "%s/EmissionsSpace_%.3e.txt", Bubble->Results->dir, Bubble->Results->Emissions->SpaceLocation[l].RadialLocation);
#endif

          FILE *file_ptr;
          if (write == APECSS_RESULTS_WRITE)  // Create new file and write all data in one go.
          {
            file_ptr = fopen(path, "w");
            fprintf(file_ptr, "# time p u");
            if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE) fprintf(file_ptr, " c");
            fprintf(file_ptr, " pInf \n");
          }
          else  // Append file.
          {
            if (!(file_ptr = fopen(path, "r")))  // File does not exist yet -> write header.
            {
              file_ptr = fopen(path, "a");
              fprintf(file_ptr, "# time p u");
              if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE) fprintf(file_ptr, " c");
              fprintf(file_ptr, " pInf \n");
            }
            else  // File exists already, append.
            {
              fclose(file_ptr);
              file_ptr = fopen(path, "a");
            }
          }

#if defined(APECSS_PRECISION_LONGDOUBLE)
          for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++)
          {
            fprintf(file_ptr, "%.*Le %.*Le %.*Le", Bubble->Results->digits, Bubble->Results->Emissions->SpaceLocation[l].t[i], Bubble->Results->digits,
                    Bubble->Results->Emissions->SpaceLocation[l].p[i], Bubble->Results->digits, Bubble->Results->Emissions->SpaceLocation[l].u[i]);

            if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE)
              fprintf(file_ptr, " %.*Le", Bubble->Results->digits, Bubble->Results->Emissions->SpaceLocation[l].c[i]);

            fprintf(file_ptr, " %.*Le \n", Bubble->Results->digits, Bubble->Results->Emissions->SpaceLocation[l].pInf[i]);
          }
#else
          for (register int i = 0; i < Bubble->Results->Emissions->SpaceLocation[l].n; i++)
          {
            fprintf(file_ptr, "%.*e %.*e %.*e", Bubble->Results->digits, Bubble->Results->Emissions->SpaceLocation[l].t[i], Bubble->Results->digits,
                    Bubble->Results->Emissions->SpaceLocation[l].p[i], Bubble->Results->digits, Bubble->Results->Emissions->SpaceLocation[l].u[i]);

            if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE)
              fprintf(file_ptr, " %.*e", Bubble->Results->digits, Bubble->Results->Emissions->SpaceLocation[l].c[i]);

            fprintf(file_ptr, " %.*e \n", Bubble->Results->digits, Bubble->Results->Emissions->SpaceLocation[l].pInf[i]);
          }
#endif

          fclose(file_ptr);
        }
      }
    }

    apecss_results_emissionsspace_free(Bubble);
  }

  return (0);
}

int apecss_results_emissionsspace_free(struct APECSS_Bubble *Bubble)
{
  for (register int l = 0; l < Bubble->Results->Emissions->nSpaceLocations; l++)
  {
    if (Bubble->Results->Emissions->SpaceLocation[l].nAllocated)
    {
      free(Bubble->Results->Emissions->SpaceLocation[l].t);
      Bubble->Results->Emissions->SpaceLocation[l].t = NULL;
      free(Bubble->Results->Emissions->SpaceLocation[l].p);
      Bubble->Results->Emissions->SpaceLocation[l].p = NULL;
      free(Bubble->Results->Emissions->SpaceLocation[l].u);
      Bubble->Results->Emissions->SpaceLocation[l].u = NULL;

      if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE)
      {
        free(Bubble->Results->Emissions->SpaceLocation[l].c);
        Bubble->Results->Emissions->SpaceLocation[l].c = NULL;
      }

      free(Bubble->Results->Emissions->SpaceLocation[l].pInf);
      Bubble->Results->Emissions->SpaceLocation[l].pInf = NULL;

      Bubble->Results->Emissions->SpaceLocation[l].n = 0;
      Bubble->Results->Emissions->SpaceLocation[l].nAllocated = 0;
    }
  }

  return (0);
}

// -------------------------------------------------------------------
// EMISSIONS  (specific emission id)
// -------------------------------------------------------------------
// Functions handling the results associated with the acoustics
// emissions of user-defined nodes.
// -------------------------------------------------------------------

int apecss_results_emissionsnodespecific_storenone(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble)
{
  return (0);
}

int apecss_results_emissionsnodespecific_storeall(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble)
{
  for (register int i = 0; i < Bubble->Results->Emissions->nNodes; i++)  // Loop over all nodes that are to stored
  {
    // The node must fulfill the following criteria:
    // 1) It must be the actual node that was specified by the user or it was emitted afterwards, in case the specified node was discarded
    // 2) It must not be at the end of the linked list, as condition could not be evaluated
    // 3) The forward neighbor needs to had been emitted befor the specified node
    if (Node->id >= Bubble->Results->Emissions->Node[i].id && Node->forward != NULL && Node->forward->id < Bubble->Results->Emissions->Node[i].id)
    {
      int n = Bubble->Results->Emissions->Node[i].n;
      Bubble->Results->Emissions->Node[i].real_id[n] = Node->id;
      Bubble->Results->Emissions->Node[i].t[n] = Bubble->t;
      Bubble->Results->Emissions->Node[i].r[n] = Node->r;
      Bubble->Results->Emissions->Node[i].p[n] = Node->p;
      Bubble->Results->Emissions->Node[i].u[n] = Node->u;
      Bubble->Results->Emissions->Node[i].c[n] = c;
      Bubble->Results->Emissions->Node[i].pInf[n] = pinf;

      Bubble->Results->Emissions->Node[i].n++;
    }
  }

  return (0);
}

int apecss_results_emissionsnodespecific_allocnone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_results_emissionsnodespecific_allocall(struct APECSS_Bubble *Bubble)
{
  for (register int i = 0; i < Bubble->Results->Emissions->nNodes; i++)
  {
    int n = Bubble->Results->Emissions->Node[i].n;

    if (n == Bubble->Results->Emissions->Node[i].nAllocated)  // Extend the length of the result arrays
    {
      apecss_results_emissionsnode_allocnode(&(*Bubble).Results->Emissions->Node[i]);
    }
    else if (Bubble->Results->Emissions->Node[i].n > Bubble->Results->Emissions->Node[i].nAllocated)
    {
      char str[APECSS_STRINGLENGTH_SPRINTF];
      sprintf(str, "Array allocation for emission node %i is invalid (Result number: %i, Allocated %i)", Bubble->Results->Emissions->Node[i].id, n,
              Bubble->Results->Emissions->Node[i].nAllocated);
      apecss_erroronscreen(-1, str);
    }
  }

  return (0);
}

int apecss_results_emissionsnodespecific_write(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Results->Emissions != NULL)
  {
    for (register int i = 0; i < Bubble->Results->Emissions->nNodes; i++)
    {
      if (Bubble->Results->Emissions->Node[i].n)
      {
        char path[APECSS_STRINGLENGTH_SPRINTF_LONG];
        sprintf(path, "%s/EmissionsNode_%i.txt", Bubble->Results->dir, Bubble->Results->Emissions->Node[i].id);

        apecss_results_emissionsnode_writenode(&(*Bubble).Results->Emissions->Node[i], path, Bubble->Results->digits);
      }
    }

    apecss_results_emissionsnodespecific_free(Bubble);
  }

  return (0);
}

int apecss_results_emissionsnodespecific_free(struct APECSS_Bubble *Bubble)
{
  for (register int i = 0; i < Bubble->Results->Emissions->nNodes; i++) apecss_results_emissionsnode_freenode(&(*Bubble).Results->Emissions->Node[i]);

  return (0);
}

// -------------------------------------------------------------------
// EMISSIONS  (min/max values)
// -------------------------------------------------------------------
// Functions handling the results associated with the extrema of the
// acoustics emissions recorded in a user-defined oscillation period.
// -------------------------------------------------------------------

int apecss_results_emissionsnodeminmax_identifynone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_results_emissionsnodeminmax_identifyall(struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT excitation_period = Bubble->tEnd;
  if (Bubble->Excitation != NULL) excitation_period = 1.0 / Bubble->Excitation->f;

  APECSS_FLOAT t = Bubble->t + Bubble->dt;
  APECSS_FLOAT period = (APECSS_FLOAT) Bubble->Results->Emissions->MinMaxPeriod;

  if ((t >= (period - 1.0) * excitation_period) && (t < period * excitation_period))  // Check if this is the specified period
  {
    if (Bubble->R < Bubble->Results->Emissions->Rmin)
    {
      Bubble->Results->Emissions->Rmin = Bubble->R;
      apecss_results_emissionsnode_freenode(Bubble->Results->Emissions->Node_Rmin);
      apecss_results_emissionsnode_allocnode(Bubble->Results->Emissions->Node_Rmin);
      Bubble->Results->Emissions->Node_Rmin->id = Bubble->dtNumber;
    }

    if (Bubble->U < Bubble->Results->Emissions->Umin)
    {
      Bubble->Results->Emissions->Umin = Bubble->U;
      apecss_results_emissionsnode_freenode(Bubble->Results->Emissions->Node_Umin);
      apecss_results_emissionsnode_allocnode(Bubble->Results->Emissions->Node_Umin);
      Bubble->Results->Emissions->Node_Umin->id = Bubble->dtNumber;
    }

    APECSS_FLOAT pLnew = Bubble->Liquid->get_pressure_bubblewall(Bubble->ODEsSol, Bubble->t, Bubble);

    if (pLnew > Bubble->Results->Emissions->pLmax)
    {
      Bubble->Results->Emissions->pLmax = pLnew;
      apecss_results_emissionsnode_freenode(Bubble->Results->Emissions->Node_pLmax);
      apecss_results_emissionsnode_allocnode(Bubble->Results->Emissions->Node_pLmax);
      Bubble->Results->Emissions->Node_pLmax->id = Bubble->dtNumber;
    }
  }

  return (0);
}

int apecss_results_emissionsnodeminmax_allocnone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_results_emissionsnodeminmax_allocall(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Results->Emissions->Node_Rmin->n == Bubble->Results->Emissions->Node_Rmin->nAllocated)  // Extend the length of the result arrays
    apecss_results_emissionsnode_allocnode(Bubble->Results->Emissions->Node_Rmin);

  if (Bubble->Results->Emissions->Node_Umin->n == Bubble->Results->Emissions->Node_Umin->nAllocated)  // Extend the length of the result arrays
    apecss_results_emissionsnode_allocnode(Bubble->Results->Emissions->Node_Umin);

  if (Bubble->Results->Emissions->Node_pLmax->n == Bubble->Results->Emissions->Node_pLmax->nAllocated)  // Extend the length of the result arrays
    apecss_results_emissionsnode_allocnode(Bubble->Results->Emissions->Node_pLmax);

  return (0);
}

int apecss_results_emissionsnodeminmax_storenone(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble)
{
  return (0);
}

int apecss_results_emissionsnodeminmax_storeall(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble)
{
  // The node must fulfill the following criteria:
  // 1) It must be the actual node that was specified by the user or it was emitted afterwards, in case the specified node was discarded
  // 2) It must not be at the end of the linked list, as condition could not be evaluated
  // 3) The forward neighbor needs to had been emitted befor the specified node
  if (Node->id >= Bubble->Results->Emissions->Node_Rmin->id && Node->forward != NULL && Node->forward->id < Bubble->Results->Emissions->Node_Rmin->id)
  {
    int n = Bubble->Results->Emissions->Node_Rmin->n;
    Bubble->Results->Emissions->Node_Rmin->real_id[n] = Node->id;
    Bubble->Results->Emissions->Node_Rmin->t[n] = Bubble->t;
    Bubble->Results->Emissions->Node_Rmin->r[n] = Node->r;
    Bubble->Results->Emissions->Node_Rmin->p[n] = Node->p;
    Bubble->Results->Emissions->Node_Rmin->u[n] = Node->u;
    Bubble->Results->Emissions->Node_Rmin->c[n] = c;
    Bubble->Results->Emissions->Node_Rmin->pInf[n] = pinf;

    Bubble->Results->Emissions->Node_Rmin->n++;
  }

  if (Node->id >= Bubble->Results->Emissions->Node_Umin->id && Node->forward != NULL && Node->forward->id < Bubble->Results->Emissions->Node_Umin->id)
  {
    int n = Bubble->Results->Emissions->Node_Umin->n;
    Bubble->Results->Emissions->Node_Umin->real_id[n] = Node->id;
    Bubble->Results->Emissions->Node_Umin->t[n] = Bubble->t;
    Bubble->Results->Emissions->Node_Umin->r[n] = Node->r;
    Bubble->Results->Emissions->Node_Umin->p[n] = Node->p;
    Bubble->Results->Emissions->Node_Umin->u[n] = Node->u;
    Bubble->Results->Emissions->Node_Umin->c[n] = c;
    Bubble->Results->Emissions->Node_Umin->pInf[n] = pinf;

    Bubble->Results->Emissions->Node_Umin->n++;
  }

  if (Node->id >= Bubble->Results->Emissions->Node_pLmax->id && Node->forward != NULL && Node->forward->id < Bubble->Results->Emissions->Node_pLmax->id)
  {
    int n = Bubble->Results->Emissions->Node_pLmax->n;
    Bubble->Results->Emissions->Node_pLmax->real_id[n] = Node->id;
    Bubble->Results->Emissions->Node_pLmax->t[n] = Bubble->t;
    Bubble->Results->Emissions->Node_pLmax->r[n] = Node->r;
    Bubble->Results->Emissions->Node_pLmax->p[n] = Node->p;
    Bubble->Results->Emissions->Node_pLmax->u[n] = Node->u;
    Bubble->Results->Emissions->Node_pLmax->c[n] = c;
    Bubble->Results->Emissions->Node_pLmax->pInf[n] = pinf;

    Bubble->Results->Emissions->Node_pLmax->n++;
  }

  return (0);
}

int apecss_results_emissionsnodeminmax_write(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Results->Emissions != NULL && Bubble->Results->Emissions->MinMaxPeriod)
  {
    if (Bubble->Results->Emissions->Node_Rmin->n)
    {
      char path[APECSS_STRINGLENGTH_SPRINTF_LONG];
      sprintf(path, "%s/EmissionsNode_Rmin.txt", Bubble->Results->dir);

      apecss_results_emissionsnode_writenode(Bubble->Results->Emissions->Node_Rmin, path, Bubble->Results->digits);
      apecss_results_emissionsnode_freenode(Bubble->Results->Emissions->Node_Rmin);
    }

    if (Bubble->Results->Emissions->Node_Umin->n)
    {
      char path[APECSS_STRINGLENGTH_SPRINTF_LONG];
      sprintf(path, "%s/EmissionsNode_Umin.txt", Bubble->Results->dir);

      apecss_results_emissionsnode_writenode(Bubble->Results->Emissions->Node_Umin, path, Bubble->Results->digits);
      apecss_results_emissionsnode_freenode(Bubble->Results->Emissions->Node_Umin);
    }

    if (Bubble->Results->Emissions->Node_pLmax->n)
    {
      char path[APECSS_STRINGLENGTH_SPRINTF_LONG];
      sprintf(path, "%s/EmissionsNode_pLmax.txt", Bubble->Results->dir);

      apecss_results_emissionsnode_writenode(Bubble->Results->Emissions->Node_pLmax, path, Bubble->Results->digits);
      apecss_results_emissionsnode_freenode(Bubble->Results->Emissions->Node_pLmax);
    }
  }

  return (0);
}

// -------------------------------------------------------------------
// EMISSIONS  (general ops for nodes)
// -------------------------------------------------------------------
// Functions with general operation associated with the results of the
// acoustic emissions.
// -------------------------------------------------------------------

int apecss_results_emissionsnode_initializenode(struct APECSS_ResultsEmissionsNode *Node)
{
  Node->n = 0;
  Node->nAllocated = 0;
  Node->real_id = NULL;
  Node->r = NULL;
  Node->t = NULL;
  Node->p = NULL;
  Node->u = NULL;
  Node->c = NULL;
  Node->pInf = NULL;

  return (0);
}

int apecss_results_emissionsnode_allocnone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_results_emissionsnode_allocall(struct APECSS_Bubble *Bubble)
{
  Bubble->Results->Emissions->allocation_nodespecific(Bubble);
  Bubble->Results->Emissions->allocation_nodeminmax(Bubble);

  return (0);
}

int apecss_results_emissionsnode_allocnode(struct APECSS_ResultsEmissionsNode *Node)
{
  Node->nAllocated = Node->n + APECSS_DATA_ALLOC_INCREMENT;

  if (Node->real_id != NULL)
  {
    int *tempint;
    tempint = malloc(Node->n * sizeof(int));

    for (register int i = 0; i < Node->n; i++) tempint[i] = Node->real_id[i];
    free(Node->real_id);
    Node->real_id = malloc(Node->nAllocated * sizeof(int));
    for (register int i = 0; i < Node->n; i++) Node->real_id[i] = tempint[i];

    free(tempint);

    APECSS_FLOAT *temp;
    temp = malloc(Node->n * sizeof(APECSS_FLOAT));

    for (register int i = 0; i < Node->n; i++) temp[i] = Node->r[i];
    free(Node->r);
    Node->r = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    for (register int i = 0; i < Node->n; i++) Node->r[i] = temp[i];

    for (register int i = 0; i < Node->n; i++) temp[i] = Node->t[i];
    free(Node->t);
    Node->t = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    for (register int i = 0; i < Node->n; i++) Node->t[i] = temp[i];

    for (register int i = 0; i < Node->n; i++) temp[i] = Node->p[i];
    free(Node->p);
    Node->p = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    for (register int i = 0; i < Node->n; i++) Node->p[i] = temp[i];

    for (register int i = 0; i < Node->n; i++) temp[i] = Node->u[i];
    free(Node->u);
    Node->u = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    for (register int i = 0; i < Node->n; i++) Node->u[i] = temp[i];

    for (register int i = 0; i < Node->n; i++) temp[i] = Node->c[i];
    free(Node->c);
    Node->c = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    for (register int i = 0; i < Node->n; i++) Node->c[i] = temp[i];

    for (register int i = 0; i < Node->n; i++) temp[i] = Node->pInf[i];
    free(Node->pInf);
    Node->pInf = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    for (register int i = 0; i < Node->n; i++) Node->pInf[i] = temp[i];

    free(temp);
  }
  else
  {
    Node->real_id = malloc(Node->nAllocated * sizeof(int));
    Node->r = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    Node->t = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    Node->p = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    Node->u = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    Node->c = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
    Node->pInf = malloc(Node->nAllocated * sizeof(APECSS_FLOAT));
  }

  return (0);
}

int apecss_results_emissionsnode_storenone(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble) { return (0); }

int apecss_results_emissionsnode_storeall(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble)
{
  Bubble->Results->Emissions->store_nodespecific(Node, c, pinf, Bubble);
  Bubble->Results->Emissions->store_nodeminmax(Node, c, pinf, Bubble);

  return (0);
}

int apecss_results_emissionsnode_writenode(struct APECSS_ResultsEmissionsNode *Node, char *path, int digits)
{
  FILE *file_ptr;
  file_ptr = fopen(path, "w");

  fprintf(file_ptr, "# real-id r t p u c pinf \n");

  for (register int i = 0; i < Node->n; i++)
  {
    fprintf(file_ptr, "%i", Node->real_id[i]);
    fprintf(file_ptr, " ");

#if defined(APECSS_PRECISION_LONGDOUBLE)
    fprintf(file_ptr, "%.*Le %.*Le %.*Le %.*Le %.*Le %.*Le \n", digits, Node->r[i], digits, Node->t[i], digits, Node->p[i], digits, Node->u[i], digits,
            Node->c[i], digits, Node->pInf[i]);
#else
    fprintf(file_ptr, "%.*e %.*e %.*e %.*e %.*e %.*e \n", digits, Node->r[i], digits, Node->t[i], digits, Node->p[i], digits, Node->u[i], digits, Node->c[i],
            digits, Node->pInf[i]);
#endif
  }

  fclose(file_ptr);

  return (0);
}

int apecss_results_emissionsnode_freenode(struct APECSS_ResultsEmissionsNode *Node)
{
  if (Node->nAllocated)
  {
    free(Node->real_id);
    Node->real_id = NULL;
    free(Node->r);
    Node->r = NULL;
    free(Node->t);
    Node->t = NULL;
    free(Node->p);
    Node->p = NULL;
    free(Node->u);
    Node->u = NULL;
    free(Node->c);
    Node->c = NULL;
    free(Node->pInf);
    Node->pInf = NULL;

    Node->n = 0;
    Node->nAllocated = 0;
  }

  return (0);
}