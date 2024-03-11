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
// OPTIONS
// -------------------------------------------------------------------
// Functions initializing, processing and handling the options and the
// related data structures of a bubble.
// -------------------------------------------------------------------

int apecss_bubble_initializestruct(struct APECSS_Bubble *Bubble)
{
  // Initialize pointers to arrays
  Bubble->ODEsSol = NULL;
  Bubble->ODEsSolOld = NULL;
  Bubble->k2 = NULL;
  Bubble->k3 = NULL;
  Bubble->k4 = NULL;
  Bubble->k5 = NULL;
  Bubble->k6 = NULL;
  Bubble->k7 = NULL;
  Bubble->kLast = NULL;

  // Initialize pointer to structures
  Bubble->NumericsODE = NULL;
  Bubble->Gas = NULL;
  Bubble->Liquid = NULL;
  Bubble->Interface = NULL;
  Bubble->Excitation = NULL;
  Bubble->Emissions = NULL;
  Bubble->Results = NULL;

  // Initialize pointers to functions
  Bubble->ode = NULL;
  Bubble->get_pressure_infinity = NULL;
  Bubble->get_pressurederivative_infinity = NULL;
  Bubble->emissions_initialize = NULL;
  Bubble->emissions_update = NULL;
  Bubble->emissions_free = NULL;
  Bubble->progress_initial = NULL;
  Bubble->progress_update = NULL;
  Bubble->progress_final = NULL;
  Bubble->results_rayleighplesset_store = NULL;
  Bubble->results_emissionstime_write = NULL;
  Bubble->results_emissionstime_check = NULL;
  Bubble->results_emissionsspace_store = NULL;
  Bubble->results_emissionsnode_alloc = NULL;
  Bubble->results_emissionsnode_store = NULL;
  Bubble->results_emissionsnodeminmax_identify = NULL;

  return (0);
}

int apecss_bubble_setdefaultoptions(struct APECSS_Bubble *Bubble)
{
  // Time
  Bubble->tStart = 0.0;
  Bubble->tEnd = 1.0;
  Bubble->dt = 1.0e-10;

  // Governing RP model
  Bubble->RPModel = APECSS_BUBBLEMODEL_RP;
  Bubble->dimensionality = 2.0;

  // Bubble properties
  Bubble->R0 = 1.0;
  Bubble->R = Bubble->R0;
  Bubble->U0 = 0.0;
  Bubble->U = Bubble->U0;
  Bubble->p0 = 1.0e5;
  Bubble->pG0 = -APECSS_LARGE;  // Set by the user or when the bubble is initialized
  Bubble->rhoG0 = 1.0;
  Bubble->r_hc = 0.0;

  // ODEs
  Bubble->nODEs = 2;  // ODEs for the bubble radius and bubble wall velocity
  Bubble->nUserODEs = 0;  // Number of ODEs defined by the user
  Bubble->dtNumber = 0;  // Time-step count
  Bubble->nSubIter = 25;  // Number of maximum sub-iterations in a given time-step

  // Interface
  Bubble->Rbuck = 1.0e-6;
  Bubble->Rrupt = 1.0e-6;
  Bubble->GompertzB = 0.0;

  // Set default pointers to functions
  Bubble->get_pressure_infinity = apecss_bubble_pressure_infinity_noexcitation;
  Bubble->get_pressurederivative_infinity = apecss_bubble_pressurederivative_infinity_noexcitation;
  Bubble->get_dimensionalradius = apecss_bubble_dimensionalradius_spherical;
  Bubble->emissions_initialize = apecss_emissions_initializenone;
  Bubble->emissions_update = apecss_emissions_updatenone;
  Bubble->emissions_free = apecss_emissions_freenone;
  Bubble->progress_initial = apecss_bubble_solver_progress_initialnone;
  Bubble->progress_update = apecss_bubble_solver_progress_updatenone;
  Bubble->progress_final = apecss_bubble_solver_progress_finalnone;
  Bubble->results_rayleighplesset_store = apecss_results_rayleighplesset_storenone;
  Bubble->results_emissionstime_write = apecss_results_emissionstime_writenone;
  Bubble->results_emissionstime_check = apecss_results_emissionstime_checknone;
  Bubble->results_emissionsspace_store = apecss_results_emissionsspace_storenone;
  Bubble->results_emissionsnode_alloc = apecss_results_emissionsnode_allocnone;
  Bubble->results_emissionsnode_store = apecss_results_emissionsnode_storenone;
  Bubble->results_emissionsnodeminmax_identify = apecss_results_emissionsnodeminmax_identifynone;

  return (0);
}

int apecss_bubble_readoptions(struct APECSS_Bubble *Bubble, char *OptionsDir)
{
  int l = 0;
  int line = 0;
  FILE *OptionsFile;
  char str[APECSS_STRINGLENGTH_SPRINTF];
  char option[APECSS_STRINGLENGTH], option2[APECSS_STRINGLENGTH], option3[APECSS_STRINGLENGTH];
  int StatusFile = 1;
  int StatusSection = 1;

  if ((OptionsFile = fopen(OptionsDir, "r")) == (FILE *) NULL)
  {
    sprintf(str, "File %s cannot be opened for reading.\n", OptionsDir);
    apecss_erroronscreen(1, str);
  }

  while ((l = apecss_readoneoption(OptionsFile, option)) != EOF && StatusFile == 1)
  {
    line += l;

    if (strncasecmp(option, "bubble", 6) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else if (strncasecmp(option2, "rpmodel", 7) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "gilmore", 7) == 0)
            {
              Bubble->RPModel = APECSS_BUBBLEMODEL_GILMORE;
            }
            else if (strncasecmp(option3, "km", 2) == 0)
            {
              Bubble->RPModel = APECSS_BUBBLEMODEL_KELLERMIKSIS;
            }
            else if (strncasecmp(option3, "rpar", 4) == 0)
            {
              Bubble->RPModel = APECSS_BUBBLEMODEL_RP_ACOUSTICRADIATION;
            }
            else if (strncasecmp(option3, "rp", 2) == 0)
            {
              Bubble->RPModel = APECSS_BUBBLEMODEL_RP;
            }
          }
        }
        else if (strncasecmp(option2, "emissions", 9) == 0)
        {
          if (Bubble->Emissions == NULL) apecss_emissions_initializestruct(Bubble);

          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;

            if (strncasecmp(option3, "ic", 2) == 0 || strncasecmp(option3, "incompressible", 14) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_INCOMPRESSIBLE;
            }
            else if (strncasecmp(option3, "fsic", 3) == 0 || strncasecmp(option3, "finitespeedincompressible", 25) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_FINITE_SPEED_INCOMPRESSIBLE;
            }
            else if (strncasecmp(option3, "qa", 2) == 0 || strncasecmp(option3, "quasiacoustic", 13) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_QUASIACOUSTIC;
            }
            else if (strncasecmp(option3, "ev", 2) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_EV;
            }
            else if (strncasecmp(option3, "tiv", 3) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_TIV;
            }

            l = apecss_readoneoption(OptionsFile, option3);
            Bubble->Emissions->CutOffDistance = APECSS_STRINGTOFLOAT(option3);
          }
        }
        else if (strncasecmp(option2, "emissionintegration", 19) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;

            if (strncasecmp(option3, "euler", 5) == 0)
            {
              Bubble->Emissions->Scheme = APECSS_EMISSION_INTEGRATE_EULER;
            }
            else if (strncasecmp(option3, "rk4", 3) == 0)
            {
              Bubble->Emissions->Scheme = APECSS_EMISSION_INTEGRATE_RK4;
            }
          }
        }
        else if (strncasecmp(option2, "kbitertolerance", 15) == 0)
        {
          if (Bubble->Emissions == NULL) apecss_emissions_initializestruct(Bubble);

          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Emissions->KB_IterTolerance = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "dimensionality", 14) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "plane", 6) == 0)
            {
              Bubble->dimensionality = 0.0;
              Bubble->get_dimensionalradius = apecss_bubble_dimensionalradius_planar;
            }
            else if (strncasecmp(option3, "cylinder", 8) == 0)
            {
              Bubble->dimensionality = 1.0;
              Bubble->get_dimensionalradius = apecss_bubble_dimensionalradius_cylindrical;
            }
            else
            {
              Bubble->dimensionality = 2.0;
              Bubble->get_dimensionalradius = apecss_bubble_dimensionalradius_spherical;
            }
          }
        }
        else if (strncasecmp(option2, "initialradius", 13) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->R0 = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "pressureambient", 15) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->p0 = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "initialgaspressure", 18) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->pG0 = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "hardcoreradius", 14) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->r_hc = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of BUBBLE is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    else if (strncasecmp(option, "results", 7) == 0)
    {
      if (Bubble->Results == NULL) apecss_results_initializestruct(Bubble);

      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else if (strncasecmp(option2, "bubble", 6) == 0)
        {
          if (Bubble->Results->RayleighPlesset == NULL) apecss_results_rayleighplesset_initializestruct(Bubble);
        }
        else if (strncasecmp(option2, "outputfreqrp", 12) == 0)
        {
          if (Bubble->Results->RayleighPlesset == NULL) apecss_results_rayleighplesset_initializestruct(Bubble);

          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Results->RayleighPlesset->freq = atoi(option3);
        }
        else if (strncasecmp(option2, "outputdigits", 12) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Results->digits = atoi(option3);
        }
        else if (strncasecmp(option2, "outputpath", 10) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            sprintf(Bubble->Results->dir, "%s", option3);

            struct stat st = {0};
            if (stat(Bubble->Results->dir, &st) == -1) mkdir(Bubble->Results->dir, 0700);
          }
        }
        else if (strncasecmp(option2, "outputfreqEmissionsSpace", 24) == 0)
        {
          if (Bubble->Results->Emissions == NULL) apecss_results_emissions_initializestruct(Bubble);

          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Results->Emissions->freqSpaceLocations = atoi(option3);
        }
        else if (strncasecmp(option2, "emissionstime", 13) == 0)
        {
          if (Bubble->Results->Emissions == NULL) apecss_results_emissions_initializestruct(Bubble);

          l = apecss_readoneoption(OptionsFile, option3);
          APECSS_FLOAT t = APECSS_STRINGTOFLOAT(option3);

          if (Bubble->Results->Emissions->nTimeInstances)
          {
            APECSS_FLOAT *temp = malloc(Bubble->Results->Emissions->nTimeInstances * sizeof(APECSS_FLOAT));
            for (register int l = 0; l < Bubble->Results->Emissions->nTimeInstances; l++) temp[l] = Bubble->Results->Emissions->TimeInstances[l];
            free(Bubble->Results->Emissions->TimeInstances);

            Bubble->Results->Emissions->nTimeInstances++;
            Bubble->Results->Emissions->TimeInstances = malloc(Bubble->Results->Emissions->nTimeInstances * sizeof(APECSS_FLOAT));
            for (register int l = 0; l < Bubble->Results->Emissions->nTimeInstances - 1; l++) Bubble->Results->Emissions->TimeInstances[l] = temp[l];
            free(temp);
          }
          else
          {
            Bubble->Results->Emissions->nTimeInstances++;
            Bubble->Results->Emissions->TimeInstances = malloc(Bubble->Results->Emissions->nTimeInstances * sizeof(APECSS_FLOAT));
          }

          Bubble->Results->Emissions->TimeInstances[Bubble->Results->Emissions->nTimeInstances - 1] = t;
        }
        else if (strncasecmp(option2, "emissionsspace", 14) == 0)
        {
          if (Bubble->Results->Emissions == NULL) apecss_results_emissions_initializestruct(Bubble);

          l = apecss_readoneoption(OptionsFile, option3);
          APECSS_FLOAT r = APECSS_STRINGTOFLOAT(option3);

          if (Bubble->Results->Emissions->nSpaceLocations)
          {
            struct APECSS_ResultsEmissionsSpace *temp;
            temp = (struct APECSS_ResultsEmissionsSpace *) malloc(Bubble->Results->Emissions->nSpaceLocations * sizeof(struct APECSS_ResultsEmissionsSpace));
            for (register int l = 0; l < Bubble->Results->Emissions->nSpaceLocations; l++) temp[l] = Bubble->Results->Emissions->SpaceLocation[l];
            free(Bubble->Results->Emissions->SpaceLocation);

            Bubble->Results->Emissions->nSpaceLocations++;
            Bubble->Results->Emissions->SpaceLocation =
                (struct APECSS_ResultsEmissionsSpace *) malloc(Bubble->Results->Emissions->nSpaceLocations * sizeof(struct APECSS_ResultsEmissionsSpace));
            for (register int l = 0; l < Bubble->Results->Emissions->nSpaceLocations - 1; l++) Bubble->Results->Emissions->SpaceLocation[l] = temp[l];
            free(temp);
          }
          else
          {
            Bubble->Results->Emissions->nSpaceLocations++;
            Bubble->Results->Emissions->SpaceLocation =
                (struct APECSS_ResultsEmissionsSpace *) malloc(Bubble->Results->Emissions->nSpaceLocations * sizeof(struct APECSS_ResultsEmissionsSpace));
          }

          Bubble->Results->Emissions->SpaceLocation[Bubble->Results->Emissions->nSpaceLocations - 1].RadialLocation = r;
          Bubble->Results->Emissions->SpaceLocation[Bubble->Results->Emissions->nSpaceLocations - 1].n = 0;
          Bubble->Results->Emissions->SpaceLocation[Bubble->Results->Emissions->nSpaceLocations - 1].nAllocated = 0;
          Bubble->Results->Emissions->SpaceLocation[Bubble->Results->Emissions->nSpaceLocations - 1].t = NULL;
          Bubble->Results->Emissions->SpaceLocation[Bubble->Results->Emissions->nSpaceLocations - 1].p = NULL;
          Bubble->Results->Emissions->SpaceLocation[Bubble->Results->Emissions->nSpaceLocations - 1].u = NULL;
          Bubble->Results->Emissions->SpaceLocation[Bubble->Results->Emissions->nSpaceLocations - 1].c = NULL;
          Bubble->Results->Emissions->SpaceLocation[Bubble->Results->Emissions->nSpaceLocations - 1].pInf = NULL;
        }
        else if (strncasecmp(option2, "emissionsnode", 13) == 0)
        {
          if (Bubble->Results->Emissions == NULL) apecss_results_emissions_initializestruct(Bubble);

          l = apecss_readoneoption(OptionsFile, option3);
          int id = atoi(option3);

          if (Bubble->Results->Emissions->nNodes)
          {
            struct APECSS_ResultsEmissionsNode *temp;
            temp = (struct APECSS_ResultsEmissionsNode *) malloc(Bubble->Results->Emissions->nNodes * sizeof(struct APECSS_ResultsEmissionsNode));
            for (register int l = 0; l < Bubble->Results->Emissions->nNodes; l++) temp[l] = Bubble->Results->Emissions->Node[l];
            free(Bubble->Results->Emissions->Node);

            Bubble->Results->Emissions->nNodes++;
            Bubble->Results->Emissions->Node =
                (struct APECSS_ResultsEmissionsNode *) malloc(Bubble->Results->Emissions->nNodes * sizeof(struct APECSS_ResultsEmissionsNode));
            for (register int l = 0; l < Bubble->Results->Emissions->nNodes - 1; l++) Bubble->Results->Emissions->Node[l] = temp[l];
            free(temp);
          }
          else
          {
            Bubble->Results->Emissions->nNodes++;
            Bubble->Results->Emissions->Node =
                (struct APECSS_ResultsEmissionsNode *) malloc(Bubble->Results->Emissions->nNodes * sizeof(struct APECSS_ResultsEmissionsNode));
          }

          Bubble->Results->Emissions->Node[Bubble->Results->Emissions->nNodes - 1].id = id;
          apecss_results_emissionsnode_initializenode(&(*Bubble).Results->Emissions->Node[Bubble->Results->Emissions->nNodes - 1]);
        }
        else if (strncasecmp(option2, "emissionsminmax", 15) == 0)
        {
          if (Bubble->Results->Emissions == NULL) apecss_results_emissions_initializestruct(Bubble);

          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Results->Emissions->MinMaxPeriod = atoi(option3);

          Bubble->Results->Emissions->Rmin = APECSS_LARGE;
          Bubble->Results->Emissions->Node_Rmin = (struct APECSS_ResultsEmissionsNode *) malloc(sizeof(struct APECSS_ResultsEmissionsNode));
          apecss_results_emissionsnode_initializenode(Bubble->Results->Emissions->Node_Rmin);

          Bubble->Results->Emissions->Umin = APECSS_LARGE;
          Bubble->Results->Emissions->Node_Umin = (struct APECSS_ResultsEmissionsNode *) malloc(sizeof(struct APECSS_ResultsEmissionsNode));
          apecss_results_emissionsnode_initializenode(Bubble->Results->Emissions->Node_Umin);

          Bubble->Results->Emissions->pLmax = -APECSS_LARGE;
          Bubble->Results->Emissions->Node_pLmax = (struct APECSS_ResultsEmissionsNode *) malloc(sizeof(struct APECSS_ResultsEmissionsNode));
          apecss_results_emissionsnode_initializenode(Bubble->Results->Emissions->Node_pLmax);
        }
        else
        {
          sprintf(str, "An unknown option of RESULTS is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    else if (strncasecmp(option, "gas", 3) == 0 || strncasecmp(option, "liquid", 6) == 0 || strncasecmp(option, "interface", 9) == 0 ||
             strncasecmp(option, "odesolver", 9) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else
        {
          // Nothing to be done here.
        }
      }
    }
    else
    {
      sprintf(str, "An unknown Section is given: %s, line %i", option, line);
      apecss_erroronscreen(1, str);
    }
  }

  fclose(OptionsFile);

  return (0);
}

int apecss_bubble_processoptions(struct APECSS_Bubble *Bubble)
{
  // Increment the number of ODEs if additional ODEs for viscoelasticity are to be solved
  if (Bubble->Liquid->Type == APECSS_LIQUID_ZENER || Bubble->Liquid->Type == APECSS_LIQUID_OLDROYDB) Bubble->nODEs += 2;

  // Allocate the solution array and the array for the function pointers of the ODEs
  Bubble->ode = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->ODEsSol = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->ODEsSolOld = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  for (register int i = 0; i < Bubble->nODEs; i++) Bubble->ODEsSol[i] = 0.0;

  // Alocate the arrays for the intermediate solutions of the ODE for the Runge-Kutta solver
  Bubble->k2 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k3 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k4 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k5 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k6 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->k7 = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));
  Bubble->kLast = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));

  // ---------------------------------------
  // RP model

  if (Bubble->RPModel == APECSS_BUBBLEMODEL_GILMORE)
    Bubble->ode[0] = apecss_rp_gilmorevelocity_ode;
  else if (Bubble->RPModel == APECSS_BUBBLEMODEL_KELLERMIKSIS)
    Bubble->ode[0] = apecss_rp_kellermiksisvelocity_ode;
  else if (Bubble->RPModel == APECSS_BUBBLEMODEL_RP_ACOUSTICRADIATION)
    Bubble->ode[0] = apecss_rp_rayleighplessetacousticrationvelocity_ode;
  else
    Bubble->ode[0] = apecss_rp_rayleighplessetvelocity_ode;

  Bubble->ode[1] = apecss_rp_bubbleradius_ode;

  // ---------------------------------------
  // Bubble->nODEs is now reset (!!!), as it has to reflect the number of actually defined ODEs,
  // e.g. so that additional user-defined ODEs can easily be added.

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
  // Emissions linked list

  if (Bubble->Emissions != NULL)
  {
    if (Bubble->Emissions->Type != APECSS_EMISSION_INCOMPRESSIBLE)
    {
      Bubble->emissions_initialize = apecss_emissions_initializelinkedlist;
      Bubble->emissions_update = apecss_emissions_updatelinkedlist;
      Bubble->emissions_free = apecss_emissions_freelinkedlist;

      if (Bubble->Emissions->Type == APECSS_EMISSION_FINITE_SPEED_INCOMPRESSIBLE)
      {
        Bubble->Emissions->advance = apecss_emissions_advance_finitespeedincompressible;
        Bubble->Emissions->compute_f = apecss_emissions_f_finitespeedincompressible;
      }
      else if (Bubble->Emissions->Type == APECSS_EMISSION_QUASIACOUSTIC)
      {
        Bubble->Emissions->advance = apecss_emissions_advance_quasiacoustic;
        Bubble->Emissions->compute_f = apecss_emissions_f_quasiacoustic;
      }
      else if (Bubble->Emissions->Type & APECSS_EMISSION_KIRKWOODBETHE)
      {
        if (Bubble->Liquid->EoS == APECSS_LIQUID_TAIT)
        {
          Bubble->Emissions->advance = apecss_emissions_advance_kirkwoodbethe_tait;

          if (Bubble->Emissions->Type == APECSS_EMISSION_EV)
          {
            Bubble->Emissions->compute_f = apecss_emissions_f_kirkwoodbethe;

            if (Bubble->Emissions->Scheme == APECSS_EMISSION_INTEGRATE_RK4)
              Bubble->Emissions->integrate_along_characteristic = apecss_emissions_integrate_ev_tait_rk4;
            else
              Bubble->Emissions->integrate_along_characteristic = apecss_emissions_integrate_ev_tait_euler;
          }
          else if (Bubble->Emissions->Type == APECSS_EMISSION_TIV)
          {
            Bubble->Emissions->compute_f = apecss_emissions_f_zero;

            if (Bubble->Emissions->Scheme == APECSS_EMISSION_INTEGRATE_RK4)
              Bubble->Emissions->integrate_along_characteristic = apecss_emissions_integrate_tiv_tait_rk4;
            else
              Bubble->Emissions->integrate_along_characteristic = apecss_emissions_integrate_tiv_tait_euler;
          }
        }
        else if (Bubble->Liquid->EoS == APECSS_LIQUID_NASG)
        {
          Bubble->Emissions->advance = apecss_emissions_advance_kirkwoodbethe_general;

          if (Bubble->Emissions->Type == APECSS_EMISSION_EV)
          {
            Bubble->Emissions->compute_f = apecss_emissions_f_kirkwoodbethe;

            if (Bubble->Emissions->Scheme == APECSS_EMISSION_INTEGRATE_RK4)
              Bubble->Emissions->integrate_along_characteristic = apecss_emissions_integrate_ev_general_rk4;
            else
              Bubble->Emissions->integrate_along_characteristic = apecss_emissions_integrate_ev_general_euler;
          }
          else if (Bubble->Emissions->Type == APECSS_EMISSION_TIV)
          {
            Bubble->Emissions->compute_f = apecss_emissions_f_zero;

            if (Bubble->Emissions->Scheme == APECSS_EMISSION_INTEGRATE_RK4)
              Bubble->Emissions->integrate_along_characteristic = apecss_emissions_integrate_tiv_general_rk4;
            else
              Bubble->Emissions->integrate_along_characteristic = apecss_emissions_integrate_tiv_general_euler;
          }
        }
        else
          apecss_erroronscreen(-1, "Unknown equation of state defined for the liquid for the emissions.");
      }
    }
  }
  else
  {
    Bubble->emissions_initialize = apecss_emissions_initializenone;
    Bubble->emissions_update = apecss_emissions_updatenone;
    Bubble->emissions_free = apecss_emissions_freenone;
  }

  // ---------------------------------------
  // Excitation

  if (Bubble->Excitation != NULL)
  {
    if (Bubble->Excitation->type == APECSS_EXCITATION_SIN)
    {
      Bubble->get_pressure_infinity = apecss_bubble_pressure_infinity_sinexcitation;
      Bubble->get_pressurederivative_infinity = apecss_bubble_pressurederivative_infinity_sinexcitation;
    }
    else
    {
      Bubble->get_pressure_infinity = apecss_bubble_pressure_infinity_noexcitation;
      Bubble->get_pressurederivative_infinity = apecss_bubble_pressurederivative_infinity_noexcitation;
    }
  }
  else
  {
    Bubble->get_pressure_infinity = apecss_bubble_pressure_infinity_noexcitation;
    Bubble->get_pressurederivative_infinity = apecss_bubble_pressurederivative_infinity_noexcitation;
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

      Bubble->results_emissionsnode_alloc = apecss_results_emissionsnode_allocnone;
      Bubble->results_emissionsnode_store = apecss_results_emissionsnode_storenone;

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

  return (0);
}

int apecss_bubble_initialize(struct APECSS_Bubble *Bubble)
{
  for (register int i = 0; i < Bubble->nODEs; i++) Bubble->ODEsSol[i] = 0.0;

  // ---------------------------------------
  // Hardcore radius

  if (Bubble->Gas->EoS == APECSS_GAS_HC && Bubble->Gas->mmol > 0.0 && Bubble->Gas->dmol > 0.0)
  {
    // Compute the hardcore radius based on the molecule properties
    APECSS_FLOAT mgref = Bubble->Gas->rhoref * 4.0 * APECSS_PI * APECSS_POW3(Bubble->R0) / 3.0;
    Bubble->r_hc = APECSS_POW(APECSS_AVOGADRO * mgref * 4.0 * APECSS_POW3(0.5 * Bubble->Gas->dmol) / (Bubble->Gas->mmol), APECSS_ONETHIRD);
  }

  // ---------------------------------------
  // Interface

  if (Bubble->Interface->LipidCoatingModel & APECSS_LIPIDCOATING_MARMOTTANT)
  {
    if (Bubble->pG0 < -Bubble->Gas->B) Bubble->pG0 = Bubble->p0 + 2.0 * Bubble->Interface->sigma0 / Bubble->R0;

    Bubble->Rbuck = Bubble->R0 / APECSS_SQRT(1.0 + Bubble->Interface->sigma0 / Bubble->Interface->Elasticity);
    Bubble->Rrupt = Bubble->Rbuck * APECSS_SQRT(1.0 + Bubble->Interface->sigma / Bubble->Interface->Elasticity);

    if (Bubble->Interface->LipidCoatingModel & APECSS_LIPIDCOATING_GOMPERTZFUNCTION)
    {
      Bubble->GompertzB =
          -APECSS_LOG(Bubble->Interface->sigma0 / Bubble->Interface->sigma) / APECSS_EXP(Bubble->Interface->GompertzC * (1.0 - Bubble->R0 / Bubble->Rbuck));
    }
  }
  else
  {
    // Set the initial gas pressure to the Laplace pressure
    if (Bubble->pG0 < -Bubble->Gas->B) Bubble->pG0 = Bubble->p0 + 2.0 * Bubble->Interface->sigma / Bubble->R0;
  }

  // ---------------------------------------
  // Initial gas pressure assumed to result from polytropic expansion/compression from the reference state

  if (Bubble->Gas->EoS == APECSS_GAS_NASG)
  {
    APECSS_FLOAT peff = APECSS_POW((Bubble->pG0 + Bubble->Gas->B), (1.0 / Bubble->Gas->Gamma));
    Bubble->rhoG0 = Bubble->Gas->Kref * peff / (1.0 + Bubble->Gas->b * Bubble->Gas->Kref * peff);
  }
  else
  {
    Bubble->rhoG0 = Bubble->Gas->rhoref * APECSS_POW((Bubble->pG0 / Bubble->Gas->pref), (1.0 / Bubble->Gas->Gamma));
  }

  // ---------------------------------------
  // Initialize the intermediate solution arrays

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

  return (0);
}

int apecss_bubble_freestruct(struct APECSS_Bubble *Bubble)
{
  free(Bubble->ODEsSol);
  Bubble->ODEsSol = NULL;
  free(Bubble->ODEsSolOld);
  Bubble->ODEsSolOld = NULL;
  free(Bubble->ode);
  Bubble->ode = NULL;

  free(Bubble->k2);
  Bubble->k2 = NULL;
  free(Bubble->k3);
  Bubble->k3 = NULL;
  free(Bubble->k4);
  Bubble->k4 = NULL;
  free(Bubble->k5);
  Bubble->k5 = NULL;
  free(Bubble->k6);
  Bubble->k6 = NULL;
  free(Bubble->k7);
  Bubble->k7 = NULL;
  free(Bubble->kLast);
  Bubble->kLast = NULL;

  if (Bubble->Emissions != NULL)
  {
    free(Bubble->Emissions);
    Bubble->Emissions = NULL;
  }

  if (Bubble->Results != NULL)
  {
    if (Bubble->Results->RayleighPlesset != NULL)
    {
      free(Bubble->Results->RayleighPlesset);
      Bubble->Results->RayleighPlesset = NULL;
    }

    if (Bubble->Results->Emissions != NULL)
    {
      if (Bubble->Results->Emissions->TimeInstances)
      {
        free(Bubble->Results->Emissions->TimeInstances);
        Bubble->Results->Emissions->TimeInstances = NULL;
      }

      if (Bubble->Results->Emissions->SpaceLocation != NULL)
      {
        free(Bubble->Results->Emissions->SpaceLocation);
        Bubble->Results->Emissions->SpaceLocation = NULL;
        Bubble->Results->Emissions->nSpaceLocations = 0;
      }

      if (Bubble->Results->Emissions->Node != NULL)
      {
        free(Bubble->Results->Emissions->Node);
        Bubble->Results->Emissions->Node = NULL;
        Bubble->Results->Emissions->nNodes = 0;
      }

      if (Bubble->Results->Emissions->Node_Rmin != NULL)
      {
        free(Bubble->Results->Emissions->Node_Rmin);
        Bubble->Results->Emissions->Node_Rmin = NULL;
      }

      if (Bubble->Results->Emissions->Node_Umin != NULL)
      {
        free(Bubble->Results->Emissions->Node_Umin);
        Bubble->Results->Emissions->Node_Umin = NULL;
      }

      if (Bubble->Results->Emissions->Node_pLmax != NULL)
      {
        free(Bubble->Results->Emissions->Node_pLmax);
        Bubble->Results->Emissions->Node_pLmax = NULL;
      }

      free(Bubble->Results->Emissions);
      Bubble->Results->Emissions = NULL;
    }

    free(Bubble->Results);
    Bubble->Results = NULL;
  }

  return (0);
}

// -------------------------------------------------------------------
// SOLVE
// -------------------------------------------------------------------
// Functions handling the solver for the bubble dynamics.
// -------------------------------------------------------------------

int apecss_bubble_solver_initialize(struct APECSS_Bubble *Bubble)
{
  // Initialize the counters
  Bubble->dtNumber = 0;
  Bubble->nSubIter = 0;

  // Set start time
  Bubble->t = Bubble->tStart;

  // Emissions (if applicable)
  Bubble->emissions_initialize(Bubble);

  // Store the initial results (if applicable)
  Bubble->results_rayleighplesset_store(Bubble);

  // Initialise the last-step solution of the RK54 scheme
  for (register int i = 0; i < Bubble->nODEs; i++) Bubble->kLast[i] = Bubble->ode[i](Bubble->ODEsSol, Bubble->tStart, Bubble);

  // Initialize the solution error
  Bubble->err = Bubble->NumericsODE->tol;

  // Initialize progress bar (if applicable)
  Bubble->progress_initial();
  Bubble->progress = 0;

  return (0);
}

int apecss_bubble_solver_run(APECSS_FLOAT tend, struct APECSS_Bubble *Bubble)
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

int apecss_bubble_solver_finalize(struct APECSS_Bubble *Bubble)
{
  Bubble->progress_final();
  Bubble->emissions_free(Bubble);

  return (0);
}

// -------------------------------------------------------------------
// PRESSURE AT INFINITY
// -------------------------------------------------------------------
// Functions defining the pressure at infinity and its derivatives.
// -------------------------------------------------------------------
// The functions are chosen in apecss_bubble_processoptions() and
// associated with the function pointers:
// - Bubble->get_pressure_infinity()
// - Bubble->get_pressurederivative_infinity()
// -------------------------------------------------------------------

APECSS_FLOAT apecss_bubble_pressure_infinity_noexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble) { return (Bubble->p0); }

APECSS_FLOAT apecss_bubble_pressure_infinity_sinexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (Bubble->p0 - Bubble->Excitation->dp * APECSS_SIN(2.0 * APECSS_PI * Bubble->Excitation->f * t));
}

APECSS_FLOAT apecss_bubble_pressurederivative_infinity_noexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble) { return (0.0); }

APECSS_FLOAT apecss_bubble_pressurederivative_infinity_sinexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble)
{
  return (-Bubble->Excitation->dp * 2.0 * APECSS_PI * Bubble->Excitation->f * APECSS_COS(2.0 * APECSS_PI * Bubble->Excitation->f * t));
}

// -------------------------------------------------------------------
// DIMENSIONAL RADIUS
// -------------------------------------------------------------------
// Functions defining the dimensional radius associated with the
// radial kinetic enthalpy of the liquid.
// -------------------------------------------------------------------
// The functions are chosen in apecss_bubble_processoptions() and
// associated with the function pointers:
// - Bubble->get_dimensionalradius()
// ------------------------------------------------------------------- +

APECSS_FLOAT apecss_bubble_dimensionalradius_planar(APECSS_FLOAT r) { return (1.0); }

APECSS_FLOAT apecss_bubble_dimensionalradius_cylindrical(APECSS_FLOAT r) { return (APECSS_SQRT(r)); }

APECSS_FLOAT apecss_bubble_dimensionalradius_spherical(APECSS_FLOAT r) { return (r); }

// -------------------------------------------------------------------
// PROGRESS SCREEN
// -------------------------------------------------------------------
// Functions handling the optional progress screen.
// -------------------------------------------------------------------
// The functions are chosen in apecss_bubble_processoptions() and
// associated with the function pointers:
// - Bubble->progress_initial()
// - Bubble->progress_update()
// - Bubble->progress_final()
// -------------------------------------------------------------------

int apecss_bubble_solver_progress_initialnone() { return (0); }

int apecss_bubble_solver_progress_initialscreen()
{
  fprintf(stderr, "| APECSS | Progress %%: ");
  return (0);
}

int apecss_bubble_solver_progress_updatenone(int *prog, APECSS_FLOAT elapsedtime, APECSS_FLOAT totaltime) { return (0); }

int apecss_bubble_solver_progress_updatescreen(int *prog, APECSS_FLOAT elapsedtime, APECSS_FLOAT totaltime)
{
  if (elapsedtime > (APECSS_FLOAT) *prog * 0.02 * totaltime)
  {
    if (!(*prog % 5))
      fprintf(stderr, "%i", *prog * 2);
    else
      fprintf(stderr, ".");

    (*prog)++;
  }

  return (0);
}

int apecss_bubble_solver_progress_finalnone() { return (0); }

int apecss_bubble_solver_progress_finalscreen()
{
  fprintf(stderr, "100\n");
  return (0);
}