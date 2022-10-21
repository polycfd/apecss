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

#include <sys/stat.h>
#include "apecss.h"

int apecss_options_setdefault(struct APECSS_Bubble* Bubble)
{
  // Call all functions to set default values
  apecss_gas_setdefaultoptions(Bubble);
  apecss_interface_setdefaultoptions(Bubble);
  apecss_liquid_setdefaultoptions(Bubble);
  apecss_bubble_setdefaultoptions(Bubble);

  return (0);
}

int apecss_options_process(struct APECSS_Bubble* Bubble)
{
  // Call all functions to process the provided options
  apecss_gas_processoptions(Bubble);
  apecss_interface_processoptions(Bubble);
  apecss_liquid_processoptions(Bubble);
  apecss_bubble_processoptions(Bubble);

  return (0);
}

int apecss_options_readfile(struct APECSS_Bubble* Bubble, char* OptionsDir)
{
  int l = 0;
  int line = 0;
  FILE* OptionsFile;
  char str[APECSS_STRINGLENGTH_SPRINTF];
  char option[APECSS_STRINGLENGTH], option2[APECSS_STRINGLENGTH], option3[APECSS_STRINGLENGTH];
  int StatusFile = 1;
  int StatusSection = 1;

  if ((OptionsFile = fopen(OptionsDir, "r")) == (FILE*) NULL)
  {
    sprintf(str, "File %s cannot be opened for reading.\n", OptionsDir);
    apecss_erroronscreen(1, str);
  }

  while ((l = apecss_readoneoption(OptionsFile, option)) != EOF && StatusFile == 1)
  {
    line += l;

    /****************************************************************/
    /* BUBBLE                                                       */
    /****************************************************************/
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

            if (strncasecmp(option3, "incompressible", 14) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_INCOMPRESSIBLE;
            }
            else if (strncasecmp(option3, "fti", 3) == 0 || strncasecmp(option3, "finitetimeincompressible", 24) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_FINITE_TIME_INCOMPRESSIBLE;
            }
            else if (strncasecmp(option3, "qa", 2) == 0 || strncasecmp(option3, "quasiacoustic", 13) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_QUASIACOUSTIC;
            }
            else if (strncasecmp(option3, "kb", 2) == 0 || strncasecmp(option3, "kirkwoodbethe", 13) == 0)
            {
              Bubble->Emissions->Type = APECSS_EMISSION_KIRKWOODBETHE;
            }

            l = apecss_readoneoption(OptionsFile, option3);
            Bubble->Emissions->CutOffDistance = APECSS_STRINGTOFLOAT(option3);
          }
        }
        else if (strncasecmp(option2, "kbitertolerance", 15) == 0)
        {
          if (Bubble->Emissions == NULL) apecss_emissions_initializestruct(Bubble);

          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Emissions->KB_IterTolerance = APECSS_STRINGTOFLOAT(option3);
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
        else
        {
          sprintf(str, "An unknown option of BUBBLE is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    /****************************************************************/
    /* GAS                                                          */
    /****************************************************************/
    else if (strncasecmp(option, "gas", 3) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else if (strncasecmp(option2, "eos", 3) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "ig", 2) == 0)
            {
              Bubble->Gas->EoS = APECSS_GAS_IG;
            }
            else if (strncasecmp(option3, "hc", 2) == 0)
            {
              Bubble->Gas->EoS = APECSS_GAS_HC;
            }
            else if (strncasecmp(option3, "nasg", 4) == 0)
            {
              Bubble->Gas->EoS = APECSS_GAS_NASG;
            }
          }
        }
        else if (strncasecmp(option2, "polytropicexponent", 18) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Gas->Gamma = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencepressure", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Gas->pref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencedensity", 16) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Gas->rhoref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "covolume", 8) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Gas->b = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "hardcoreradius", 14) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Gas->h = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "taitpressureconst", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Gas->B = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "molecularweight", 15) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Gas->mmol = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "moleculardiameter", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Gas->dmol = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of GAS is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    /****************************************************************/
    /* LIQUID                                                       */
    /****************************************************************/
    else if (strncasecmp(option, "liquid", 6) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else if (strncasecmp(option2, "liquidtype", 10) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "newtonian", 9) == 0)
            {
              Bubble->Liquid->Type = APECSS_LIQUID_NEWTONIAN;
            }
            else if (strncasecmp(option3, "kelvinvoigt", 11) == 0)
            {
              Bubble->Liquid->Type = APECSS_LIQUID_KELVINVOIGT;
            }
            else if (strncasecmp(option3, "zener", 5) == 0)
            {
              Bubble->Liquid->Type = APECSS_LIQUID_ZENER;
            }
            else if (strncasecmp(option3, "oldroydb", 8) == 0)
            {
              Bubble->Liquid->Type = APECSS_LIQUID_OLDROYDB;
            }
          }
        }
        else if (strncasecmp(option2, "eos", 3) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "tait", 4) == 0)
            {
              Bubble->Liquid->EoS = APECSS_LIQUID_TAIT;
            }
            else if (strncasecmp(option3, "nasg", 4) == 0)
            {
              Bubble->Liquid->EoS = APECSS_LIQUID_NASG;
            }
          }
        }
        else if (strncasecmp(option2, "polytropicexponent", 18) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->Gamma = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencepressure", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->pref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencedensity", 16) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->rhoref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "referencesoundspeed", 19) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->cref = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "covolume", 8) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->b = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "taitpressureconst", 17) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->B = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "viscosity", 9) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->mu = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "shearmodulus", 12) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->G = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "polymerviscosity", 16) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->eta = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "relaxationtime", 14) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Liquid->lambda = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of LIQUID is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    /****************************************************************/
    /* INTERFACE                                                    */
    /****************************************************************/
    else if (strncasecmp(option, "interface", 9) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else if (strncasecmp(option2, "surfacetensioncoeff", 19) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Interface->sigma = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "lipidcoatingmodel", 17) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "none", 4) == 0)
            {
              Bubble->Interface->LipidCoatingModel = APECSS_LIPIDCOATING_NONE;
            }
            else if (strncasecmp(option3, "marmottant", 10) == 0)
            {
              Bubble->Interface->LipidCoatingModel = APECSS_LIPIDCOATING_MARMOTTANT;
            }
            else if (strncasecmp(option3, "gompertz-marmottant", 19) == 0)
            {
              Bubble->Interface->LipidCoatingModel = APECSS_LIPIDCOATING_MARMOTTANT + APECSS_LIPIDCOATING_GOMPERTZFUNCTION;
            }
          }
        }
        else if (strncasecmp(option2, "sigmainit", 9) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Interface->sigma0 = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "elasticity", 10) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Interface->Elasticity = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "dilatationalviscosity", 21) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->Interface->Viscosity = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of INTERFACE is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    /****************************************************************/
    /* RESULTS                                                      */
    /****************************************************************/
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
            APECSS_FLOAT* temp = malloc(Bubble->Results->Emissions->nTimeInstances * sizeof(APECSS_FLOAT));
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
            struct APECSS_ResultsEmissionsSpace* temp;
            temp = (struct APECSS_ResultsEmissionsSpace*) malloc(Bubble->Results->Emissions->nSpaceLocations * sizeof(struct APECSS_ResultsEmissionsSpace));
            for (register int l = 0; l < Bubble->Results->Emissions->nSpaceLocations; l++) temp[l] = Bubble->Results->Emissions->SpaceLocation[l];
            free(Bubble->Results->Emissions->SpaceLocation);

            Bubble->Results->Emissions->nSpaceLocations++;
            Bubble->Results->Emissions->SpaceLocation =
                (struct APECSS_ResultsEmissionsSpace*) malloc(Bubble->Results->Emissions->nSpaceLocations * sizeof(struct APECSS_ResultsEmissionsSpace));
            for (register int l = 0; l < Bubble->Results->Emissions->nSpaceLocations - 1; l++) Bubble->Results->Emissions->SpaceLocation[l] = temp[l];
            free(temp);
          }
          else
          {
            Bubble->Results->Emissions->nSpaceLocations++;
            Bubble->Results->Emissions->SpaceLocation =
                (struct APECSS_ResultsEmissionsSpace*) malloc(Bubble->Results->Emissions->nSpaceLocations * sizeof(struct APECSS_ResultsEmissionsSpace));
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
            struct APECSS_ResultsEmissionsNode* temp;
            temp = (struct APECSS_ResultsEmissionsNode*) malloc(Bubble->Results->Emissions->nNodes * sizeof(struct APECSS_ResultsEmissionsNode));
            for (register int l = 0; l < Bubble->Results->Emissions->nNodes; l++) temp[l] = Bubble->Results->Emissions->Node[l];
            free(Bubble->Results->Emissions->Node);

            Bubble->Results->Emissions->nNodes++;
            Bubble->Results->Emissions->Node =
                (struct APECSS_ResultsEmissionsNode*) malloc(Bubble->Results->Emissions->nNodes * sizeof(struct APECSS_ResultsEmissionsNode));
            for (register int l = 0; l < Bubble->Results->Emissions->nNodes - 1; l++) Bubble->Results->Emissions->Node[l] = temp[l];
            free(temp);
          }
          else
          {
            Bubble->Results->Emissions->nNodes++;
            Bubble->Results->Emissions->Node =
                (struct APECSS_ResultsEmissionsNode*) malloc(Bubble->Results->Emissions->nNodes * sizeof(struct APECSS_ResultsEmissionsNode));
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
          Bubble->Results->Emissions->Node_Rmin = (struct APECSS_ResultsEmissionsNode*) malloc(sizeof(struct APECSS_ResultsEmissionsNode));
          apecss_results_emissionsnode_initializenode(Bubble->Results->Emissions->Node_Rmin);

          Bubble->Results->Emissions->Umin = APECSS_LARGE;
          Bubble->Results->Emissions->Node_Umin = (struct APECSS_ResultsEmissionsNode*) malloc(sizeof(struct APECSS_ResultsEmissionsNode));
          apecss_results_emissionsnode_initializenode(Bubble->Results->Emissions->Node_Umin);

          Bubble->Results->Emissions->pLmax = -APECSS_LARGE;
          Bubble->Results->Emissions->Node_pLmax = (struct APECSS_ResultsEmissionsNode*) malloc(sizeof(struct APECSS_ResultsEmissionsNode));
          apecss_results_emissionsnode_initializenode(Bubble->Results->Emissions->Node_pLmax);
        }
        else
        {
          sprintf(str, "An unknown option of RESULTS is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    /****************************************************************/
    /* ODESOLVER                                                    */
    /****************************************************************/
    else if (strncasecmp(option, "odesolver", 9) == 0)
    {
      StatusSection = 1;
      while (StatusSection == 1 && (l = apecss_readoneoption(OptionsFile, option2)) != EOF)
      {
        line += l;
        if (strncasecmp(option2, "END", 3) == 0)
        {
          StatusSection = 0;
        }
        else if (strncasecmp(option2, "RK", 2) == 0)
        {
          if ((l = apecss_readoneoption(OptionsFile, option3)) == EOF)
          {
            StatusSection = 0;
            StatusFile = 0;
          }
          else
          {
            line += l - 1;
            if (strncasecmp(option3, "7M", 2) == 0)
            {
              Bubble->NumericsODE->RKtype = APECSS_RK54_7M;
            }
            else if (strncasecmp(option3, "7S", 2) == 0)
            {
              Bubble->NumericsODE->RKtype = APECSS_RK54_7S;
            }
          }
        }
        else if (strncasecmp(option2, "tolerance", 9) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->NumericsODE->tol = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "maxSubIterations", 16) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->NumericsODE->maxSubIter = atoi(option3);
        }
        else if (strncasecmp(option2, "minTimeStep", 11) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->NumericsODE->dtMin = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "maxTimeStep", 11) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          Bubble->NumericsODE->dtMax = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of ODESOLVER is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    /****************************************************************/
    /*                                                              */
    /****************************************************************/
    else
    {
      sprintf(str, "An unknown Section is given: %s, line %i", option, line);
      apecss_erroronscreen(1, str);
    }
  }

  fclose(OptionsFile);

  return (0);
}

int apecss_readoneoption(FILE* OptionsFile, char* option)
{
  char ch, ch2;
  int line;
  char str[APECSS_STRINGLENGTH_SPRINTF];

  line = 1;
  while ((ch = getc(OptionsFile)) != EOF)
  {
    if (ch != '\n' && ch != '#')
    {
      while (ch == ' ' || ch == '\t' || ch == '\n')
      {
        ch = getc(OptionsFile);
      }

      memcpy(option, &ch, sizeof(char));
      ch2 = getc(OptionsFile);
      ungetc(ch2, OptionsFile);

      if (ch2 != ' ' && ch2 != '\n')
      {
        apecss_linegetskip(str, OptionsFile);
        memcpy(option + 1, str, sizeof(str));
      }
      else
      {
        str[0] = '\0';
        option[1] = '\0';
      }

      return (line);
    }

    if (ch != '\n')
    {
      apecss_lineget(str, OptionsFile);
    }

    line++;
  }

  return (EOF);
}

int apecss_lineget(char* str, FILE* OptionsFile)
{
  char ch;
  int index = 0;

  while ((ch = getc(OptionsFile)) != EOF)
  {
    str[index++] = ch;
    if (index == APECSS_STRINGLENGTH_SPRINTF) apecss_erroronscreen(1, "String in file too long!");

    if (ch == '\n')
    {
      str[index] = '\0';
      return (index - 1);
    }
  }

  return (0);
}

int apecss_linegetskip(char* str, FILE* OptionsFile)
{
  char ch;
  int index = 0;

  while ((ch = getc(OptionsFile)) != EOF)
  {
    if (ch != ' ' && ch != '\n')
    {
      str[index++] = ch;
      if (index == APECSS_STRINGLENGTH_SPRINTF) apecss_erroronscreen(1, "String in file too long!");
    }
    if ((ch == '\n' || ch == ' ') && index > 0)
    {
      str[index] = '\0';
      return (index);
    }
  }

  return (0);
}