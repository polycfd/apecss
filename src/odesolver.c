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

#include "apecss.h"

// -------------------------------------------------------------------
// OPTIONS
// -------------------------------------------------------------------
// Functions initializing, processing and handling the options of the
// ODE solver.
// -------------------------------------------------------------------

int apecss_odesolver_setdefaultoptions(struct APECSS_NumericsODE *NumericsODE)
{
  NumericsODE->RKtype = APECSS_RK54_7M;
  NumericsODE->tol = 1.0e-10;
  NumericsODE->maxSubIter = 20;
  NumericsODE->dtMin = 1.0e-13;
  NumericsODE->dtMax = 1.0e-6;
  NumericsODE->minScale = 0.5;
  NumericsODE->maxScale = 2.0;
  NumericsODE->control_coeff_alpha = 0.9;
  NumericsODE->control_coeff_q = 0.1;

  return (0);
}

int apecss_odesolver_readoptions(struct APECSS_NumericsODE *NumericsODE, char *OptionsDir)
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

    if (strncasecmp(option, "odesolver", 9) == 0)
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
              NumericsODE->RKtype = APECSS_RK54_7M;
            }
            else if (strncasecmp(option3, "7S", 2) == 0)
            {
              NumericsODE->RKtype = APECSS_RK54_7S;
            }
          }
        }
        else if (strncasecmp(option2, "tolerance", 9) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          NumericsODE->tol = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "maxSubIterations", 16) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          NumericsODE->maxSubIter = atoi(option3);
        }
        else if (strncasecmp(option2, "minTimeStep", 11) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          NumericsODE->dtMin = APECSS_STRINGTOFLOAT(option3);
        }
        else if (strncasecmp(option2, "maxTimeStep", 11) == 0)
        {
          l = apecss_readoneoption(OptionsFile, option3);
          NumericsODE->dtMax = APECSS_STRINGTOFLOAT(option3);
        }
        else
        {
          sprintf(str, "An unknown option of ODESOLVER is given: %s, line %i", option2, line);
          apecss_erroronscreen(1, str);
        }
      }
    }
    else if (strncasecmp(option, "bubble", 6) == 0 || strncasecmp(option, "gas", 3) == 0 || strncasecmp(option, "liquid", 6) == 0 ||
             strncasecmp(option, "interface", 9) == 0 || strncasecmp(option, "results", 7) == 0)
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

int apecss_odesolver_processoptions(struct APECSS_NumericsODE *ODEs)
{
  if (ODEs->RKtype == APECSS_RK54_7M)
  {
    // RK5(4)7M coefficients of Dormand & Prince (1980), Table 2
    // RK matrix
    ODEs->a21 = 1.0 / 5.0;
    ODEs->a31 = 3.0 / 40.0;
    ODEs->a32 = 9.0 / 40.0;
    ODEs->a41 = 44.0 / 45.0;
    ODEs->a42 = -56.0 / 15.0;
    ODEs->a43 = 32.0 / 9.0;
    ODEs->a51 = 19372.0 / 6561.0;
    ODEs->a52 = -25360.0 / 2187.0;
    ODEs->a53 = 64448.0 / 6561.0;
    ODEs->a54 = -212.0 / 729.0;
    ODEs->a61 = 9017.0 / 3168.0;
    ODEs->a62 = -355.0 / 33.0;
    ODEs->a63 = 46732.0 / 5247.0;
    ODEs->a64 = 49.0 / 176.0;
    ODEs->a65 = -5103.0 / 18656.0;
    ODEs->a71 = 35.0 / 384.0;
    ODEs->a72 = 0.0;
    ODEs->a73 = 500.0 / 1113.0;
    ODEs->a74 = 125.0 / 192.0;
    ODEs->a75 = -2187.0 / 6784.0;
    ODEs->a76 = 11.0 / 84.0;

    // 5th-order weights
    ODEs->b1 = 35.0 / 384.0;
    ODEs->b3 = 500.0 / 1113.0;
    ODEs->b4 = 125.0 / 192.0;
    ODEs->b5 = -2187.0 / 6784.0;
    ODEs->b6 = 11.0 / 84.0;
    ODEs->b7 = 0.0;

    // 4th-order weights
    ODEs->bs1 = 5179.0 / 57600.0;
    ODEs->bs3 = 7571.0 / 16695.0;
    ODEs->bs4 = 393.0 / 640.0;
    ODEs->bs5 = -92097.0 / 339200.0;
    ODEs->bs6 = 187.0 / 2100.0;
    ODEs->bs7 = 1.0 / 40.0;

    // Nodes
    ODEs->c2 = 1.0 / 5.0;
    ODEs->c3 = 3.0 / 10.0;
    ODEs->c4 = 4.0 / 5.0;
    ODEs->c5 = 8.0 / 9.0;
    ODEs->c6 = 1.0;
    ODEs->c7 = 1.0;
  }
  else if (ODEs->RKtype == APECSS_RK54_7S)
  {
    // RK5(4)7S coefficients of Dormand & Prince (1980), Table 3
    // RK matrix
    ODEs->a21 = 2.0 / 9.0;
    ODEs->a31 = 1.0 / 12.0;
    ODEs->a32 = 1.0 / 4.0;
    ODEs->a41 = 55.0 / 324.0;
    ODEs->a42 = -25.0 / 108.0;
    ODEs->a43 = 50.0 / 81.0;
    ODEs->a51 = 83.0 / 330.0;
    ODEs->a52 = -13.0 / 22.0;
    ODEs->a53 = 61.0 / 66.0;
    ODEs->a54 = 9.0 / 110.0;
    ODEs->a61 = -19.0 / 28.0;
    ODEs->a62 = 9.0 / 4.0;
    ODEs->a63 = 1.0 / 7.0;
    ODEs->a64 = -27.0 / 7.0;
    ODEs->a65 = 22.0 / 7.0;
    ODEs->a71 = 19.0 / 200.0;
    ODEs->a72 = 0.0;
    ODEs->a73 = 3.0 / 5.0;
    ODEs->a74 = -243.0 / 400.0;
    ODEs->a75 = 33.0 / 40.0;
    ODEs->a76 = 7.0 / 80.0;

    // 5th-order weights
    ODEs->b1 = 19.0 / 200.0;
    ODEs->b3 = 3.0 / 5.0;
    ODEs->b4 = -243.0 / 400.0;
    ODEs->b5 = 33.0 / 40.0;
    ODEs->b6 = 7.0 / 80.0;
    ODEs->b7 = 0.0;

    // 4th-order weights
    ODEs->bs1 = 431.0 / 5000.0;
    ODEs->bs3 = 333.0 / 500.0;
    ODEs->bs4 = -7857.0 / 10000.0;
    ODEs->bs5 = 957.0 / 1000.0;
    ODEs->bs6 = 193.0 / 2000.0;
    ODEs->bs7 = -1.0 / 50.0;

    // Nodes
    ODEs->c2 = 2.0 / 9.0;
    ODEs->c3 = 1.0 / 3.0;
    ODEs->c4 = 5.0 / 9.0;
    ODEs->c5 = 2.0 / 3.0;
    ODEs->c6 = 1.0;
    ODEs->c7 = 1.0;
  }

  // Embedded error coefficients
  ODEs->e1 = ODEs->b1 - ODEs->bs1;
  ODEs->e3 = ODEs->b3 - ODEs->bs3;
  ODEs->e4 = ODEs->b4 - ODEs->bs4;
  ODEs->e5 = ODEs->b5 - ODEs->bs5;
  ODEs->e6 = ODEs->b6 - ODEs->bs6;
  ODEs->e7 = ODEs->b7 - ODEs->bs7;

  return (0);
}

// -------------------------------------------------------------------
// SOLVER
// -------------------------------------------------------------------
// Function with the actual ODE solver, based on the RK5(4) embedded
// Runge-Kutta scheme of Dormand and Prince (1980).
// -------------------------------------------------------------------

APECSS_FLOAT apecss_odesolver(struct APECSS_Bubble *Bubble)
{
  APECSS_FLOAT *SolTemp, t;
  SolTemp = malloc(Bubble->nODEs * sizeof(APECSS_FLOAT));

  // Step 1 ---------------------------------------------------
  // (First Same As Last)

  // Step 2 ---------------------------------------------------
  t = Bubble->t + Bubble->NumericsODE->c2 * Bubble->dt;

  for (register int i = 0; i < Bubble->nODEs; i++) SolTemp[i] = Bubble->ODEsSol[i] + Bubble->dt * Bubble->NumericsODE->a21 * Bubble->kLast[i];

  Bubble->k2[0] = Bubble->ode[0](SolTemp, t, Bubble);
  Bubble->k2[1] = SolTemp[0];
  for (register int i = 2; i < Bubble->nODEs; i++) Bubble->k2[i] = Bubble->ode[i](SolTemp, t, Bubble);

  // Step 3 ---------------------------------------------------
  t = Bubble->t + Bubble->NumericsODE->c3 * Bubble->dt;

  for (register int i = 0; i < Bubble->nODEs; i++)
    SolTemp[i] = Bubble->ODEsSol[i] + Bubble->dt * (Bubble->NumericsODE->a31 * Bubble->kLast[i] + Bubble->NumericsODE->a32 * Bubble->k2[i]);

  Bubble->k3[0] = Bubble->ode[0](SolTemp, t, Bubble);
  Bubble->k3[1] = SolTemp[0];
  for (register int i = 2; i < Bubble->nODEs; i++) Bubble->k3[i] = Bubble->ode[i](SolTemp, t, Bubble);

  // Step 4 ---------------------------------------------------
  t = Bubble->t + Bubble->NumericsODE->c4 * Bubble->dt;

  for (register int i = 0; i < Bubble->nODEs; i++)
    SolTemp[i] = Bubble->ODEsSol[i] + Bubble->dt * (Bubble->NumericsODE->a41 * Bubble->kLast[i] + Bubble->NumericsODE->a42 * Bubble->k2[i] +
                                                    Bubble->NumericsODE->a43 * Bubble->k3[i]);

  Bubble->k4[0] = Bubble->ode[0](SolTemp, t, Bubble);
  Bubble->k4[1] = SolTemp[0];
  for (register int i = 2; i < Bubble->nODEs; i++) Bubble->k4[i] = Bubble->ode[i](SolTemp, t, Bubble);

  // Step 5 ---------------------------------------------------
  t = Bubble->t + Bubble->NumericsODE->c5 * Bubble->dt;

  for (register int i = 0; i < Bubble->nODEs; i++)
    SolTemp[i] = Bubble->ODEsSol[i] + Bubble->dt * (Bubble->NumericsODE->a51 * Bubble->kLast[i] + Bubble->NumericsODE->a52 * Bubble->k2[i] +
                                                    Bubble->NumericsODE->a53 * Bubble->k3[i] + Bubble->NumericsODE->a54 * Bubble->k4[i]);

  Bubble->k5[0] = Bubble->ode[0](SolTemp, t, Bubble);
  Bubble->k5[1] = SolTemp[0];
  for (register int i = 2; i < Bubble->nODEs; i++) Bubble->k5[i] = Bubble->ode[i](SolTemp, t, Bubble);

  // Step 6 ---------------------------------------------------
  t = Bubble->t + Bubble->NumericsODE->c6 * Bubble->dt;

  for (register int i = 0; i < Bubble->nODEs; i++)
    SolTemp[i] = Bubble->ODEsSol[i] +
                 Bubble->dt * (Bubble->NumericsODE->a61 * Bubble->kLast[i] + Bubble->NumericsODE->a62 * Bubble->k2[i] +
                               Bubble->NumericsODE->a63 * Bubble->k3[i] + Bubble->NumericsODE->a64 * Bubble->k4[i] + Bubble->NumericsODE->a65 * Bubble->k5[i]);

  Bubble->k6[0] = Bubble->ode[0](SolTemp, t, Bubble);
  Bubble->k6[1] = SolTemp[0];
  for (register int i = 2; i < Bubble->nODEs; i++) Bubble->k6[i] = Bubble->ode[i](SolTemp, t, Bubble);

  // Step 7 ---------------------------------------------------
  t = Bubble->t + Bubble->NumericsODE->c7 * Bubble->dt;

  for (register int i = 0; i < Bubble->nODEs; i++)
    SolTemp[i] = Bubble->ODEsSol[i] + Bubble->dt * (Bubble->NumericsODE->a71 * Bubble->kLast[i] + Bubble->NumericsODE->a72 * Bubble->k2[i] +
                                                    Bubble->NumericsODE->a73 * Bubble->k3[i] + Bubble->NumericsODE->a74 * Bubble->k4[i] +
                                                    Bubble->NumericsODE->a75 * Bubble->k5[i] + Bubble->NumericsODE->a76 * Bubble->k6[i]);

  Bubble->k7[0] = Bubble->ode[0](SolTemp, t, Bubble);
  Bubble->k7[1] = SolTemp[0];
  for (register int i = 2; i < Bubble->nODEs; i++) Bubble->k7[i] = Bubble->ode[i](SolTemp, t, Bubble);

  free(SolTemp);

  // New solution ---------------------------------------------
  for (register int i = 0; i < Bubble->nODEs; i++)
    Bubble->ODEsSol[i] +=
        Bubble->dt * (Bubble->NumericsODE->b1 * Bubble->kLast[i] + Bubble->NumericsODE->b3 * Bubble->k3[i] + Bubble->NumericsODE->b4 * Bubble->k4[i] +
                      Bubble->NumericsODE->b5 * Bubble->k5[i] + Bubble->NumericsODE->b6 * Bubble->k6[i] + Bubble->NumericsODE->b7 * Bubble->k7[i]);

  // Compute solution error -----------------------------------
  APECSS_FLOAT err = 0.0;
  for (register int i = 0; i < Bubble->nODEs; i++)
    err = APECSS_MAX(err,
                     APECSS_ABS(Bubble->NumericsODE->e1 * Bubble->kLast[i] + Bubble->NumericsODE->e3 * Bubble->k3[i] + Bubble->NumericsODE->e4 * Bubble->k4[i] +
                                Bubble->NumericsODE->e5 * Bubble->k5[i] + Bubble->NumericsODE->e6 * Bubble->k6[i] + Bubble->NumericsODE->e7 * Bubble->k7[i]) *
                         Bubble->dt / (APECSS_ABS(Bubble->ODEsSol[i]) + APECSS_SMALL));

  return (err);
}

// -------------------------------------------------------------------
// TIME STEP
// -------------------------------------------------------------------
// Function setting the time-step based on the solution error of the
// previous solve.
// -------------------------------------------------------------------

int apecss_odesolver_settimestep(struct APECSS_NumericsODE *ODEs, APECSS_FLOAT err, APECSS_FLOAT timetoend, APECSS_FLOAT *dt)
{
  // Scale by which the current time step size is multiplied
  APECSS_FLOAT scale = (ODEs->control_coeff_alpha) * APECSS_POW(ODEs->tol / err, ODEs->control_coeff_q);

  // Apply upper and lower limit to scale
  scale = APECSS_MIN(APECSS_MAX(scale, ODEs->minScale), ODEs->maxScale);

  // Apply scale **/
  *dt *= scale;

  // Apply upper and lower limit to time step size
  *dt = APECSS_MIN3(APECSS_MAX(*dt, ODEs->dtMin), ODEs->dtMax, timetoend);

  return (0);
}
