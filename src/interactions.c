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
// Functions computing different types of interactions between cavita-
// tion bubbles inside a cluster
// -------------------------------------------------------------------

int apecss_interactions_instantaneous(struct APECSS_Bubble *Bubbles[])
{
  // Compute interactions by considering them instantaneous (fully incompressible approach)
  int nBubbles = Bubbles[0]->Interaction->nBubbles;

  // All bubbles are supposed to be in the same liquid
  struct APECSS_Liquid *Liquid = Bubbles[0]->Liquid;

  // Preliminary step to gather bubble wall acceleration values
  // APECSS_FLOAT Bubbles_A[nBubbles];
  // for (register int i = 0; i < nBubbles; i++)
  // {
  //   APECSS_FLOAT acceleration = Bubbles[i]->ode[0](Bubbles[i]->ODEsSol, Bubbles[i]->t, Bubbles[i]);
  //   Bubbles_A[i] = acceleration;
  // }

  // // Reset pressure contributions of the neighbours
  // for (register int i = 0; i < nBubbles; i++)
  // {
  //   Bubbles[i]->Interaction->dp_neighbor = 0.0;
  // }

  for (register int i = 0; i < nBubbles; i++)
  {
    APECSS_FLOAT x_i = Bubbles[i]->Interaction->location[0];
    APECSS_FLOAT y_i = Bubbles[i]->Interaction->location[1];
    APECSS_FLOAT z_i = Bubbles[i]->Interaction->location[2];

    Bubbles[i]->Interaction->dp_neighbor = 0.0;

    for (register int j = 0; j < nBubbles; j++)
    {
      if (j != i)
      {
        APECSS_FLOAT x_j = Bubbles[j]->Interaction->location[0];
        APECSS_FLOAT y_j = Bubbles[j]->Interaction->location[1];
        APECSS_FLOAT z_j = Bubbles[j]->Interaction->location[2];

        APECSS_FLOAT interbubble_dist = APECSS_SQRT(APECSS_POW2(x_i - x_j) + APECSS_POW2(y_i - y_j) + APECSS_POW2(z_i - z_j));
        // printf("Bubble %d vs Bubble %d, %e %e %e", i, j, interbubble_dist, Bubbles[i]->ode[0](Bubbles[i]->ODEsSol, Bubbles[i]->t, Bubbles[i]),
        //        Bubbles[j]->ode[0](Bubbles[j]->ODEsSol, Bubbles[j]->t, Bubbles[j]));
        // printf("Bubble %d vs Bubble %d, %e %e %e", i, j, interbubble_dist, Bubbles_A[i], Bubbles_A[j]);

        // APECSS_FLOAT dp = Liquid->rhoref * ((2.0 * Bubbles[j]->R * APECSS_POW2(Bubbles[j]->U) +
        //                                      APECSS_POW2(Bubbles[j]->R) * Bubbles[j]->ode[0](Bubbles[j]->ODEsSol, Bubbles[j]->t, Bubbles[j])) *
        //                                         (1 / interbubble_dist) -
        //                                     (APECSS_POW4(Bubbles[j]->R) * APECSS_POW2(Bubbles[j]->U)) * (1 / (2 * APECSS_POW4(interbubble_dist))));

        // APECSS_FLOAT dp =
        //     Liquid->rhoref * (((2.0 * Bubbles[j]->R * APECSS_POW2(Bubbles[j]->U) + APECSS_POW2(Bubbles[j]->R) * Bubbles_A[j]) * (1 / interbubble_dist)) -
        //                       ((APECSS_POW4(Bubbles[j]->R) * APECSS_POW2(Bubbles[j]->U)) * (1 / (2 * APECSS_POW4(interbubble_dist)))));

        APECSS_FLOAT dp = Liquid->rhoref * 2.0 * Bubbles[j]->R * APECSS_POW2(Bubbles[j]->U) * (1 / interbubble_dist);
        dp += Liquid->rhoref * APECSS_POW2(Bubbles[j]->R) * Bubbles[j]->ode[0](Bubbles[j]->ODEsSol, Bubbles[j]->t, Bubbles[j]) * (1 / interbubble_dist);
        dp += Liquid->rhoref * APECSS_POW4(Bubbles[j]->R) * APECSS_POW2(Bubbles[j]->U) * (1 / (2 * APECSS_POW4(interbubble_dist)));

        // if ((i == 0) && (j == 1))
        //   printf(" %e %e %e", 2.0 * Bubbles[j]->R * APECSS_POW2(Bubbles[j]->U) * (1 / interbubble_dist),
        //          APECSS_POW2(Bubbles[j]->R) * Bubbles_A[j] * (1 / interbubble_dist),
        //          (APECSS_POW4(Bubbles[j]->R) * APECSS_POW2(Bubbles[j]->U)) * (1 / (2 * APECSS_POW4(interbubble_dist))));

        Bubbles[i]->Interaction->dp_neighbor += dp;
        // printf(" %e %e\n", Bubbles[i]->ode[0](Bubbles[i]->ODEsSol, Bubbles[i]->t, Bubbles[i]),
        //        Bubbles[j]->ode[0](Bubbles[j]->ODEsSol, Bubbles[j]->t, Bubbles[j]));
        // printf(" %e\n", Bubbles[i]->Interaction->dp_neighbor);
      }
    }
  }

  return (0);
}

int apecss_interactions_quasi_acoustic(struct APECSS_Bubble *Bubbles[])
{
  // Compute interactions using emission nodes in the quasi acoustic assumption
  int nBubbles = Bubbles[0]->Interaction->nBubbles;

  // All bubbles are supposed to be in the same liquid
  struct APECSS_Liquid *Liquid = Bubbles[0]->Liquid;

  // Reset pressure contributions of the neighbours
  for (register int i = 0; i < nBubbles; i++)
  {
    Bubbles[i]->Interaction->dp_neighbor = 0.0;
  }

  for (register int i = 0; i < nBubbles; i++)
  {
    APECSS_FLOAT sumU = 0.0;
    APECSS_FLOAT sumG = 0.0;

    APECSS_FLOAT x_i = Bubbles[i]->Interaction->location[0];
    APECSS_FLOAT y_i = Bubbles[i]->Interaction->location[1];
    APECSS_FLOAT z_i = Bubbles[i]->Interaction->location[2];

    for (register int j = 0; j < nBubbles; j++)
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

        APECSS_FLOAT r_b = Bubbles[i]->R;

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

        if (nodecount_bubble != 0)
        {
          APECSS_FLOAT inv_nodecount = 1 / (APECSS_FLOAT) nodecount_bubble;
          APECSS_FLOAT sumU_temp = sumU_bubble * inv_nodecount;
          APECSS_FLOAT sumG_temp = sumG_bubble * inv_nodecount;

          sumU += (sumU_temp);
          sumG += (sumG_temp);
        }
      }
    }
    Bubbles[i]->Interaction->dp_neighbor += Liquid->rhoref * (sumG - 0.5 * APECSS_POW2(sumU));
    // printf(" %e\n", Bubbles[i]->Interaction->dp_neighbor);
  }

  return (0);
}

int apecss_interactions_cutoffdistance(struct APECSS_Bubble *Bubbles[])
{
  // Compute the proper cut off distance for each bubbles in a cluster depending on the location of their neighbors
  int nBubbles = Bubbles[0]->Interaction->nBubbles;

  for (register int i = 0; i < nBubbles; i++)
  {
    APECSS_FLOAT x_i = Bubbles[i]->Interaction->location[0];
    APECSS_FLOAT y_i = Bubbles[i]->Interaction->location[1];
    APECSS_FLOAT z_i = Bubbles[i]->Interaction->location[2];
    APECSS_FLOAT dist = 0.0;

    for (register int j = 0; j < nBubbles; j++)
    {
      APECSS_FLOAT x_j = Bubbles[j]->Interaction->location[0];
      APECSS_FLOAT y_j = Bubbles[j]->Interaction->location[1];
      APECSS_FLOAT z_j = Bubbles[j]->Interaction->location[2];

      APECSS_FLOAT dist_ij = APECSS_SQRT(APECSS_POW2(x_i - x_j) + APECSS_POW2(y_i - y_j) + APECSS_POW2(z_i - z_j));

      if (dist_ij > dist) dist = dist_ij;
    }
    if (Bubbles[i]->Emissions != NULL) Bubbles[i]->Emissions->CutOffDistance = 1.05 * dist;
  }

  return (0);
}