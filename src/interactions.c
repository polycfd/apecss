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

int apecss_instantaneous_interactions(struct APECSS_Bubble *Bubbles[])
{
  // Compute interactions uby considering them instantaneous (fully incompressible approach)
  int nBubbles = Bubbles[0]->Interaction->nBubbles;

  // All bubbles are supposed to be in the same liquid
  struct APECSS_Liquid *Liquid = Bubbles[0]->Liquid;

  for (register int i = 0; i < nBubbles; i++)
  {
    APECSS_FLOAT x_i = Bubbles[i]->Interaction->location[0];
    APECSS_FLOAT y_i = Bubbles[i]->Interaction->location[1];
    APECSS_FLOAT z_i = Bubbles[i]->Interaction->location[2];

    for (register int j = 0; j < nBubbles; j++)
    {
      if (j != i)
      {
        APECSS_FLOAT x_j = Bubbles[j]->Interaction->location[0];
        APECSS_FLOAT y_j = Bubbles[j]->Interaction->location[1];
        APECSS_FLOAT z_j = Bubbles[j]->Interaction->location[2];

        APECSS_FLOAT interbubble_dist = APECSS_SQRT(APECSS_POW2(x_i - x_j) + APECSS_POW2(y_i - y_j) + APECSS_POW2(z_i - z_j));

        APECSS_FLOAT dp = Liquid->rhoref * ((2.0 * Bubbles[j]->R * APECSS_POW2(Bubbles[j]->U) +
                                             APECSS_POW2(Bubbles[j]->R) * Bubbles[j]->ode[0](Bubbles[j]->ODEsSol, Bubbles[j]->t, Bubbles[j])) *
                                                (1 / interbubble_dist) -
                                            (APECSS_POW4(Bubbles[j]->R) * APECSS_POW2(Bubbles[j]->U)) * (1 / (2 * APECSS_POW4(interbubble_dist))));

        Bubbles[i]->Interaction->dp_neighbor += dp;
      }
    }
  }

  return (0);
}

int apecss_quasi_acoustic_interactions(struct APECSS_Bubble *Bubbles[])
{
  // Compute interactions using emission nodes in the quasi acoustic assumption
  int nBubbles = Bubbles[0]->Interaction->nBubbles;

  // All bubbles are supposed to be in the same liquid
  struct APECSS_Liquid *Liquid = Bubbles[0]->Liquid;

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
  }

  return (0);
}