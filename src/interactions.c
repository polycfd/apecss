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

int apecss_instantaneous_interactions(struct APECSS_Bubble *Bubble_1, struct APECSS_Bubble *Bubble_2)
{
  // Compute interactions by considering them instantaneous (fully incompressible approach)
  struct APECSS_Interaction *Interaction_1 = Bubble_1->Interaction;
  struct APECSS_Interaction *Interaction_2 = Bubble_2->Interaction;

  // Both bubbles are supposed to be in the same liquid
  struct APECSS_Liquid *Liquid = Bubble_1->Liquid;

  APECSS_FLOAT loc_1[3] = Interaction_1->location;
  APECSS_FLOAT loc_2[3] = Interaction_2->location;

  APECSS_FLOAT interbubble_dist = APECSS_SQRT(APECSS_POW2(loc_1[0] - loc_2[0]) + APECSS_POW2(loc_1[1] - loc_2[1]) + APECSS_POW2(loc_1[2] - loc_2[2]));

  APECSS_FLOAT dp_1_2 =
      Liquid->rhoref * ((2.0 * Bubble_1->R * APECSS_POW2(Bubble_1->U) + APECSS_POW2(Bubble_1->R) * Bubble_1->ode[0](Bubble_1->ODEsSol, Bubble_1->t, Bubble_1)) *
                            (1 / interbubble_dist) -
                        (APECSS_POW4(Bubble_1->R) * APECSS_POW2(Bubble_1->U)) * (1 / (2 * APECSS_POW4(interbubble_dist))));

  APECSS_FLOAT dp_2_1 =
      Liquid->rhoref * ((2.0 * Bubble_2->R * APECSS_POW2(Bubble_2->U) + APECSS_POW2(Bubble_2->R) * Bubble_2->ode[0](Bubble_2->ODEsSol, Bubble_2->t, Bubble_2)) *
                            (1 / interbubble_dist) -
                        (APECSS_POW4(Bubble_2->R) * APECSS_POW2(Bubble_2->U)) * (1 / (2 * APECSS_POW4(interbubble_dist))));

  Interaction_1->dp_neighbor += dp_2_1;
  Interaction_2->dp_neighbor += dp_1_2;

  return (0);
}

int apecss_quasi_acoustic_interactions(struct APECSS_Bubble *Bubbles[])
{
  // Compute interactions using emission nodes in the quasi acoustic assumption
  int nBubbles = Bubbles[0]->Interaction->nBubbles;

  struct APECSS_Interaction *Interactions[nBubbles];
  for (register int i = 0; i < nBubbles; i++)
  {
    Interactions[i] = Bubbles[i]->Interaction;
  }

  // All bubbles are supposed to be in the same liquid
  struct APECSS_Liquid *Liquid = Bubbles[0]->Liquid;

  for (register int i = 0; i < nBubbles; i++)
  {
    APECSS_FLOAT sumU = 0.0;
    APECSS_FLOAT sumG = 0.0;

    for (register int j = 0; j < nBubbles; j++)
    {
      if (j != i)
      {
        APECSS_FLOAT sumU_bubble = 0.0;
        APECSS_FLOAT sumG_bubble = 0.0;
        int nodecount_bubble = 0;

        APECSS_FLOAT loc_i[3] = Interactions[i]->location;
        APECSS_FLOAT loc_j[3] = Interactions[j]->location;

        APECSS_FLOAT interbubble_dist = APECSS_SQRT(APECSS_POW2(loc_i[0] - loc_j[0]) + APECSS_POW2(loc_i[1] - loc_j[1]) + APECSS_POW2(loc_i[2] - loc_j[2]));

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

        APECSS_FLOAT dp = 0.0;
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
    Interactions[i]->dp_neighbor += Liquid->rhoref * (sumG - 0.5 * APECSS_POW2(sumU));
  }

  return (0);
}