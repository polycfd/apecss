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
// INITIALIZE / FREE
// -------------------------------------------------------------------

int apecss_emissions_initializestruct(struct APECSS_Bubble *Bubble)
{
  Bubble->Emissions = (struct APECSS_Emissions *) malloc(sizeof(struct APECSS_Emissions));
  Bubble->Emissions->Type = APECSS_EMISSION_NONE;
  Bubble->Emissions->CutOffDistance = 1.0e-3;
  Bubble->Emissions->KB_IterTolerance = 1.0;
  Bubble->Emissions->nNodes = 0;
  Bubble->Emissions->FirstNode = NULL;
  Bubble->Emissions->LastNode = NULL;

  return (0);
}

int apecss_emissions_initializenone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_emissions_initializelinkedlist(struct APECSS_Bubble *Bubble)
{
  Bubble->Emissions->FirstNode = NULL;
  Bubble->Emissions->LastNode = NULL;
  Bubble->Emissions->nNodes = 0;

  apecss_emissions_addnode(Bubble);

  return (0);
}

int apecss_emissions_freenone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_emissions_freelinkedlist(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Next = Bubble->Emissions->FirstNode;
  struct APECSS_EmissionNode *Current;

  while (Next != NULL)
  {
    Current = Next;
    Next = Current->forward;
    free(Current);
  }

  Bubble->Emissions->nNodes = 0;

  return (0);
}

// -------------------------------------------------------------------
// UPDATE
// -------------------------------------------------------------------

int apecss_emissions_updatenone(struct APECSS_Bubble *Bubble) { return (0); }

int apecss_emissions_updatelinkedlist(struct APECSS_Bubble *Bubble)
{
  if (Bubble->Emissions->nNodes) Bubble->Emissions->advance(Bubble);
  apecss_emissions_addnode(Bubble);
  if (Bubble->Emissions->LastNode->r > Bubble->Emissions->CutOffDistance) apecss_emissions_removenode(Bubble);

  return (0);
}

int apecss_emissions_addnode(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *New = NULL;
  New = (struct APECSS_EmissionNode *) malloc(sizeof(struct APECSS_EmissionNode));

  APECSS_FLOAT pL = Bubble->Liquid->get_pressure_bubblewall(Bubble->ODEsSol, Bubble->t, Bubble);
  APECSS_FLOAT pinf = Bubble->Liquid->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT rhoL = Bubble->Liquid->get_density(pL, Bubble->Liquid);
  APECSS_FLOAT rhoinf = Bubble->Liquid->get_density(pinf, Bubble->Liquid);

  New->g = Bubble->R *
           (Bubble->Liquid->get_enthalpy(pL, rhoL, Bubble->Liquid) - Bubble->Liquid->get_enthalpy(pinf, rhoinf, Bubble->Liquid) + 0.5 * APECSS_POW2(Bubble->U));
  New->f = APECSS_POW2(Bubble->R) * Bubble->U -
           Bubble->R * New->g / (Bubble->Liquid->get_soundspeed(pL, rhoL, Bubble->Liquid) + Bubble->Emissions->get_advectingvelocity(Bubble->U));

  New->r = Bubble->R;
  New->p = pL;
  New->u = Bubble->U;

  if (Bubble->Emissions->nNodes)
  {
    New->forward = Bubble->Emissions->FirstNode;
    Bubble->Emissions->FirstNode->backward = New;
  }
  else
  {
    New->forward = NULL;
    Bubble->Emissions->LastNode = New;
  }

  New->backward = NULL;
  Bubble->Emissions->FirstNode = New;

  Bubble->Emissions->nNodes += 1;
  New->id = Bubble->dtNumber;

  return (0);
}

int apecss_emissions_removenode(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Obsolete;

  do
  {
    Obsolete = Bubble->Emissions->LastNode;

    if (Bubble->Emissions->nNodes > 1)
    {
      Obsolete->backward->forward = NULL;
      Bubble->Emissions->LastNode = Obsolete->backward;
      Bubble->Emissions->nNodes -= 1;
      free(Obsolete);
    }
    else
    {
      Bubble->Emissions->FirstNode = NULL;
      Bubble->Emissions->LastNode = NULL;
      Bubble->Emissions->nNodes = 0;
      free(Obsolete);
      break;
    }
  } while (Bubble->Emissions->LastNode->r > Bubble->Emissions->CutOffDistance);

  return (0);
}

// -------------------------------------------------------------------
// ADVANCE EMISSION NODES
// -------------------------------------------------------------------

int apecss_emissions_advance_finitetimeincompressible(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Current = Bubble->Emissions->FirstNode;

  APECSS_FLOAT pinf = Bubble->Liquid->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT rho = Bubble->Liquid->rhoref;
  APECSS_FLOAT c = Bubble->Liquid->cref;
  APECSS_FLOAT dr = c * Bubble->dt;

  APECSS_FLOAT A = 2.0 * Bubble->R * APECSS_POW2(Bubble->U) + Bubble->ode[0](Bubble->ODEsSol, Bubble->t, Bubble) * APECSS_POW2(Bubble->R);
  APECSS_FLOAT B = APECSS_POW4(Bubble->R) * APECSS_POW2(Bubble->U);
  APECSS_FLOAT C = APECSS_POW2(Bubble->R) * Bubble->U;

  while (Current != NULL)
  {
    Current->r += dr;
    Current->p = pinf + rho * (A / Current->r - B / (2.0 * APECSS_POW4(Current->r)));
    Current->u = C / APECSS_POW2(Current->r);

    Bubble->results_emissionsnode_store(Current, c, pinf, Bubble);
    Current = Current->forward;
  }

  return (0);
}

int apecss_emissions_advance_quasiacoustic(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Current = Bubble->Emissions->FirstNode;

  APECSS_FLOAT pinf = Bubble->Liquid->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT rho = Bubble->Liquid->rhoref;
  APECSS_FLOAT c = Bubble->Liquid->cref;
  APECSS_FLOAT dr = c * Bubble->dt;

  while (Current != NULL)
  {
    Current->r += dr;
    Current->p = pinf + rho * (Current->g / Current->r - 0.5 * APECSS_POW2(Current->u));
    Current->u = Current->f / APECSS_POW2(Current->r) + Current->g / (Current->r * c);

    Bubble->results_emissionsnode_store(Current, c, pinf, Bubble);
    Current = Current->forward;
  }

  return (0);
}

int apecss_emissions_advance_kirkwoodbethe(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Current = Bubble->Emissions->LastNode;

  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT pinf = Bubble->Liquid->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT psi_inf = Gamma / Bubble->Liquid->get_density(pinf, Bubble->Liquid);
  APECSS_FLOAT dt = Bubble->dt;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  while (Current != NULL)
  {
    int NodeDiscarded = 0;
    APECSS_FLOAT rho, c, psi, Speed, p_prev;
    APECSS_FLOAT r_old = Current->r;  // Previous position of this node

    do
    {
      rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
      psi = Gamma / rho;
      c = Bubble->Liquid->get_soundspeed(Current->p, rho, Bubble->Liquid);
      Speed = c + Current->u;
      p_prev = Current->p;

      Current->r = r_old + Speed * dt;
      Current->u = Current->f / APECSS_POW2(Current->r) + Current->g / (Current->r * Speed);
      Current->p = ((psi_inf - b) * pinf + (psi_inf - psi) * B + (Gamma - 1.0) * (Current->g / Current->r - 0.5 * APECSS_POW2(Current->u))) / (psi - b);

      if (Current->p < -B)  // Check for unphysical pressure values
      {
        if (Current->forward != NULL)
          Current->forward->backward = Current->backward;
        else
          Bubble->Emissions->LastNode = Current->backward;

        if (Current->backward != NULL)
          Current->backward->forward = Current->forward;
        else
          Bubble->Emissions->FirstNode = Current->forward;

        struct APECSS_EmissionNode *Obsolete = Current;
        Current = Current->backward;
        free(Obsolete);
        Bubble->Emissions->nNodes -= 1;

        NodeDiscarded = 1;
        break;
      }
    } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

    if (!NodeDiscarded)
    {
      if (Current->forward != NULL && Current->r > Current->forward->r)  // Check for shock formation
      {
        Current->forward->backward = Current->backward;

        if (Current->backward != NULL)
          Current->backward->forward = Current->forward;
        else
          Bubble->Emissions->FirstNode = Current->forward;

        struct APECSS_EmissionNode *Obsolete = Current;
        Current = Current->backward;
        free(Obsolete);
        Bubble->Emissions->nNodes -= 1;
      }
      else if (Current->backward != NULL && Current->r < Current->backward->r)  // Check for shock formation
      {
        Current->backward->forward = Current->forward;

        if (Current->forward != NULL)
          Current->forward->backward = Current->backward;
        else
          Bubble->Emissions->LastNode = Current->backward;

        struct APECSS_EmissionNode *Obsolete = Current;
        Current = Current->backward;
        free(Obsolete);
        Bubble->Emissions->nNodes -= 1;
      }
      else if (Current->backward == NULL && Current->r < Bubble->R)  // Check for shock formation
      {
        Bubble->Emissions->FirstNode = Current->forward;

        if (Current->forward != NULL)
          Current->forward->backward = Current->backward;
        else
          Bubble->Emissions->LastNode = Current->backward;

        struct APECSS_EmissionNode *Obsolete = Current;
        Current = NULL;
        free(Obsolete);
        Bubble->Emissions->nNodes -= 1;
      }
      else
      {
        Bubble->results_emissionsnode_store(Current, c, pinf, Bubble);
        Current = Current->backward;
      }
    }
  }

  return (0);
}

// -------------------------------------------------------------------
// ADVECTING VELOCITY
// -------------------------------------------------------------------

APECSS_FLOAT apecss_emissions_getadvectingvelocity_returnzero(APECSS_FLOAT u) { return (0.0); }

APECSS_FLOAT apecss_emissions_getadvectingvelocity_returnvelocity(APECSS_FLOAT u) { return (u); }