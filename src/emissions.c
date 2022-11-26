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
  Bubble->Emissions->Scheme = APECSS_EMISSION_INTEGRATE_EULER;
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
  // Allocate the new node
  struct APECSS_EmissionNode *New = NULL;
  New = (struct APECSS_EmissionNode *) malloc(sizeof(struct APECSS_EmissionNode));

  // Values used multiple times
  APECSS_FLOAT pL = Bubble->Liquid->get_pressure_bubblewall(Bubble->ODEsSol, Bubble->t, Bubble);
  APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT rhoL = Bubble->Liquid->get_density(pL, Bubble->Liquid);
  APECSS_FLOAT rhoinf = Bubble->Liquid->get_density(pinf, Bubble->Liquid);
  APECSS_FLOAT hL = Bubble->Liquid->get_enthalpy(pL, rhoL, Bubble->Liquid);

  // Compute the invariants f and g
  New->g = Bubble->R * (hL - Bubble->Liquid->get_enthalpy(pinf, rhoinf, Bubble->Liquid) + 0.5 * APECSS_POW2(Bubble->U));
  New->f = APECSS_POW2(Bubble->R) * Bubble->U -
           Bubble->R * New->g / (Bubble->Liquid->get_soundspeed(pL, rhoL, Bubble->Liquid) + Bubble->Emissions->get_advectingvelocity(Bubble->U));

  // Apply the boundary conditions
  New->r = Bubble->R;
  New->u = Bubble->U;
  New->h = hL;
  New->p = pL;

  // Set the neighbor information
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

  // Count the new node in
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

  // Values used multiple times
  APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT rho = Bubble->Liquid->rhoref;
  APECSS_FLOAT c = Bubble->Liquid->cref;
  APECSS_FLOAT dr = c * Bubble->dt;

  // Define constants as shortcuts
  APECSS_FLOAT A = 2.0 * Bubble->R * APECSS_POW2(Bubble->U) + Bubble->ode[0](Bubble->ODEsSol, Bubble->t, Bubble) * APECSS_POW2(Bubble->R);
  APECSS_FLOAT B = APECSS_POW4(Bubble->R) * APECSS_POW2(Bubble->U);
  APECSS_FLOAT C = APECSS_POW2(Bubble->R) * Bubble->U;

  while (Current != NULL)
  {
    Current->r += dr;
    Current->u = C / APECSS_POW2(Current->r);
    Current->p = pinf + rho * (A / Current->r - B / (2.0 * APECSS_POW4(Current->r)));

    // Store data (if applicable)
    Bubble->results_emissionsnode_store(Current, c, pinf, Bubble);

    // Move to the next node
    Current = Current->forward;
  }

  return (0);
}

int apecss_emissions_advance_quasiacoustic(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Current = Bubble->Emissions->FirstNode;

  // Values used multiple times
  APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT rho = Bubble->Liquid->rhoref;
  APECSS_FLOAT c = Bubble->Liquid->cref;
  APECSS_FLOAT dr = c * Bubble->dt;

  while (Current != NULL)
  {
    Current->r += dr;
    Current->p = pinf + rho * (Current->g / Current->r - 0.5 * APECSS_POW2(Current->u));
    Current->u = Current->f / APECSS_POW2(Current->r) + Current->g / (Current->r * c);

    // Store data (if applicable)
    Bubble->results_emissionsnode_store(Current, c, pinf, Bubble);

    // Move to the next node
    Current = Current->forward;
  }

  return (0);
}

int apecss_emissions_advance_kirkwoodbethe_tait(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Current = Bubble->Emissions->LastNode;

  // Constants used multiple times
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT h_fac = (Gamma - 1.0) * Bubble->Liquid->rhoref / (Gamma * APECSS_POW(Bubble->Liquid->pref + B, (1.0 / Gamma)));
  APECSS_FLOAT h_exp = 1.0 / (1.0 - (1.0 / Gamma));

  // Assuming the time-step is small compared to the timescale over which the ambient/driving pressure changes,
  // we take the reference state at infinity as constant.
  APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT hinf = apecss_liquid_enthalpy_nasg(pinf, Bubble->Liquid->get_density(pinf, Bubble->Liquid), Bubble->Liquid);

  while (Current != NULL)
  {
    if (Bubble->Emissions->integrate_along_characteristic(Bubble, Current, hinf))
    {
      if (Current->forward != NULL && Current->r > Current->forward->r)  // Check for shock formation (outward moving node)
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
      else if (Current->backward != NULL && Current->r < Current->backward->r)  // Check for shock formation (inward moving node)
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
      else if (Current->backward == NULL && Current->r < Bubble->R)  // Check for shock formation (inward moving node adjacent to the bubble wall)
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
        // Evaluate pressure
        Current->p = APECSS_POW(h_fac * Current->h, h_exp) - B;

        // Store data (if applicable)
        Bubble->results_emissionsnode_store(Current, APECSS_SQRT((Gamma - 1.0) * Current->h), pinf, Bubble);

        // Move to the next node
        Current = Current->backward;
      }
    }
    else  // Unphysical enthalpy obtained during integration, node is removed
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
    }
  }

  return (0);
}

int apecss_emissions_advance_kirkwoodbethe_general(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Current = Bubble->Emissions->LastNode;

  // Assuming the time-step is small compared to the timescale over which the ambient/driving pressure changes,
  // we take the reference state at infinity as constant.
  APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT hinf = apecss_liquid_enthalpy_nasg(pinf, Bubble->Liquid->get_density(pinf, Bubble->Liquid), Bubble->Liquid);

  while (Current != NULL)
  {
    if (Bubble->Emissions->integrate_along_characteristic(Bubble, Current, hinf))
    {
      if (Current->forward != NULL && Current->r > Current->forward->r)  // Check for shock formation (outward moving node)
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
      else if (Current->backward != NULL && Current->r < Current->backward->r)  // Check for shock formation (inward moving node)
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
      else if (Current->backward == NULL && Current->r < Bubble->R)  // Check for shock formation (inward moving node adjacent to the bubble wall)
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
        // Pressure has already been evaluated iteratively while integrating along the characteristic

        // Store data (if applicable)
        Bubble->results_emissionsnode_store(
            Current, Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid), pinf, Bubble);

        // Move to the next node
        Current = Current->backward;
      }
    }
    else
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
    }
  }

  return (0);
}

// -------------------------------------------------------------------
// INTEGRATION OF POSITION AND VELOCITY ALONG CHARACTERISTIC
// -------------------------------------------------------------------

int apecss_emissions_integrate_ekb_tait_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Simplified method based on the Kirkwood-Bethe hypothesis of Denner & Schenke, with u treated explicitly.

  APECSS_FLOAT c = APECSS_SQRT((Bubble->Liquid->Gamma - 1.0) * Current->h);
  APECSS_FLOAT u = Current->u;

  Current->r += Bubble->dt * (c + u);
  Current->u = Current->f / APECSS_POW2(Current->r) + Current->g / (Current->r * (c + Current->u));
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  if (Current->h < 0)
    return (0);
  else
    return (1);
}

int apecss_emissions_integrate_ekb_tait_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Simplified method based on the Kirkwood-Bethe hypothesis of Denner & Schenke, with u treated explicitly.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT dt = Bubble->dt;

  APECSS_FLOAT c = APECSS_SQRT((Gamma - 1.0) * Current->h);

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c + u;

  APECSS_FLOAT h = (hinf + Current->g / r - 0.5 * APECSS_POW2(u));
  if (h < 0) return (0);

  c = APECSS_SQRT((Gamma - 1.0) * h);

  r = Current->r + 0.5 * dt * kr1;
  u = Current->f / APECSS_POW2(r) + Current->g / (r * (c + u));
  APECSS_FLOAT kr2 = c + u;

  h = (hinf + Current->g / r - 0.5 * APECSS_POW2(u));
  if (h < 0) return (0);

  c = APECSS_SQRT((Gamma - 1.0) * h);

  r = Current->r + 0.5 * dt * kr2;
  u = Current->f / APECSS_POW2(r) + Current->g / (r * (c + u));
  APECSS_FLOAT kr3 = c + u;

  h = (hinf + Current->g / r - 0.5 * APECSS_POW2(u));
  if (h < 0) return (0);

  c = APECSS_SQRT((Gamma - 1.0) * h);

  r = Current->r + dt * kr3;
  u = Current->f / APECSS_POW2(r) + Current->g / (r * (c + u));
  APECSS_FLOAT kr4 = c + u;

  Current->r += dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->u = Current->f / APECSS_POW2(Current->r) + Current->g / (Current->r * (c + Current->u));
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  if (Current->h < 0)
    return (0);
  else
    return (1);
}

int apecss_emissions_integrate_gfc_tait_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Fully-compressible second-order approximation of Gilmore (1952), see Section IV.C, with integration of u over r.

  APECSS_FLOAT c = APECSS_SQRT((Bubble->Liquid->Gamma - 1.0) * Current->h);
  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;

  APECSS_FLOAT dr = Bubble->dt * (c + u);
  Current->r += dr;
  Current->u += dr * (-2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u)));
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  return (1);
}

int apecss_emissions_integrate_gfc_tait_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Fully-compressible second-order approximation of Gilmore (1952), see Section IV.C, with integration of u over r.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT dt = Bubble->dt;

  APECSS_FLOAT c = APECSS_SQRT((Gamma - 1.0) * Current->h);

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c + u;
  APECSS_FLOAT ku1 = -2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / r - 0.5 * APECSS_POW2(u)));

  r = Current->r + 0.5 * dt * kr1;
  u = Current->u + 0.5 * dt * kr1 * ku1;
  APECSS_FLOAT kr2 = c + u;
  APECSS_FLOAT ku2 = -2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / r - 0.5 * APECSS_POW2(u)));

  r = Current->r + 0.5 * dt * kr2;
  u = Current->u + 0.5 * dt * kr2 * ku2;
  APECSS_FLOAT kr3 = c + u;
  APECSS_FLOAT ku3 = -2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / r - 0.5 * APECSS_POW2(u)));

  r = Current->r + dt * kr3;
  u = Current->u + dt * kr3 * ku3;
  APECSS_FLOAT kr4 = c + u;
  APECSS_FLOAT ku4 = -2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u));

  APECSS_FLOAT dr = dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->r += dr;
  Current->u += dr * (ku1 + 2.0 * (ku2 + ku3) + ku4) * APECSS_ONESIXTH;
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  return (1);
}

int apecss_emissions_integrate_hpe_tait_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Method of Hickling & Plesset (1963), see Appendix I, and Ebeling (1978), with integration of u over t.

  APECSS_FLOAT c = APECSS_SQRT((Bubble->Liquid->Gamma - 1.0) * Current->h);
  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;

  Current->r += Bubble->dt * (c + u);
  Current->u += Bubble->dt * (((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u)));
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  return (1);
}

int apecss_emissions_integrate_hpe_tait_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Method of Hickling & Plesset (1963), see Appendix I, and Ebeling (1978), with integration of u over t.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT dt = Bubble->dt;

  APECSS_FLOAT c = APECSS_SQRT((Gamma - 1.0) * Current->h);

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c + u;
  APECSS_FLOAT ku1 = ((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / r - 0.5 * APECSS_POW2(u)));

  r = Current->r + 0.5 * dt * kr1;
  u = Current->u + 0.5 * dt * ku1;
  APECSS_FLOAT kr2 = c + u;
  APECSS_FLOAT ku2 = ((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / r - 0.5 * APECSS_POW2(u)));

  r = Current->r + 0.5 * dt * kr2;
  u = Current->u + 0.5 * dt * ku2;
  APECSS_FLOAT kr3 = c + u;
  APECSS_FLOAT ku3 = ((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / r - 0.5 * APECSS_POW2(u)));

  r = Current->r + dt * kr3;
  u = Current->u + dt * ku3;
  APECSS_FLOAT kr4 = c + u;
  APECSS_FLOAT ku4 = ((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u));

  Current->r += dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->u += dt * (ku1 + 2.0 * (ku2 + ku3) + ku4) * APECSS_ONESIXTH;
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  return (1);
}

int apecss_emissions_integrate_ekb_general_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Simplified method based on the Kirkwood-Bethe hypothesis of Denner & Schenke, with u treated explicitly.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);
  APECSS_FLOAT u = Current->u;

  Current->r += Bubble->dt * (c + u);
  Current->u = Current->f / APECSS_POW2(Current->r) + Current->g / (Current->r * (c + Current->u));
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  APECSS_FLOAT p_prev;
  do
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (Current->p < -B) return (0);  // Check for unphysical pressure values
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

int apecss_emissions_integrate_ekb_general_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Simplified method based on the Kirkwood-Bethe hypothesis of Denner & Schenke, with u treated explicitly.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT dt = Bubble->dt;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c + u;

  APECSS_FLOAT h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  APECSS_FLOAT p = Current->p, p_prev;
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (p < -B) return (0);  // Check for unphysical pressure values
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));

  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + 0.5 * dt * kr1;
  u = Current->f / APECSS_POW2(r) + Current->g / (r * (c + u));
  APECSS_FLOAT kr2 = c + u;

  h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  do
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (p < -B) return (0);  // Check for unphysical pressure values
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));

  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + 0.5 * dt * kr2;
  u = Current->f / APECSS_POW2(r) + Current->g / (r * (c + u));
  APECSS_FLOAT kr3 = c + u;

  h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  do
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (p < -B) return (0);  // Check for unphysical pressure values
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));

  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + dt * kr3;
  u = Current->f / APECSS_POW2(r) + Current->g / (r * (c + u));
  APECSS_FLOAT kr4 = c + u;

  Current->r += dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->u = Current->f / APECSS_POW2(Current->r) + Current->g / (Current->r * (c + Current->u));
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  do
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (Current->p < -B) return (0);  // Check for unphysical pressure values
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

int apecss_emissions_integrate_gfc_general_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Fully-compressible second-order approximation of Gilmore (1952), see Section IV.C, with integration of u over r.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);
  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;

  APECSS_FLOAT dr = Bubble->dt * (c + u);
  Current->r += dr;
  Current->u += dr * (-2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u)));
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  APECSS_FLOAT p_prev;
  do
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

int apecss_emissions_integrate_gfc_general_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Fully-compressible second-order approximation of Gilmore (1952), see Section IV.C, with integration of u over r.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT dt = Bubble->dt;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c + u;
  APECSS_FLOAT ku1 = -2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u));

  APECSS_FLOAT h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  APECSS_FLOAT p = Current->p, p_prev;
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + 0.5 * dt * kr1;
  u = Current->u + 0.5 * dt * kr1 * ku1;
  APECSS_FLOAT kr2 = c + u;
  APECSS_FLOAT ku2 = -2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u));

  h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  do
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + 0.5 * dt * kr2;
  u = Current->u + 0.5 * dt * kr2 * ku2;
  APECSS_FLOAT kr3 = c + u;
  APECSS_FLOAT ku3 = -2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u));

  h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  do
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + dt * kr3;
  u = Current->u + dt * kr3 * ku3;
  APECSS_FLOAT kr4 = c + u;
  APECSS_FLOAT ku4 = -2.0 * (u + APECSS_POW3(u) / (APECSS_POW2(c) - APECSS_POW2(u))) / r + Current->g / (APECSS_POW2(r) * (c - u));

  APECSS_FLOAT dr = dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->r += dr;
  Current->u += dr * (ku1 + 2.0 * (ku2 + ku3) + ku4) * APECSS_ONESIXTH;
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  do
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

int apecss_emissions_integrate_hpe_general_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Method of Hickling & Plesset (1963), see Appendix I, and Ebeling (1978), with integration of u over t.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);
  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;

  Current->r += Bubble->dt * (c + u);
  Current->u += Bubble->dt * (((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u)));
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  APECSS_FLOAT p_prev;
  do
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

int apecss_emissions_integrate_hpe_general_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Method of Hickling & Plesset (1963), see Appendix I, and Ebeling (1978), with integration of u over t.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT dt = Bubble->dt;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c + u;
  APECSS_FLOAT ku1 = ((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u));

  APECSS_FLOAT h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  APECSS_FLOAT p = Current->p, p_prev;
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + 0.5 * dt * kr1;
  u = Current->u + 0.5 * dt * ku1;
  APECSS_FLOAT kr2 = c + u;
  APECSS_FLOAT ku2 = ((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u));

  h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  do
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + 0.5 * dt * kr2;
  u = Current->u + 0.5 * dt * ku2;
  APECSS_FLOAT kr3 = c + u;
  APECSS_FLOAT ku3 = ((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u));

  h = hinf + Current->g / r - 0.5 * APECSS_POW2(u);
  do
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  r = Current->r + dt * kr3;
  u = Current->u + dt * ku3;
  APECSS_FLOAT kr4 = c + u;
  APECSS_FLOAT ku4 = ((c + u) * Current->g / r - 2.0 * APECSS_POW2(c) * u) / (r * (c - u));

  Current->r += dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->u += dt * (ku1 + 2.0 * (ku2 + ku3) + ku4) * APECSS_ONESIXTH;
  Current->h = hinf + Current->g / Current->r - 0.5 * APECSS_POW2(Current->u);

  do
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

// -------------------------------------------------------------------
// ADVECTING VELOCITY
// -------------------------------------------------------------------

APECSS_FLOAT apecss_emissions_getadvectingvelocity_returnzero(APECSS_FLOAT u) { return (0.0); }

APECSS_FLOAT apecss_emissions_getadvectingvelocity_returnvelocity(APECSS_FLOAT u) { return (u); }