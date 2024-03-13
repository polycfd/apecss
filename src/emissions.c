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
// INITIALIZE / FREE
// -------------------------------------------------------------------
// Functions dealing with the memory management of the emissions.
// -------------------------------------------------------------------

int apecss_emissions_initializestruct(struct APECSS_Bubble *Bubble)
{
  Bubble->Emissions = (struct APECSS_Emissions *) malloc(sizeof(struct APECSS_Emissions));
  Bubble->Emissions->Type = APECSS_EMISSION_NONE;
  Bubble->Emissions->Scheme = APECSS_EMISSION_INTEGRATE_RK4;
  Bubble->Emissions->CutOffDistance = 1.0e-3;
  Bubble->Emissions->KB_IterTolerance = 1.0;
  Bubble->Emissions->nNodes = 0;
  Bubble->Emissions->FirstNode = NULL;
  Bubble->Emissions->LastNode = NULL;
  Bubble->Emissions->advance = NULL;
  Bubble->Emissions->compute_f = NULL;
  Bubble->Emissions->integrate_along_characteristic = NULL;

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
// Functions updating the linked list of emission nodes.
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
  APECSS_FLOAT rhoL = Bubble->Liquid->get_density(pL, Bubble->Liquid);
  APECSS_FLOAT hL = Bubble->Liquid->get_enthalpy(pL, rhoL, Bubble->Liquid);
  APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT hinf = Bubble->Liquid->get_enthalpy(pinf, Bubble->Liquid->get_density(pinf, Bubble->Liquid), Bubble->Liquid);

  // Compute the invariants f and g
  New->g = Bubble->get_dimensionalradius(Bubble->R) * (hL - hinf + 0.5 * APECSS_POW2(Bubble->U));
  New->f = Bubble->Emissions->compute_f(Bubble, New->g, pL, rhoL);

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
// Functions advancing the emission nodes in space dependent on the
// chosen model for the emissions.
// -------------------------------------------------------------------
// The appropriate function is chosen in apecss_bubble_processoptions()
// and hooked up to the function pointer Bubble->Emissions->advance().
// -------------------------------------------------------------------

int apecss_emissions_advance_finitespeedincompressible(struct APECSS_Bubble *Bubble)
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
    Current->u = Current->f / APECSS_POW2(Current->r);
    Current->p = pinf + rho * (Current->g / Current->r - 0.5 * APECSS_POW2(Current->u));

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
    Current->u = Current->f / APECSS_POW2(Current->r) + Current->g / (Current->r * c);
    Current->p = pinf + rho * (Current->g / Current->r - 0.5 * APECSS_POW2(Current->u));

    // Store data (if applicable)
    Bubble->results_emissionsnode_store(Current, c, pinf, Bubble);

    // Move to the next node
    Current = Current->forward;
  }

  return (0);
}

// int apecss_emissions_advance_kirkwoodbethe_tait(struct APECSS_Bubble *Bubble)
// {
//   struct APECSS_EmissionNode *Current = Bubble->Emissions->LastNode;

//   // Constants used multiple times
//   APECSS_FLOAT B = Bubble->Liquid->B;
//   APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
//   APECSS_FLOAT h_fac = (Gamma - 1.0) * Bubble->Liquid->rhoref / (Gamma * APECSS_POW(Bubble->Liquid->pref + B, (1.0 / Gamma)));
//   APECSS_FLOAT h_exp = 1.0 / (1.0 - (1.0 / Gamma));

//   // Assuming the time-step is small compared to the timescale over which the ambient/driving pressure changes,
//   // we take the reference state at infinity as constant.
//   APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
//   APECSS_FLOAT hinf = apecss_liquid_enthalpy_nasg(pinf, Bubble->Liquid->get_density(pinf, Bubble->Liquid), Bubble->Liquid);

//   while (Current != NULL)
//   {
//     if (Bubble->Emissions->integrate_along_characteristic(Bubble, Current, hinf))
//     {
//       if (Current->forward != NULL && Current->r > Current->forward->r)  // Check for shock formation (outward moving node)
//       {
//         Current->forward->backward = Current->backward;

//         if (Current->backward != NULL)
//           Current->backward->forward = Current->forward;
//         else
//           Bubble->Emissions->FirstNode = Current->forward;

//         struct APECSS_EmissionNode *Obsolete = Current;
//         Current = Current->backward;
//         free(Obsolete);
//         Bubble->Emissions->nNodes -= 1;
//       }
//       else if (Current->backward != NULL && Current->r < Current->backward->r)  // Check for shock formation (inward moving node)
//       {
//         Current->backward->forward = Current->forward;

//         if (Current->forward != NULL)
//           Current->forward->backward = Current->backward;
//         else
//           Bubble->Emissions->LastNode = Current->backward;

//         struct APECSS_EmissionNode *Obsolete = Current;
//         Current = Current->backward;
//         free(Obsolete);
//         Bubble->Emissions->nNodes -= 1;
//       }
//       else if (Current->backward == NULL && Current->r < Bubble->R)  // Check for shock formation (inward moving node adjacent to the bubble wall)
//       {
//         Bubble->Emissions->FirstNode = Current->forward;

//         if (Current->forward != NULL)
//           Current->forward->backward = Current->backward;
//         else
//           Bubble->Emissions->LastNode = Current->backward;

//         struct APECSS_EmissionNode *Obsolete = Current;
//         Current = NULL;
//         free(Obsolete);
//         Bubble->Emissions->nNodes -= 1;
//       }
//       else
//       {
//         // Evaluate pressure
//         Current->p = APECSS_POW(h_fac * Current->h, h_exp) - B;

//         // Store data (if applicable)
//         Bubble->results_emissionsnode_store(Current, APECSS_SQRT((Gamma - 1.0) * Current->h), pinf, Bubble);

//         // Move to the next node
//         Current = Current->backward;
//       }
//     }
//     else  // Physically implausible enthalpy/pressure detected, node is discarded
//     {
//       if (Current->forward != NULL)
//         Current->forward->backward = Current->backward;
//       else
//         Bubble->Emissions->LastNode = Current->backward;

//       if (Current->backward != NULL)
//         Current->backward->forward = Current->forward;
//       else
//         Bubble->Emissions->FirstNode = Current->forward;

//       struct APECSS_EmissionNode *Obsolete = Current;
//       Current = Current->backward;
//       free(Obsolete);
//       Bubble->Emissions->nNodes -= 1;
//     }
//   }

//   return (0);
// }

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

  // ---------------------------------------
  // Integrate along the outgoing characteristics

  while (Current != NULL)
  {
    if (Bubble->Emissions->integrate_along_characteristic(Bubble, Current, hinf))
    {
      // Evaluate pressure
      Current->p = APECSS_POW(h_fac * Current->h, h_exp) - B;

      // Move to the next node
      Current = Current->backward;
    }
    else  // Physically implausible enthalpy/pressure detected, node is discarded
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

  // ---------------------------------------
  // Shock treatment (if necessary)

  // Current = Bubble->Emissions->FirstNode;
  // APECSS_FLOAT mass0 = 0.0;
  // while (Current != NULL)
  // {
  //   if (Current->forward != NULL)
  //   {
  //     APECSS_FLOAT rhoavg = 0.5 * (Bubble->Liquid->get_density(Current->p, Bubble->Liquid) + Bubble->Liquid->get_density(Current->forward->p, Bubble->Liquid));
  //     mass0 += rhoavg * (APECSS_POW3(Current->forward->r) - APECSS_POW3(Current->r));
  //   }

  //   Current = Current->forward;
  // }

  int check_for_shocks = 1;
  while (check_for_shocks)
  {
    Current = Bubble->Emissions->FirstNode;
    check_for_shocks = 0;

    while (Current != NULL)
    {
      if (Current->forward != NULL && Current->r > Current->forward->r)  // Check for shock formation
      {
        check_for_shocks = 1;  // Requires another complete subsequent iteration to make sure all multivalued solutions have been or are treated

        // Define new properties of the forward node
        Current->forward->r = 0.5 * (Current->r + Current->forward->r);
        Current->forward->u = 0.5 * (Current->u + Current->forward->u);
        Current->forward->g = 0.5 * (Current->g + Current->forward->g);
        Current->forward->h = hinf + Current->forward->g / Bubble->get_dimensionalradius(Current->forward->r) - 0.5 * APECSS_POW2(Current->forward->u);
        Current->forward->p = APECSS_POW(h_fac * Current->h, h_exp) - B;
        Current->forward->f = Bubble->Emissions->compute_f(Bubble, Current->g, Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid));

        // Current node is obsolete and will be deleted
        struct APECSS_EmissionNode *Obsolete = Current;

        // Redefine the neighbor information
        Current->forward->backward = Current->backward;
        if (Current->backward != NULL)
          Current->backward->forward = Current->forward;
        else
          Bubble->Emissions->FirstNode = Current->forward;

        // Move to the next node
        Current = Current->forward;

        // Delete obsolete node
        free(Obsolete);
        Bubble->Emissions->nNodes -= 1;
      }
      else
      {
        // Move to the next node without action
        Current = Current->forward;
      }
    }
  }

  // Current = Bubble->Emissions->FirstNode;
  // APECSS_FLOAT mass1 = 0.0;
  // while (Current != NULL)
  // {
  //   if (Current->forward != NULL)
  //   {
  //     APECSS_FLOAT rhoavg = 0.5 * (Bubble->Liquid->get_density(Current->p, Bubble->Liquid) + Bubble->Liquid->get_density(Current->forward->p, Bubble->Liquid));
  //     mass1 += rhoavg * (APECSS_POW3(Current->forward->r) - APECSS_POW3(Current->r));
  //   }

  //   Current = Current->forward;
  // }

  // ---------------------------------------
  // Store data (if applicable)

  if (Bubble->Results->Emissions->nNodes + Bubble->Results->Emissions->MinMaxPeriod)
  {
    Current = Bubble->Emissions->FirstNode;

    while (Current != NULL)
    {
      Bubble->results_emissionsnode_store(Current, APECSS_SQRT((Gamma - 1.0) * Current->h), pinf, Bubble);
      Current = Current->forward;
    }
  }

  return (0);
}

// int apecss_emissions_advance_kirkwoodbethe_general(struct APECSS_Bubble *Bubble)
// {
//   struct APECSS_EmissionNode *Current = Bubble->Emissions->LastNode;

//   // Assuming the time-step is small compared to the timescale over which the ambient/driving pressure changes,
//   // we take the reference state at infinity as constant.
//   APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
//   APECSS_FLOAT hinf = apecss_liquid_enthalpy_nasg(pinf, Bubble->Liquid->get_density(pinf, Bubble->Liquid), Bubble->Liquid);

//   while (Current != NULL)
//   {
//     if (Bubble->Emissions->integrate_along_characteristic(Bubble, Current, hinf))
//     {
//       if (Current->forward != NULL && Current->r > Current->forward->r)  // Check for shock formation (outward moving node)
//       {
//         Current->forward->backward = Current->backward;

//         if (Current->backward != NULL)
//           Current->backward->forward = Current->forward;
//         else
//           Bubble->Emissions->FirstNode = Current->forward;

//         struct APECSS_EmissionNode *Obsolete = Current;
//         Current = Current->backward;
//         free(Obsolete);
//         Bubble->Emissions->nNodes -= 1;
//       }
//       else if (Current->backward != NULL && Current->r < Current->backward->r)  // Check for shock formation (inward moving node)
//       {
//         Current->backward->forward = Current->forward;

//         if (Current->forward != NULL)
//           Current->forward->backward = Current->backward;
//         else
//           Bubble->Emissions->LastNode = Current->backward;

//         struct APECSS_EmissionNode *Obsolete = Current;
//         Current = Current->backward;
//         free(Obsolete);
//         Bubble->Emissions->nNodes -= 1;
//       }
//       else if (Current->backward == NULL && Current->r < Bubble->R)  // Check for shock formation (inward moving node adjacent to the bubble wall)
//       {
//         Bubble->Emissions->FirstNode = Current->forward;

//         if (Current->forward != NULL)
//           Current->forward->backward = Current->backward;
//         else
//           Bubble->Emissions->LastNode = Current->backward;

//         struct APECSS_EmissionNode *Obsolete = Current;
//         Current = NULL;
//         free(Obsolete);
//         Bubble->Emissions->nNodes -= 1;
//       }
//       else
//       {
//         // Pressure has already been evaluated iteratively while integrating along the characteristic

//         // Store data (if applicable)
//         Bubble->results_emissionsnode_store(
//             Current, Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid), pinf, Bubble);

//         // Move to the next node
//         Current = Current->backward;
//       }
//     }
//     else  // Physically implausible enthalpy/pressure detected, node is discarded
//     {
//       if (Current->forward != NULL)
//         Current->forward->backward = Current->backward;
//       else
//         Bubble->Emissions->LastNode = Current->backward;

//       if (Current->backward != NULL)
//         Current->backward->forward = Current->forward;
//       else
//         Bubble->Emissions->FirstNode = Current->forward;

//       struct APECSS_EmissionNode *Obsolete = Current;
//       Current = Current->backward;
//       free(Obsolete);
//       Bubble->Emissions->nNodes -= 1;
//     }
//   }

//   return (0);
// }

int apecss_emissions_advance_kirkwoodbethe_general(struct APECSS_Bubble *Bubble)
{
  struct APECSS_EmissionNode *Current = Bubble->Emissions->LastNode;

  // Assuming the time-step is small compared to the timescale over which the ambient/driving pressure changes,
  // we take the reference state at infinity as constant.
  APECSS_FLOAT pinf = Bubble->get_pressure_infinity(Bubble->t, Bubble);
  APECSS_FLOAT hinf = apecss_liquid_enthalpy_nasg(pinf, Bubble->Liquid->get_density(pinf, Bubble->Liquid), Bubble->Liquid);

  // ---------------------------------------
  // Integrate along the outgoing characteristic

  while (Current != NULL)
  {
    if (Bubble->Emissions->integrate_along_characteristic(Bubble, Current, hinf))
    {
      // Pressure has already been evaluated iteratively while integrating along the characteristic

      // Move to the next node
      Current = Current->backward;
    }
    else  // Physically implausible enthalpy/pressure detected, node is discarded. This generally only occurs with the explicit velocity (EV) integration.
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

  // ---------------------------------------
  // Shock treatment (if necessary)

  int check_for_shocks = 1;
  while (check_for_shocks)
  {
    Current = Bubble->Emissions->FirstNode;
    check_for_shocks = 0;

    while (Current != NULL)
    {
      if (Current->forward != NULL && Current->r > Current->forward->r)  // Check for shock formation
      {
        check_for_shocks = 1;  // Requires another complete subsequent iteration to make sure all multivalued solutions have been or are treated

        // Define new properties of the forward node
        Current->forward->r = 0.5 * (Current->r + Current->forward->r);
        Current->forward->h = 0.5 * (Current->h + Current->forward->h);
        Current->forward->u = 0.5 * (Current->u + Current->forward->u);
        Current->forward->g = Current->forward->r * (Current->forward->h - hinf + 0.5 * APECSS_POW2(Current->forward->u));

        // Current node is obsolete and will be deleted
        struct APECSS_EmissionNode *Obsolete = Current;

        // Redefine the neighbor information
        Current->forward->backward = Current->backward;
        if (Current->backward != NULL)
          Current->backward->forward = Current->forward;
        else
          Bubble->Emissions->FirstNode = Current->forward;

        // Move to the next node
        Current = Current->forward;

        // Delete obsolete node
        free(Obsolete);
        Bubble->Emissions->nNodes -= 1;
      }
      else
      {
        // Move to the next node without action
        Current = Current->forward;
      }
    }
  }

  // ---------------------------------------
  // Compute the pressure and (if applicable) the invariant f

  Current = Bubble->Emissions->FirstNode;
  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  while (Current != NULL)
  {
    APECSS_FLOAT p_prev;
    do  // Compute pressure iteratively
    {
      p_prev = Current->p;
      APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
      Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
    } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

    Current->f = Bubble->Emissions->compute_f(Bubble, Current->g, Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid));

    Current = Current->forward;
  }

  // ---------------------------------------
  // Store data (if applicable)

  if (Bubble->Results->Emissions->nNodes + Bubble->Results->Emissions->MinMaxPeriod > 0)  // Results of at least one specific node are to be stored.
  {
    Current = Bubble->Emissions->FirstNode;

    while (Current != NULL)
    {
      Bubble->results_emissionsnode_store(
          Current, Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid), pinf, Bubble);
      Current = Current->forward;
    }
  }

  return (0);
}

// -------------------------------------------------------------------
// INTEGRATION OF POSITION AND VELOCITY ALONG CHARACTERISTIC
// -------------------------------------------------------------------
// Functions integrating the position and velocity of the emission
// nodes dependent on the type of velocity integration and scheme.
// -------------------------------------------------------------------
// The appropriate function is chosen in apecss_bubble_processoptions()
// and hooked up to the function pointer
// Bubble->Emissions->integrate_along_characteristic().
// -------------------------------------------------------------------

int apecss_emissions_integrate_ev_tait_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Simplified method based on the Kirkwood-Bethe hypothesis, with u treated explicitly.

  APECSS_FLOAT c = APECSS_SQRT((Bubble->Liquid->Gamma - 1.0) * Current->h);
  APECSS_FLOAT u = Current->u;

  Current->r += Bubble->dt * (c + u);
  Current->u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(Current->r) * Current->r) +
               Current->g / (Bubble->get_dimensionalradius(Current->r) * (c + u));
  Current->h = hinf + Current->g / Bubble->get_dimensionalradius(Current->r) - 0.5 * APECSS_POW2(Current->u);

  if (Current->h < 0)  // Flag a physically implausible solution
    return (0);
  else
    return (1);
}

int apecss_emissions_integrate_ev_tait_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Simplified method based on the Kirkwood-Bethe hypothesis, with u treated explicitly.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT dt = Bubble->dt;

  APECSS_FLOAT c_old = APECSS_SQRT((Gamma - 1.0) * Current->h);

  // Step 1

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c_old + u;

  APECSS_FLOAT h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  if (h < 0) return (0);  // Flag a physically implausible solution

  APECSS_FLOAT c = APECSS_SQRT((Gamma - 1.0) * h);

  // Step 2

  r = Current->r + 0.5 * dt * kr1;
  u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(r) * r) + Current->g / (Bubble->get_dimensionalradius(r) * (c + u));
  APECSS_FLOAT kr2 = c + u;

  h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  if (h < 0) return (0);  // Flag a physically implausible solution

  c = APECSS_SQRT((Gamma - 1.0) * h);

  // Step 3

  r = Current->r + 0.5 * dt * kr2;
  u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(r) * r) + Current->g / (Bubble->get_dimensionalradius(r) * (c + u));
  APECSS_FLOAT kr3 = c + u;

  h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  if (h < 0) return (0);  // Flag a physically implausible solution

  c = APECSS_SQRT((Gamma - 1.0) * h);

  // Step 4

  r = Current->r + dt * kr3;
  u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(r) * r) + Current->g / (Bubble->get_dimensionalradius(r) * (c + u));
  APECSS_FLOAT kr4 = c + u;

  // Solve

  Current->r += dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(Current->r) * Current->r) +
               Current->g / (Bubble->get_dimensionalradius(Current->r) * (c_old + Current->u));
  Current->h = hinf + Current->g / Bubble->get_dimensionalradius(Current->r) - 0.5 * APECSS_POW2(Current->u);

  if (Current->h < 0)  // Flag a physically implausible solution
    return (0);
  else
    return (1);
}

int apecss_emissions_integrate_tiv_tait_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Method of Hickling & Plesset (1963), see Appendix I, and Ebeling (1978), with integration of u over t.

  APECSS_FLOAT c = APECSS_SQRT((Bubble->Liquid->Gamma - 1.0) * Current->h);
  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;

  Current->r += Bubble->dt * (c + u);
  Current->u += Bubble->dt * (Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u)));
  Current->h = hinf + Current->g / Bubble->get_dimensionalradius(Current->r) - 0.5 * APECSS_POW2(Current->u);

  return (1);
}

int apecss_emissions_integrate_tiv_tait_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Method of Hickling & Plesset (1963), see Appendix I, and Ebeling (1978), with integration of u over t.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT dt = Bubble->dt;

  APECSS_FLOAT c = APECSS_SQRT((Gamma - 1.0) * Current->h);

  // Step 1

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c + u;
  APECSS_FLOAT ku1 = Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u)));

  // Step 2

  r = Current->r + 0.5 * dt * kr1;
  u = Current->u + 0.5 * dt * ku1;
  APECSS_FLOAT kr2 = c + u;
  APECSS_FLOAT ku2 = Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u)));

  // Step 3

  r = Current->r + 0.5 * dt * kr2;
  u = Current->u + 0.5 * dt * ku2;
  APECSS_FLOAT kr3 = c + u;
  APECSS_FLOAT ku3 = Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u));

  c = APECSS_SQRT((Gamma - 1.0) * (hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u)));

  // Step 4

  r = Current->r + dt * kr3;
  u = Current->u + dt * ku3;
  APECSS_FLOAT kr4 = c + u;
  APECSS_FLOAT ku4 = Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u));

  // Solve

  Current->r += dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->u += dt * (ku1 + 2.0 * (ku2 + ku3) + ku4) * APECSS_ONESIXTH;
  Current->h = hinf + Current->g / Bubble->get_dimensionalradius(Current->r) - 0.5 * APECSS_POW2(Current->u);

  return (1);
}

int apecss_emissions_integrate_ev_general_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Simplified method based on the Kirkwood-Bethe hypothesis, with u treated explicitly.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);
  APECSS_FLOAT u = Current->u;

  Current->r += Bubble->dt * (c + u);
  Current->u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(Current->r) * Current->r) +
               Current->g / (Bubble->get_dimensionalradius(Current->r) * (c + Current->u));
  Current->h = hinf + Current->g / Bubble->get_dimensionalradius(Current->r) - 0.5 * APECSS_POW2(Current->u);

  APECSS_FLOAT p_prev;
  do
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (Current->p < -B) return (0);  // // Flag a physically implausible solution
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

int apecss_emissions_integrate_ev_general_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Simplified method based on the Kirkwood-Bethe hypothesis, with u treated explicitly.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT dt = Bubble->dt;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c_old = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);

  // Step 1

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c_old + u;

  APECSS_FLOAT h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  APECSS_FLOAT p = Current->p, p_prev;
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (p < -B) return (0);  // Flag a physically implausible solution
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  // Step 2

  r = Current->r + 0.5 * dt * kr1;
  u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(r) * r) + Current->g / (Bubble->get_dimensionalradius(r) * (c + u));
  APECSS_FLOAT kr2 = c + u;

  h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (p < -B) return (0);  // Flag a physically implausible solution
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));

  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  // Step 3

  r = Current->r + 0.5 * dt * kr2;
  u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(r) * r) + Current->g / (Bubble->get_dimensionalradius(r) * (c + u));
  APECSS_FLOAT kr3 = c + u;

  h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (p < -B) return (0);  // Flag a physically implausible solution
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));

  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  // Step 4

  r = Current->r + dt * kr3;
  u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(r) * r) + Current->g / (Bubble->get_dimensionalradius(r) * (c + u));
  APECSS_FLOAT kr4 = c + u;

  // Solve

  Current->r += dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->u = Bubble->dimensionality * Current->f / (2.0 * Bubble->get_dimensionalradius(Current->r) * Current->r) +
               Current->g / (Bubble->get_dimensionalradius(Current->r) * (c_old + Current->u));
  Current->h = hinf + Current->g / Bubble->get_dimensionalradius(Current->r) - 0.5 * APECSS_POW2(Current->u);

  do  // Compute pressure iteratively
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);

    if (Current->p < -B) return (0);  // Flag a physically implausible solution
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

int apecss_emissions_integrate_tiv_general_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
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
  Current->u += Bubble->dt * (Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u)));
  Current->h = hinf + Current->g / Bubble->get_dimensionalradius(Current->r) - 0.5 * APECSS_POW2(Current->u);

  APECSS_FLOAT p_prev;
  do  // Compute pressure iteratively
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

int apecss_emissions_integrate_tiv_general_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf)
{
  // Method of Hickling & Plesset (1963), see Appendix I, and Ebeling (1978), with integration of u over t.

  APECSS_FLOAT Gamma = Bubble->Liquid->Gamma;
  APECSS_FLOAT B = Bubble->Liquid->B;
  APECSS_FLOAT b = Bubble->Liquid->b;
  APECSS_FLOAT dt = Bubble->dt;
  APECSS_FLOAT tol = Bubble->Emissions->KB_IterTolerance;

  APECSS_FLOAT c = Bubble->Liquid->get_soundspeed(Current->p, Bubble->Liquid->get_density(Current->p, Bubble->Liquid), Bubble->Liquid);

  // Step 1

  APECSS_FLOAT r = Current->r;
  APECSS_FLOAT u = Current->u;
  APECSS_FLOAT kr1 = c + u;
  APECSS_FLOAT ku1 = Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u));

  APECSS_FLOAT h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  APECSS_FLOAT p = Current->p, p_prev;
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  // Step 2

  r = Current->r + 0.5 * dt * kr1;
  u = Current->u + 0.5 * dt * ku1;
  APECSS_FLOAT kr2 = c + u;
  APECSS_FLOAT ku2 = Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u));

  h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  // Step 3

  r = Current->r + 0.5 * dt * kr2;
  u = Current->u + 0.5 * dt * ku2;
  APECSS_FLOAT kr3 = c + u;
  APECSS_FLOAT ku3 = Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u));

  h = hinf + Current->g / Bubble->get_dimensionalradius(r) - 0.5 * APECSS_POW2(u);
  do  // Compute pressure iteratively
  {
    p_prev = p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(p, Bubble->Liquid);
    p = ((Gamma - 1.0) * rho * h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - p)) > tol * APECSS_ABS(p));
  c = Bubble->Liquid->get_soundspeed(p, Bubble->Liquid->get_density(p, Bubble->Liquid), Bubble->Liquid);

  // Step 4

  r = Current->r + dt * kr3;
  u = Current->u + dt * ku3;
  APECSS_FLOAT kr4 = c + u;
  APECSS_FLOAT ku4 = Bubble->dimensionality * ((c + u) * Current->g / (2.0 * Bubble->get_dimensionalradius(r)) - APECSS_POW2(c) * u) / (r * (c - u));

  // Solve

  Current->r += dt * (kr1 + 2.0 * (kr2 + kr3) + kr4) * APECSS_ONESIXTH;
  Current->u += dt * (ku1 + 2.0 * (ku2 + ku3) + ku4) * APECSS_ONESIXTH;
  Current->h = hinf + Current->g / Bubble->get_dimensionalradius(Current->r) - 0.5 * APECSS_POW2(Current->u);

  do  // Compute pressure iteratively
  {
    p_prev = Current->p;
    APECSS_FLOAT rho = Bubble->Liquid->get_density(Current->p, Bubble->Liquid);
    Current->p = ((Gamma - 1.0) * rho * Current->h - (1.0 - b * rho) * Gamma * B) / (Gamma - b * rho);
  } while (APECSS_ABS((p_prev - Current->p)) > tol * APECSS_ABS(Current->p));

  return (1);
}

// -------------------------------------------------------------------
// INVARIANT f
// -------------------------------------------------------------------
// Functions computing the invariant f dependent on the model chosen
// for the emissions.
// -------------------------------------------------------------------
// The appropriate function is chosen in apecss_bubble_processoptions()
// and hooked up to the function pointer
// Bubble->Emissions->compute_f().
// -------------------------------------------------------------------

APECSS_FLOAT apecss_emissions_f_zero(struct APECSS_Bubble *Bubble, APECSS_FLOAT g, APECSS_FLOAT pL, APECSS_FLOAT rhoL) { return (0.0); }

APECSS_FLOAT apecss_emissions_f_finitespeedincompressible(struct APECSS_Bubble *Bubble, APECSS_FLOAT g, APECSS_FLOAT pL, APECSS_FLOAT rhoL)
{
  return (APECSS_POW2(Bubble->R) * Bubble->U);
}

APECSS_FLOAT apecss_emissions_f_quasiacoustic(struct APECSS_Bubble *Bubble, APECSS_FLOAT g, APECSS_FLOAT pL, APECSS_FLOAT rhoL)
{
  return (APECSS_POW2(Bubble->R) * Bubble->U - Bubble->R * g / Bubble->Liquid->cref);
}

APECSS_FLOAT apecss_emissions_f_kirkwoodbethe(struct APECSS_Bubble *Bubble, APECSS_FLOAT g, APECSS_FLOAT pL, APECSS_FLOAT rhoL)
{
  return (2.0 *
          (Bubble->get_dimensionalradius(Bubble->R) * Bubble->R * Bubble->U -
           Bubble->R * g / (Bubble->Liquid->get_soundspeed(pL, rhoL, Bubble->Liquid) + Bubble->U)) /
          (Bubble->dimensionality + APECSS_SMALL));
}