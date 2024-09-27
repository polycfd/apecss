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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifndef APECSS_H_
#define APECSS_H_

#define APECSS_VERSION_NUM (1.7)
static const char APECSS_RELEASE_DATE[] = "27-Sep-2024";

// -------------------------------------------------------------------
// CONSTANTS & MACROS
// -------------------------------------------------------------------

// Floating point precision (default is double precision)
// #define APECSS_PRECISION_LONGDOUBLE

// Functions dependent on chosen floating point precision
#if defined(APECSS_PRECISION_LONGDOUBLE)
typedef long double APECSS_FLOAT;
#define APECSS_EPS (1.0e-18)
#define APECSS_SMALL (1.0e-40)
#define APECSS_POW(a, b) (powl(a, b))
#define APECSS_SQRT(a) (sqrtl(a))
#define APECSS_ABS(a) (fabsl(a))
#define APECSS_SIN(a) (sinl(a))
#define APECSS_COS(a) (cosl(a))
#define APECSS_EXP(a) (expl(a))
#define APECSS_LOG(a) (logl(a))
#define APECSS_CEIL(a) (ceill(a))
#define APECSS_STRINGTOFLOAT(a) (strtold(a, NULL))
#else
typedef double APECSS_FLOAT;
#define APECSS_EPS (1.0e-15)
#define APECSS_SMALL (1.0e-30)
#define APECSS_POW(a, b) (pow(a, b))
#define APECSS_SQRT(a) (sqrt(a))
#define APECSS_ABS(a) (fabs(a))
#define APECSS_SIN(a) (sin(a))
#define APECSS_COS(a) (cos(a))
#define APECSS_EXP(a) (exp(a))
#define APECSS_LOG(a) (log(a))
#define APECSS_CEIL(a) (ceil(a))
#define APECSS_STRINGTOFLOAT(a) (strtod(a, NULL))
#endif

// Simple operations
#define APECSS_POW2(x) ((x) * (x))
#define APECSS_POW3(x) ((x) * (x) * (x))
#define APECSS_POW4(x) ((x) * (x) * (x) * (x))
#define APECSS_MAX(a, b) ((a) > (b) ? (a) : (b))
#define APECSS_MAX3(a, b, c) (((a) >= (b) && (a) >= (c)) ? (a) : ((b) >= (a) && (b) >= (c)) ? (b) : (c))
#define APECSS_MIN(a, b) ((a) < (b) ? (a) : (b))
#define APECSS_MIN3(a, b, c) (((a) <= (b) && (a) <= (c)) ? (a) : ((b) <= (a) && (b) <= (c)) ? (b) : (c))

// Pre-defined constants
#define APECSS_PI (3.1415926535897932384626433832795028841972)
#define APECSS_E (2.7182818284590452353602874713526624977572)
#define APECSS_ONETHIRD (0.3333333333333333333333333333333333333333)
#define APECSS_ONESIXTH (0.1666666666666666666666666666666666666667)
#define APECSS_AVOGADRO (6.02214076e23)
#define APECSS_LN_OF_2 (0.693147180559945309)
#define APECSS_LN_OF_10 (2.302585092994045684)
#define APECSS_LARGE (1.0e15)

// Runge-Kutta scheme
#define APECSS_RK54_7M (0)  // RK5(4)7M (minimum truncation) coefficients of Dormand & Prince (1980)
#define APECSS_RK54_7S (1)  // RK5(4)7S (stability optimized) coefficients of Dormand & Prince (1980)

// Bubble model
#define APECSS_BUBBLEMODEL_RP (1)  // Standard Rayleigh-Plesset model
#define APECSS_BUBBLEMODEL_RP_ACOUSTICRADIATION (2)  // Rayleigh-Plesset model incl. acoustic radiation term
#define APECSS_BUBBLEMODEL_KELLERMIKSIS (3)  // Keller-Miksis model
#define APECSS_BUBBLEMODEL_GILMORE (4)  // Gilmore model

// Gas equations of state
#define APECSS_GAS_IG (1)  // Ideal gas EoS
#define APECSS_GAS_HC (2)  // Hard-core gas
#define APECSS_GAS_NASG (3)  // Noble-Abel-stiffend-gas EoS

// Liquid equations of state
#define APECSS_LIQUID_TAIT (1)  // Tait EoS
#define APECSS_LIQUID_NASG (2)  // Noble-Abel-stiffend-gas EoS

// Liquid models
#define APECSS_LIQUID_NEWTONIAN (1)  // Newtonian fluid
#define APECSS_LIQUID_KELVINVOIGT (2)  // Kelvin-Voigt solid
#define APECSS_LIQUID_ZENER (3)  // Zener solid, requires two additional ODEs
#define APECSS_LIQUID_OLDROYDB (4)  // Oldroyd-B liquid, requires two additional ODEs

// Lipid coating model (bit-wise)
#define APECSS_LIPIDCOATING_NONE (1)  // Clean interface
#define APECSS_LIPIDCOATING_MARMOTTANT (2)  // Lipid-coating model of Marmottant et al. (2005)
#define APECSS_LIPIDCOATING_GOMPERTZFUNCTION (4)  // Lipid-coating model of Guemmer et al. (2021)

// Excitation types
#define APECSS_EXCITATION_NONE (0)  // No excitation
#define APECSS_EXCITATION_SIN (1)  // Sinusoidal excitation

// Emission type (bit-wise)
#define APECSS_EMISSION_NONE (0)  // No emissions are computed | 0000 0000
#define APECSS_EMISSION_INCOMPRESSIBLE (1)  // Incompressible emissions, see Neppiras (1980) | 0000 0001
#define APECSS_EMISSION_FINITE_SPEED_INCOMPRESSIBLE (2)  // Incompressible emissions tracked with a finite speed of sound | 0000 0010
#define APECSS_EMISSION_QUASIACOUSTIC (4)  // Quasi-acoustic model of Trilling/Gilmore (1952) | 0000 0100
#define APECSS_EMISSION_KIRKWOODBETHE (16)  // A model based on the Kirkwood-Bethe hypothesis (EKB, GFC, HPE) is used | 0001 0000
#define APECSS_EMISSION_EV (17)  // Explicit expression for velocity of Denner & Schenke | 0001 0001
#define APECSS_EMISSION_TIV (18)  // Temporally-intergrated velocity of Hickling & Plesset (1963) | 0001 0010

// Scheme to integrate emissions along outgoing characteristic
#define APECSS_EMISSION_INTEGRATE_EULER (0)  // Euler scheme
#define APECSS_EMISSION_INTEGRATE_RK4 (1)  // Conventional fourth-order Runge-Kutta scheme

// Results
#define APECSS_RESULTS_DISCARD (0)  // Discard the results and do not write results to file.
#define APECSS_RESULTS_WRITE (1)  // Write all data to file in one go.
#define APECSS_RESULTS_APPEND (2)  // Append the results file with a new batch of results.

// Misc
#define APECSS_DATA_ALLOC_INCREMENT (10000)  // Allocation increment for the arrays in which the results are stored
#define APECSS_STRINGLENGTH (512)  // Standard string length
#define APECSS_STRINGLENGTH_SPRINTF (1024)  // String length for printing in the terminal
#define APECSS_STRINGLENGTH_SPRINTF_LONG (2048)  // Length of a long string for printing in the terminal

// Debugging aids, returning the line number and file name in which the macro is called
#define APECSS_WHERE printf("HERE - %s:%d\n", __FILE__, __LINE__);
#define APECSS_WHERE_INT(a) printf("HERE - %s:%d = %i\n", __FILE__, __LINE__, a);
#if defined(APECSS_PRECISION_LONGDOUBLE)
#define APECSS_WHERE_FLOAT(a) printf("HERE - %s:%d = %Le\n", __FILE__, __LINE__, a);
#else
#define APECSS_WHERE_FLOAT(a) printf("HERE - %s:%d = %e\n", __FILE__, __LINE__, a);
#endif

// -------------------------------------------------------------------
// VARIABLE DECLARATIONS
// -------------------------------------------------------------------

struct APECSS_Bubble;  // Dummy for the function pointer below

struct APECSS_Gas
{
  int EoS;  // Equation of state

  APECSS_FLOAT Gamma;  // Polytropic exponent
  APECSS_FLOAT B;  // Attractive pressure constant [Pa]
  APECSS_FLOAT b;  // Co-volume [m^3/kg]
  APECSS_FLOAT dmol;  // Molecular kinetic diameter [m]
  APECSS_FLOAT mmol;  // Molecular weight [kg/mol]
  APECSS_FLOAT pref;  // Reference pressure [Pa]
  APECSS_FLOAT rhoref;  // Reference density [kg/m^3]
  APECSS_FLOAT Kref;  // Constant reference coefficient of a polytropic EoS

  // Pointers to the functions describing the gas pressure and its derivative
  APECSS_FLOAT (*get_pressure)(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative)(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

  // Void pointer for additional user-defined data
  void *user_data;
};

struct APECSS_Liquid
{
  int Type;  // Type of liquid (i.e. Newtonian or viscoelastic)
  int EoS;  // Equation of state

  APECSS_FLOAT Gamma;  // Polytropic exponent
  APECSS_FLOAT B;  // Attractive pressure constant [Pa]
  APECSS_FLOAT b;  // Co-volume [m^3/kg]
  APECSS_FLOAT h;  // van der Waals hardcore radius [m]
  APECSS_FLOAT dmol;  // Molecular diameter [m]
  APECSS_FLOAT mmol;  // Molecular weight [kg/mol]
  APECSS_FLOAT mu;  // Viscosity [Pa s]
  APECSS_FLOAT pref;  // Reference pressure [Pa]
  APECSS_FLOAT rhoref;  // Reference density [kg/m^3]
  APECSS_FLOAT cref;  // Reference speed of sound [m/s^2]
  APECSS_FLOAT Kref;  // Constant reference coefficient of a polytropic EoS
  APECSS_FLOAT G;  // Shear modulus [Pa]
  APECSS_FLOAT eta;  // Polymer viscosity [Pa s]
  APECSS_FLOAT lambda;  // Relaxation time [s]

  // Pointers to the functions describing the properties of the liquid
  APECSS_FLOAT (*get_density)(APECSS_FLOAT p, struct APECSS_Liquid *Liquid);
  APECSS_FLOAT (*get_soundspeed)(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
  APECSS_FLOAT (*get_enthalpy)(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);

  // Pointers to the functions describing the liquid pressure and its derivative at the bubble wall
  APECSS_FLOAT (*get_pressure_bubblewall)(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_bubblewall_expl)(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

  // Pointers to the functions describing the pressure and its derivative due to viscous contributions
  APECSS_FLOAT (*get_pressure_viscous)(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_viscous_expl)(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_viscous_impl)(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);

  // Void pointer for additional user-defined data
  void *user_data;
};

struct APECSS_Interface
{
  APECSS_FLOAT sigma;  // Surface tension coefficient [N/m]

  // Lipid coating properties of the Marmottant and Marmottant-Gompertz models
  int LipidCoatingModel;  // Coating model for the lipid coating of the interface
  APECSS_FLOAT Elasticity;  // Surface elasticity
  APECSS_FLOAT Viscosity;  // Dilatational viscosity
  APECSS_FLOAT sigma0;  // Initial surface tension coefficient [N/m]
  APECSS_FLOAT GompertzC;  // C coefficient of the Marmottant-Gompertz model

  // Pointer to the function defining the surface tension coefficient
  APECSS_FLOAT (*get_surfacetension)(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);

  // Pointers to the functions describing the pressure and its derivative due to surface tension
  APECSS_FLOAT (*get_pressure_surfacetension)(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_surfacetension)(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);

  // Pointers to the functions describing the pressure and its derivative due to surface viscous contributions
  APECSS_FLOAT (*get_pressure_viscous)(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
  APECSS_FLOAT (*get_pressurederivative_viscous_expl)(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
  APECSS_FLOAT (*get_pressurederivative_viscous_impl)(APECSS_FLOAT R, struct APECSS_Interface *Interface);

  // Void pointer for additional user-defined data
  void *user_data;
};

struct APECSS_EmissionNode
{
  int id;  // Unique ID of the emission node
  APECSS_FLOAT r;  // Radial coordinate
  APECSS_FLOAT h;  // Enthalpy
  APECSS_FLOAT p;  // Pressure
  APECSS_FLOAT u;  // Velocity
  APECSS_FLOAT f;  // = const.
  APECSS_FLOAT g;  // = const.
  struct APECSS_EmissionNode *forward;  // Forward (i.e. outward) neighbor
  struct APECSS_EmissionNode *backward;  // Backward (i.e. inward) neighbor

  // Void pointer for additional user-defined data
  void *user_data;
};

struct APECSS_Emissions
{
  int Type;  // Model for the acoustic emissions
  int Scheme;  // Scheme used to integrate along the outgoing characteristic
  APECSS_FLOAT CutOffDistance;  // Distance above which the nodes are deleted
  APECSS_FLOAT KB_IterTolerance;  // Iteration tolerance to obtain pressure in the general Kirkwood-Bethe model

  int nNodes;  // Total number of emission nodes in the linked list
  struct APECSS_EmissionNode *FirstNode;  // First node of the linked list
  struct APECSS_EmissionNode *LastNode;  // Last node of the linked list

  // Pruning the list of emission nodes
  int pruneList;  // Flag indicating whether the linked list is pruned
  int (*prune_test)(struct APECSS_EmissionNode *Node);

  // Pointer to the function advancing the emission nodes
  int (*advance)(struct APECSS_Bubble *Bubble);

  // Pointer to the function that computes the invariant f
  APECSS_FLOAT (*compute_f)(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Node);

  // Pointer to the function integrating the radial position and velocity along the outgoing characteristic
  int (*integrate_along_characteristic)(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
};

struct APECSS_NumericsODE
{
  int RKtype;  // Runge-Kutta scheme

  APECSS_FLOAT dtMin;  // Minimum allowable time-step
  APECSS_FLOAT dtMax;  // Maximum allowable time-step
  APECSS_FLOAT minScale;  // Minimum allowable value by which the dt-controller multiplies dt
  APECSS_FLOAT maxScale;  // Maximum allowable value by which the dt-controller multiplies dt
  APECSS_FLOAT tol;  // Desired tolerance of the solution
  APECSS_FLOAT control_coeff_alpha;
  APECSS_FLOAT control_coeff_q;
  int maxSubIter;  // Maximum number of sub-iterations to control the error

  // Coefficients of the Runge-Kutta method
  APECSS_FLOAT a21;
  APECSS_FLOAT a31, a32;
  APECSS_FLOAT a41, a42, a43;
  APECSS_FLOAT a51, a52, a53, a54;
  APECSS_FLOAT a61, a62, a63, a64, a65;
  APECSS_FLOAT a71, a72, a73, a74, a75, a76;
  APECSS_FLOAT b1, b3, b4, b5, b6, b7;
  APECSS_FLOAT bs1, bs3, bs4, bs5, bs6, bs7;
  APECSS_FLOAT c2, c3, c4, c5, c6, c7;
  APECSS_FLOAT e1, e3, e4, e5, e6, e7;
};

struct APECSS_Excitation
{
  int type;  // Type of the excitation
  APECSS_FLOAT f;  // Frequency [Hz]
  APECSS_FLOAT dp;  // Maximum pressure amplitude [Pa]
};

struct APECSS_Interaction
{
  int nBubbles;  // Number of bubbles in the whole cluster considered
  APECSS_FLOAT location[3];  // 3D coordinates of the considered bubble
  APECSS_FLOAT dp_neighbor;  // Pressure induced by interactions with neighboring bubbles
  APECSS_FLOAT last_t_1;  // Previous timestep when interaction where considered
  APECSS_FLOAT last_t_2;  // Previous timestep when interaction where considered before last_t_1
  APECSS_FLOAT last_p_1;  // Pressure induced by neighboring bubbles or total pressure in the far field at last_t_1
  APECSS_FLOAT last_p_2;  // Pressure induced by neighboring bubbles or total pressure in the far field at last_t_2
};

struct APECSS_ResultsBubble
{
  int freq;  // Frequency with which the results are stored (with respect to the time-step number)
  int n;  // Number of stored results
  int nAllocated;  // Length of the arrays

  APECSS_FLOAT *t;  // Time [s]
  APECSS_FLOAT *R;  // Bubble radius [m]
  APECSS_FLOAT *U;  // Bubble wall velocity [m/s]
  APECSS_FLOAT *dt;  // Time-step size [s]
  APECSS_FLOAT *pG;  // Pressure of the gas inside the bubble [Pa]
  APECSS_FLOAT *pL;  // Pressure in the liquid at the bubble wall [Pa]
  APECSS_FLOAT *pinf;  // Pressure in the liquid at the bubble wall [Pa]
  APECSS_FLOAT *cL;  // Speed of sound in the liquid at the bubble wall [m/s]

  int nUserODEs;  // The solutions of the first nUserODEs additional user-defined ODEs are written to file
  APECSS_FLOAT **UserODEsSol;  // Solution of the user-defined ODEs to be written to file
  char **UserODEsName;  // Name of the value of the user-defined ODEs to be written to file (e.g. 'T' for temperature)
};

struct APECSS_ResultsEmissionsSpace
{
  int n;  // Number of stored results
  int nAllocated;  // Length of the arrays
  APECSS_FLOAT RadialLocation;  // Radial location at which the results are recorded [m]
  APECSS_FLOAT *t;  // Time [s]
  APECSS_FLOAT *p;  // Absolute pressure [Pa]
  APECSS_FLOAT *u;  // Velocity [m/s]
  APECSS_FLOAT *c;  // Speed of sound [m/s]
  APECSS_FLOAT *pInf;  // Pressure at infinity at the given time [Pa]
};

struct APECSS_ResultsEmissionsNode
{
  int id;  // Desired ID of the emission node that it to be written out
  int n;  // Number of stored results
  int nAllocated;  // Length of the arrays
  int *real_id;  // Actual ID of the emission node
  APECSS_FLOAT *r;  // Radial location [m]
  APECSS_FLOAT *t;  // Time [s]
  APECSS_FLOAT *p;  // Absolute pressure [Pa]
  APECSS_FLOAT *u;  // Velocity [m/s]
  APECSS_FLOAT *c;  // Speed of sound [m/s]
  APECSS_FLOAT *pInf;  // Pressure at infinity at the given time [Pa]
};

struct APECSS_ResultsEmissions
{
  // Results of the acoustic emissions at specific time instances (option: EmissionsTime)
  int nTimeInstances;  // Number of time instances at which the acoustic emissions ought to be written out
  int nextTimeInstance;  // Count of the next time instance to be written out
  APECSS_FLOAT *TimeInstances;  // Time instances that are to be written out

  // Results of the acoustic emissions recorded at specific radial locations (option: EmissionsSpace)
  int freqSpaceLocations;  // Frequency with which the results are stored (with respect to the time-step number)
  int nSpaceLocations;  // Number of radial locations
  struct APECSS_ResultsEmissionsSpace *SpaceLocation;  // Array of structures containing the results at the specified radial locations

  // Results of a specific emission node (option: EmissionsNode)
  int nNodes;  // Number of nodes
  struct APECSS_ResultsEmissionsNode *Node;  // Array of structures containing the results for the specified nodes

  // Results of the acoustic emission at min/max instances (option EmissionsMinMax)
  int MinMaxPeriod;  // Excitation period in which the min/max instances are to be written out
  struct APECSS_ResultsEmissionsNode *Node_Rmin, *Node_Umin, *Node_pLmax;  // Array of structures containing the results for the specified nodes
  APECSS_FLOAT Rmin, Umin, pLmax;  // Respective min/max values

  // Pointer to the function allocating the solution arrays
  int (*allocation_nodespecific)(struct APECSS_Bubble *Bubble);
  int (*allocation_nodeminmax)(struct APECSS_Bubble *Bubble);

  // Pointer to the function storing the solution
  int (*store_nodespecific)(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
  int (*store_nodeminmax)(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
};

struct APECSS_Results
{
  int digits;  // Number of digits for the output
  char dir[APECSS_STRINGLENGTH_SPRINTF];  // Directory path where the results are stored

  struct APECSS_ResultsBubble *RayleighPlesset;  // Structure containing the results of the Rayleigh-Plesset model (if applicable)
  struct APECSS_ResultsEmissions *Emissions;  // Structure containing the results of the acoustic emissions (if applicable)
};

struct APECSS_Bubble
{
  APECSS_FLOAT tStart;  // Start time [s]
  APECSS_FLOAT tEnd;  // End time [s]

  APECSS_FLOAT dt;  // Time-step

  int RPModel;  // Model governing the bubble dynamics
  APECSS_FLOAT dimensionality;  // Dimensionality of the bubble

  // Primary variables
  APECSS_FLOAT t;  // Time [s]
  APECSS_FLOAT R;  // Radius of the bubble [m]
  APECSS_FLOAT U;  // Velocity of the bubble wall [m/s]

  // Ambient and initial conditions
  APECSS_FLOAT R0;  // Initial radius [m]
  APECSS_FLOAT U0;  // Initial velocity [m/s]
  APECSS_FLOAT p0;  // Ambient pressure [Pa]
  APECSS_FLOAT pG0;  // Initial gas pressure [Pa]
  APECSS_FLOAT rhoG0;  // Initial gas density [kg/m^3]
  APECSS_FLOAT r_hc;  // van der Waals hardcore radius [m]

  // Interface properties associated with a lipid monolayer coating
  APECSS_FLOAT Rbuck;  // Buckling radius of the Marmottant/Marmottant-Gompertz model [m]
  APECSS_FLOAT Rrupt;  // Rupture radius of the Marmottant/Marmottant-Gompertz model [m]
  APECSS_FLOAT GompertzB;  // B coefficient of the Marmottant-Gompertz model

  // ODEs
  int nODEs;  // Total number of ODEs
  int nUserODEs;  // Number of additional user-defined ODEs
  int dtNumber;  // Time-step number
  int nSubIter;  // Total number of sub-iterations to control the error
  APECSS_FLOAT *ODEsSol;  // Solution of each ODE
  APECSS_FLOAT *ODEsSolOld;  // // Old solution of each ODE, required for sub-iterations
  APECSS_FLOAT err;  // Solution error
  APECSS_FLOAT (**ode)(APECSS_FLOAT *, APECSS_FLOAT, struct APECSS_Bubble *);  // Array of pointers to the functions containing the ODEs
  APECSS_FLOAT *k2, *k3, *k4, *k5, *k6, *k7, *kLast;  // Intermediate solutions of the Runge-Kutta solver
  struct APECSS_NumericsODE *NumericsODE;  // Structure containg the parameters of the Runge-Kutta solver

  // Properties of the fluids and interface
  struct APECSS_Gas *Gas;
  struct APECSS_Liquid *Liquid;
  struct APECSS_Interface *Interface;

  // Excitation of the bubble (if applicable)
  struct APECSS_Excitation *Excitation;

  // Interactions with neighboring bubbles (if applicable)
  struct APECSS_Interaction *Interaction;

  // Pointers to the functions describing the liquid pressure and its derivative at infinity
  APECSS_FLOAT (*get_pressure_infinity)(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_infinity)(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

  // Pointer to the function returning the radius r^{\alpha/2}, where \alpha defines the dimensionality of the bubble.
  APECSS_FLOAT (*get_dimensionalradius)(APECSS_FLOAT r);

  // Lagrangian emissions (if applicable)
  struct APECSS_Emissions *Emissions;
  int (*emissions_initialize)(struct APECSS_Bubble *Bubble);
  int (*emissions_update)(struct APECSS_Bubble *Bubble);
  int (*emissions_free)(struct APECSS_Bubble *Bubble);

  // Results (if applicable)
  struct APECSS_Results *Results;
  int (*results_rayleighplesset_store)(struct APECSS_Bubble *Bubble);
  int (*results_emissionstime_write)(struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*results_emissionstime_check)(APECSS_FLOAT tend, struct APECSS_Bubble *Bubble);
  int (*results_emissionsspace_store)(struct APECSS_Bubble *Bubble);
  int (*results_emissionsnodeminmax_identify)(struct APECSS_Bubble *Bubble);
  int (*results_emissionsnode_alloc)(struct APECSS_Bubble *Bubble);
  int (*results_emissionsnode_store)(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);

  // Progress screen (if applicable)
  int progress;
  int (*progress_initial)();
  int (*progress_update)(int *prog, APECSS_FLOAT t, APECSS_FLOAT totaltime);
  int (*progress_final)();

  // Void pointer for additional user-defined data
  void *user_data;
};

// -------------------------------------------------------------------
// FUNCTION DECLARATIONS
// -------------------------------------------------------------------

// ---------------------
// bubble.c

int apecss_bubble_initializestruct(struct APECSS_Bubble *Bubble);
int apecss_bubble_setdefaultoptions(struct APECSS_Bubble *Bubble);
int apecss_bubble_readoptions(struct APECSS_Bubble *Bubble, char *OptionsDir);
int apecss_bubble_processoptions(struct APECSS_Bubble *Bubble);
int apecss_bubble_initialize(struct APECSS_Bubble *Bubble);
int apecss_bubble_solver_initialize(struct APECSS_Bubble *Bubble);
int apecss_bubble_solver_finalize(struct APECSS_Bubble *Bubble);
int apecss_bubble_solver_run(APECSS_FLOAT tend, struct APECSS_Bubble *Bubble);
int apecss_bubble_freestruct(struct APECSS_Bubble *Bubble);
int apecss_bubble_solver_progress_initialnone();
int apecss_bubble_solver_progress_initialscreen();
int apecss_bubble_solver_progress_updatenone(int *prog, APECSS_FLOAT elapsedtime, APECSS_FLOAT totaltime);
int apecss_bubble_solver_progress_updatescreen(int *prog, APECSS_FLOAT elapsedtime, APECSS_FLOAT totaltime);
int apecss_bubble_solver_progress_finalnone();
int apecss_bubble_solver_progress_finalscreen();
APECSS_FLOAT apecss_bubble_pressure_infinity_noexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_bubble_pressure_infinity_sinexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_bubble_pressurederivative_infinity_noexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_bubble_pressurederivative_infinity_sinexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_bubble_dimensionalradius_planar(APECSS_FLOAT r);
APECSS_FLOAT apecss_bubble_dimensionalradius_cylindrical(APECSS_FLOAT r);
APECSS_FLOAT apecss_bubble_dimensionalradius_spherical(APECSS_FLOAT r);

// ---------------------
// emissions.c

int apecss_emissions_initializestruct(struct APECSS_Bubble *Bubble);
int apecss_emissions_initializenone(struct APECSS_Bubble *Bubble);
int apecss_emissions_initializelinkedlist(struct APECSS_Bubble *Bubble);
int apecss_emissions_freenone(struct APECSS_Bubble *Bubble);
int apecss_emissions_freelinkedlist(struct APECSS_Bubble *Bubble);
int apecss_emissions_updatenone(struct APECSS_Bubble *Bubble);
int apecss_emissions_updatelinkedlist(struct APECSS_Bubble *Bubble);
int apecss_emissions_addnode(struct APECSS_Bubble *Bubble);
int apecss_emissions_prunelist(struct APECSS_Bubble *Bubble);
int apecss_emissions_prune_no_node(struct APECSS_EmissionNode *Node);
int apecss_emissions_removenode(struct APECSS_Bubble *Bubble);
int apecss_emissions_advance_finitespeedincompressible(struct APECSS_Bubble *Bubble);
int apecss_emissions_advance_quasiacoustic(struct APECSS_Bubble *Bubble);
int apecss_emissions_advance_kirkwoodbethe_tait(struct APECSS_Bubble *Bubble);
int apecss_emissions_advance_kirkwoodbethe_general(struct APECSS_Bubble *Bubble);
int apecss_emissions_integrate_ev_tait_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
int apecss_emissions_integrate_ev_tait_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
int apecss_emissions_integrate_tiv_tait_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
int apecss_emissions_integrate_tiv_tait_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
int apecss_emissions_integrate_ev_general_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
int apecss_emissions_integrate_ev_general_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
int apecss_emissions_integrate_tiv_general_euler(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
int apecss_emissions_integrate_tiv_general_rk4(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Current, APECSS_FLOAT hinf);
APECSS_FLOAT apecss_emissions_f_zero(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Node);
APECSS_FLOAT apecss_emissions_f_finitespeedincompressible(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Node);
APECSS_FLOAT apecss_emissions_f_quasiacoustic(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Node);
APECSS_FLOAT apecss_emissions_f_kirkwoodbethe(struct APECSS_Bubble *Bubble, struct APECSS_EmissionNode *Node);

// ---------------------
// gas.c

int apecss_gas_setdefaultoptions(struct APECSS_Gas *Gas);
int apecss_gas_readoptions(struct APECSS_Gas *Gas, char *OptionsDir);
int apecss_gas_processoptions(struct APECSS_Gas *Gas);
APECSS_FLOAT apecss_gas_pressure_ig(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_gas_pressure_hc(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_gas_pressure_nasg(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_gas_pressurederivative_ig(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_gas_pressurederivative_hc(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_gas_pressurederivative_nasg(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_gas_density_constmass(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_gas_densityderivative_constmass(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);

// ---------------------
// interface.c

int apecss_interface_setdefaultoptions(struct APECSS_Interface *Interface);
int apecss_interface_readoptions(struct APECSS_Interface *Interface, char *OptionsDir);
int apecss_interface_processoptions(struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetension_clean(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetension_marmottant(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetension_gompertzmarmottant(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetensionderivative_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetensionderivative_gompertzmarmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetensionpressure_clean(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetensionpressure_marmottant(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetensionpressure_gompertzmarmottant(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_gompertzmarmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_pressure_viscous_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_pressure_viscous_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_pressurederivative_viscous_cleanexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_pressurederivative_viscous_marmottantexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_pressurederivative_viscous_cleanimpl(APECSS_FLOAT R, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_pressurederivative_viscous_marmottantimpl(APECSS_FLOAT R, struct APECSS_Interface *Interface);

// ---------------------
// interactions.c

int apecss_interactions_instantaneous(struct APECSS_Bubble *Bubbles[]);
int apecss_interactions_quasi_acoustic(struct APECSS_Bubble *Bubbles[]);
int apecss_interactions_cutoffdistance(struct APECSS_Bubble *Bubbles[]);

// ---------------------
// liquid.c

int apecss_liquid_setdefaultoptions(struct APECSS_Liquid *Liquid);
int apecss_liquid_readoptions(struct APECSS_Liquid *Liquid, char *OptionsDir);
int apecss_liquid_processoptions(struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_density_fixed(APECSS_FLOAT p, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_density_tait(APECSS_FLOAT p, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_density_nasg(APECSS_FLOAT p, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_soundspeed_fixed(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_soundspeed_tait(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_soundspeed_nasg(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_enthalpy_quasiacoustic(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_enthalpy_tait(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_enthalpy_nasg(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_pressure_bubblewall(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_bubblewall_kelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_bubblewall_oldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_expl(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_explkelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_exploldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_viscous(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_viscous_expl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_viscous_impl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_viscous_nonimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);

// ---------------------
// misc.c

int apecss_writeonscreen(char *str);
int apecss_erroronscreen(int num, char *message);
int apecss_infoscreen();
int apecss_helpscreen();
int apecss_readoneoption(FILE *OptionsFile, char *stuk);
int apecss_lineget(char *ssring, FILE *fp);
int apecss_linegetskip(char *ssring, FILE *fp);

// ---------------------
// odesolver.c

int apecss_odesolver_setdefaultoptions(struct APECSS_NumericsODE *NumericsODE);
int apecss_odesolver_readoptions(struct APECSS_NumericsODE *NumericsODE, char *OptionsDir);
APECSS_FLOAT apecss_odesolver(struct APECSS_Bubble *Bubble);
int apecss_odesolver_settimestep(struct APECSS_NumericsODE *ODEs, APECSS_FLOAT err, APECSS_FLOAT timetoend, APECSS_FLOAT *dt);
int apecss_odesolver_processoptions(struct APECSS_NumericsODE *ODEs);

// ---------------------
// rayleighplesset.c

APECSS_FLOAT apecss_rp_bubbleradius_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_rp_rayleighplessetvelocity_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_rp_rayleighplessetacousticrationvelocity_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_rp_kellermiksisvelocity_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_rp_gilmorevelocity_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

// ---------------------
// results.c

int apecss_results_initializestruct(struct APECSS_Bubble *Bubble);
int apecss_results_rayleighplesset_initializestruct(struct APECSS_Bubble *Bubble);
int apecss_results_emissions_initializestruct(struct APECSS_Bubble *Bubble);
int apecss_results_rayleighplesset_storenone(struct APECSS_Bubble *Bubble);
int apecss_results_rayleighplesset_storeall(struct APECSS_Bubble *Bubble);
int apecss_results_rayleighplesset_initialize(struct APECSS_Bubble *Bubble);
int apecss_results_rayleighplesset_free(struct APECSS_Bubble *Bubble);
int apecss_results_rayleighplesset_write(struct APECSS_Bubble *Bubble, int write);
int apecss_results_emissionstime_writenone(struct APECSS_Bubble *Bubble);
int apecss_results_emissionstime_writeall(struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_results_emissionstime_checknone(APECSS_FLOAT tend, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_results_emissionstime_checktime(APECSS_FLOAT tend, struct APECSS_Bubble *Bubble);
int apecss_results_emissionsspace_storenone(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsspace_storeall(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsspace_write(struct APECSS_Bubble *Bubble, int write);
int apecss_results_emissionsspace_free(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodespecific_storenone(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodespecific_storeall(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodespecific_allocnone(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodespecific_allocall(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodespecific_write(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodespecific_free(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodeminmax_identifynone(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodeminmax_identifyall(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodeminmax_allocnone(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodeminmax_allocall(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodeminmax_storenone(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodeminmax_storeall(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnodeminmax_write(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnode_initializenode(struct APECSS_ResultsEmissionsNode *Node);
int apecss_results_emissionsnode_allocnone(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnode_allocall(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnode_allocnode(struct APECSS_ResultsEmissionsNode *Node);
int apecss_results_emissionsnode_storenone(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnode_storeall(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
int apecss_results_emissionsnode_writenode(struct APECSS_ResultsEmissionsNode *Node, char *path, int digits);
int apecss_results_emissionsnode_freenode(struct APECSS_ResultsEmissionsNode *Node);

// ---------------------
// viscoelasticity.c

APECSS_FLOAT apecss_viscoelastic_zenertaurr_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_viscoelastic_zenervarsigma_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_viscoelastic_oldroydb1_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_viscoelastic_oldroydb2_ode(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

#endif  // APECSS_H_
