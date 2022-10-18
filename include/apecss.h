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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef APECSS_H_
#define APECSS_H_

#define APECSS_VERSION_NUM (0.6)
static const char APECSS_RELEASE_DATE[] = "16-Oct-2022";

// -------------------------------------------------------------------
// CONSTANTS & MACROS
// -------------------------------------------------------------------

// Floating point precision
// #define APECSS_PRECISION_LONGDOUBLE

// Data types
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

// Runge-Kutta scheme
#define APECSS_RK54_7M (0)  // RK5(4)7M (minimum truncation) coefficients of Dormand & Prince (1980)
#define APECSS_RK54_7S (1)  // RK5(4)7S (stability optimized) coefficients of Dormand & Prince (1980)

// Bubble model (bit-wise)
#define APECSS_BUBBLEMODEL_RP (1)  // Standard Rayleigh-Plesset model
#define APECSS_BUBBLEMODEL_RP_ACOUSTICRADIATION (2)  // Rayleigh-Plesset model incl. acoustic radiation term
#define APECSS_BUBBLEMODEL_KELLERMIKSIS (4)  // Keller-Miksis model
#define APECSS_BUBBLEMODEL_GILMORE (8)  // Gilmore model

// Gas equations of state (bit-wise)
#define APECSS_GAS_IG (1)  // ideal gas
#define APECSS_GAS_HC (2)  // hard-core gas
#define APECSS_GAS_NASG (4)  // Noble-Abel stiffend gas

// Liquid models (bit-wise)
#define APECSS_LIQUID_NEWTONIAN (1)  // Newtonian fluid
#define APECSS_LIQUID_KELVINVOIGT (2)  // Kelvin-Voigt solid
#define APECSS_LIQUID_ZENER (4)  // Zener solid, requires two additional ODEs
#define APECSS_LIQUID_OLDROYDB (8)  // Oldroyd-B liquid, requires two additional ODEs

// Lipid coating model (bit-wise)
#define APECSS_LIPIDCOATING_NONE (1)
#define APECSS_LIPIDCOATING_MARMOTTANT (2)
#define APECSS_LIPIDCOATING_GOMPERTZFUNCTION (4)

// Excitation types (bit-wise)
#define APECSS_EXCITATION_NONE (0)
#define APECSS_EXCITATION_SIN (1)

// Emission type (bit-wise)
#define APECSS_EMISSION_NONE (0)
#define APECSS_EMISSION_INCOMPRESSIBLE (1)
#define APECSS_EMISSION_FINITE_TIME_INCOMPRESSIBLE (2)
#define APECSS_EMISSION_QUASIACOUSTIC (4)
#define APECSS_EMISSION_KIRKWOODBETHE (8)

// Misc
#define APECSS_DATA_ALLOC_INCREMENT (10000)
#define APECSS_STRINGLENGTH (512)
#define APECSS_STRINGLENGTH_SPRINTF (1024)
#define APECSS_STRINGLENGTH_SPRINTF_LONG (2048)

// Debugging aids
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
  int EOS;

  APECSS_FLOAT Gamma;  // Polytropic exponent
  APECSS_FLOAT B;  // Attractive pressure constant [Pa]
  APECSS_FLOAT b;  // Co-volume [m^3/kg]
  APECSS_FLOAT h;  // van der Waals hardcore radius [m]
  APECSS_FLOAT dmol;  // Molecular diameter [m]
  APECSS_FLOAT mmol;  // Molecular weight [kg/mol]
  APECSS_FLOAT pref;  // Reference pressure [Pa]
  APECSS_FLOAT rhoref;  // Reference density [kg/m^3]
  APECSS_FLOAT Kref;  // Constant reference coefficient of a polytropic EOS

  // Pointers to the functions describing the gas pressure and its derivative
  APECSS_FLOAT (*get_pressure)(APECSS_FLOAT *Sol, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative)(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
};

struct APECSS_Liquid
{
  int Type;

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
  APECSS_FLOAT Kref;  // Constant reference coefficient of a polytropic EOS
  APECSS_FLOAT G;  // Shear modulus
  APECSS_FLOAT eta;  // Polymer viscosity
  APECSS_FLOAT lambda;  // Relaxation time

  APECSS_FLOAT (*get_density)(APECSS_FLOAT p, struct APECSS_Liquid *Liquid);
  APECSS_FLOAT (*get_soundspeed)(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
  APECSS_FLOAT (*get_enthalpy)(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);

  // Pointers to the functions describing the liquid pressure and its derivative
  APECSS_FLOAT (*get_pressure_bubblewall)(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_bubblewall_expl)(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);

  // Pointers to the functions describing the pressure and its derivative due to viscous contributions
  APECSS_FLOAT (*get_pressure_viscous)(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_viscous_expl)(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_viscous_impl)(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);

  // Pointers to the functions describing the state of the liquid at infinity
  APECSS_FLOAT (*get_pressure_infinity)(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_pressurederivative_infinity)(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
};

struct APECSS_Interface
{
  APECSS_FLOAT sigma;  // Surface tension coefficient [N/m]

  // Lipid coating properties of the Marmottant and Marmottant-Gompertz models
  int LipidCoatingModel;  // Coating model for the lipid coating of the interface
  APECSS_FLOAT Elasticity;  // Surface elasticity
  APECSS_FLOAT Viscosity;  // Dilatational viscosity
  APECSS_FLOAT sigma0;  // Initial surface tension coefficient [N/m]
  APECSS_FLOAT Rbuck;  // Buckling radius [m]
  APECSS_FLOAT Rrupt;  // Rupture radius [m]
  APECSS_FLOAT GompertzB;  // B coefficient of the Marmottant-Gompertz model
  APECSS_FLOAT GompertzC;  // C coefficient of the Marmottant-Gompertz model

  APECSS_FLOAT (*get_surfacetension)(APECSS_FLOAT R, struct APECSS_Interface *Interface);

  // Pointers to the functions describing the pressure and its derivative due to surface tension
  APECSS_FLOAT (*get_pressure_surfacetension)(APECSS_FLOAT R, struct APECSS_Interface *Interface);
  APECSS_FLOAT (*get_pressurederivative_surfacetension)(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
};

struct APECSS_EmissionNode
{
  int id;
  APECSS_FLOAT r;  // Radial coordinate
  APECSS_FLOAT p;  // Pressure
  APECSS_FLOAT u;  // Velocity
  APECSS_FLOAT f;  // = const.
  APECSS_FLOAT g;  // = const.
  struct APECSS_EmissionNode *forward;
  struct APECSS_EmissionNode *backward;
};

struct APECSS_Emissions
{
  int Type;
  APECSS_FLOAT CutOffDistance;
  APECSS_FLOAT KB_IterTolerance;

  int nNodes;
  struct APECSS_EmissionNode *FirstNode;  // First node of the linked list
  struct APECSS_EmissionNode *LastNode;  // Last node of the linked list
  int (*advance)(struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*get_advectingvelocity)(APECSS_FLOAT u);
};

struct APECSS_NumericsODE
{
  int RKtype;  // Runge-Kutta scheme

  APECSS_FLOAT dtMin;  // Minimum allowable time-step
  APECSS_FLOAT dtMax;  // Maximum allowable time-step
  APECSS_FLOAT minScale;  // Minimum allowable value by which the dt-controller multipliesdt
  APECSS_FLOAT maxScale;  // Maximum allowable value by which the dt-controller multipliesdt
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
  APECSS_FLOAT *cL;  // Speed of sound in the liquid at the bubble wall [m/s]

  int nUserODEs;  // The solutions of the first nUserODEs additional ODEs defined by the user are written to file
  APECSS_FLOAT **UserODEsSol;
  char **UserODEsName;
};

struct APECSS_ResultsEmissionsSpace
{
  int n;  // Number of stored results
  int nAllocated;  // Length of the arrays
  APECSS_FLOAT RadialLocation;
  APECSS_FLOAT *t;
  APECSS_FLOAT *p;
  APECSS_FLOAT *u;
  APECSS_FLOAT *c;
  APECSS_FLOAT *pInf;
};

struct APECSS_ResultsEmissionsNode
{
  int id;
  int n;  // Number of stored results
  int nAllocated;  // Length of the arrays
  int *real_id;
  APECSS_FLOAT *r;
  APECSS_FLOAT *t;
  APECSS_FLOAT *p;
  APECSS_FLOAT *u;
  APECSS_FLOAT *c;
  APECSS_FLOAT *pInf;
};

struct APECSS_ResultsEmissions
{
  int nTimeInstances;
  int nextTimeInstance;
  APECSS_FLOAT *TimeInstances;

  int freqSpaceLocations;  // Frequency with which the results are stored (with respect to the time-step number)
  int nSpaceLocations;
  struct APECSS_ResultsEmissionsSpace *SpaceLocation;

  int nNodes;
  struct APECSS_ResultsEmissionsNode *Node;

  int MinMaxPeriod;
  struct APECSS_ResultsEmissionsNode *Node_Rmin, *Node_Umin, *Node_pLmax;
  APECSS_FLOAT Rmin, Umin, pLmax;

  int (*allocation_nodespecific)(struct APECSS_Bubble *Bubble);
  int (*allocation_nodeminmax)(struct APECSS_Bubble *Bubble);

  int (*store_nodespecific)(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
  int (*store_nodeminmax)(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);
};

struct APECSS_Results
{
  int digits;  // Number of digits for the output
  char dir[APECSS_STRINGLENGTH_SPRINTF];  // Directory path where the results are stored

  struct APECSS_ResultsBubble *RayleighPlesset;
  struct APECSS_ResultsEmissions *Emissions;
};

struct APECSS_Bubble
{
  int id;  // Bubble ID
  APECSS_FLOAT tStart;  // Start time [s]
  APECSS_FLOAT tEnd;  // End time [s]

  APECSS_FLOAT dt;  // Time-step
  int dtNumber;  // Time-step number
  int nSubIter;  // Total number of sub-iterations to control the error
  APECSS_FLOAT *k2, *k3, *k4, *k5, *k6, *k7, *kLast;

  int RPModel;  // Model governing the bubble dynamics

  // Primary variables
  APECSS_FLOAT t;  // Physical time [s]
  APECSS_FLOAT R;  // Radius of the bubble [m]
  APECSS_FLOAT U;  // Velocity of the bubble wall [m/s]

  // Ambient and initial conditions
  APECSS_FLOAT p0;  // Ambient pressure [Pa]
  APECSS_FLOAT T0;  // Ambient temperature [K]
  APECSS_FLOAT R0;  // Initial radius [m]
  APECSS_FLOAT U0;  // Initial velocity [m/s]
  APECSS_FLOAT pG0;  // Initial gas pressure [Pa]
  APECSS_FLOAT rhoG0;  // Initial gas density [kg/m^3]

  // ODEs
  int nODEs;  // Total number of ODEs
  int nUserODEs;  // Number of additional user-defined ODEs
  APECSS_FLOAT *ODEsSol;  // Solution of each ODE
  APECSS_FLOAT (**ode)(APECSS_FLOAT *, APECSS_FLOAT, struct APECSS_Bubble *);  // Array of pointers to the functions containing the ODEs
  struct APECSS_NumericsODE *NumericsODE;

  // Properties of the fluids and interface
  struct APECSS_Gas *Gas;
  struct APECSS_Liquid *Liquid;
  struct APECSS_Interface *Interface;

  // Excitation of the bubble
  struct APECSS_Excitation *Excitation;

  // Lagrangian emissions
  struct APECSS_Emissions *Emissions;
  int (*emissions_initialize)(struct APECSS_Bubble *Bubble);
  int (*emissions_update)(struct APECSS_Bubble *Bubble);
  int (*emissions_free)(struct APECSS_Bubble *Bubble);

  // Results
  struct APECSS_Results *Results;
  int (*results_rayleighplesset_store)(struct APECSS_Bubble *Bubble);
  int (*results_emissionstime_write)(struct APECSS_Bubble *Bubble);
  APECSS_FLOAT (*results_emissionstime_check)(struct APECSS_Bubble *Bubble);
  int (*results_emissionsspace_store)(struct APECSS_Bubble *Bubble);
  int (*results_emissionsnodeminmax_identify)(struct APECSS_Bubble *Bubble);
  int (*results_emissionsnode_alloc)(struct APECSS_Bubble *Bubble);
  int (*results_emissionsnode_store)(struct APECSS_EmissionNode *Node, APECSS_FLOAT c, APECSS_FLOAT pinf, struct APECSS_Bubble *Bubble);

  // Progress screen
  int (*progress_initial)();
  int (*progress_update)(int *prog, APECSS_FLOAT t, APECSS_FLOAT totaltime);
  int (*progress_final)();
};

// -------------------------------------------------------------------
// FUNCTION DECLARATIONS
// -------------------------------------------------------------------

// ---------------------
// bubble.c

int apecss_bubble_setdefaultoptions(struct APECSS_Bubble *Bubble);
int apecss_bubble_processoptions(struct APECSS_Bubble *Bubble);
int apecss_bubble_initialize(struct APECSS_Bubble *Bubble);
int apecss_bubble_solve(struct APECSS_Bubble *Bubble);
int apecss_bubble_freearrays(struct APECSS_Bubble *Bubble);
int apecss_bubble_solverprogress_initialnone();
int apecss_bubble_solverprogress_initialscreen();
int apecss_bubble_solverprogress_updatenone(int *prog, APECSS_FLOAT t, APECSS_FLOAT totaltime);
int apecss_bubble_solverprogress_updatescreen(int *prog, APECSS_FLOAT t, APECSS_FLOAT totaltime);
int apecss_bubble_solverprogress_finalnone();
int apecss_bubble_solverprogress_finalscreen();

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
int apecss_emissions_removenode(struct APECSS_Bubble *Bubble);
int apecss_emissions_advance_finitetimeincompressible(struct APECSS_Bubble *Bubble);
int apecss_emissions_advance_quasiacoustic(struct APECSS_Bubble *Bubble);
int apecss_emissions_advance_kirkwoodbethe(struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_emissions_getadvectingvelocity_returnzero(APECSS_FLOAT u);
APECSS_FLOAT apecss_emissions_getadvectingvelocity_returnvelocity(APECSS_FLOAT u);

// ---------------------
// gas.c

int apecss_gas_setdefaultoptions(struct APECSS_Bubble *Bubble);
int apecss_gas_processoptions(struct APECSS_Bubble *Bubble);
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

int apecss_interface_setdefaultoptions(struct APECSS_Bubble *Bubble);
int apecss_interface_processoptions(struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_interface_surfacetension_clean(APECSS_FLOAT R, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetension_marmottant(APECSS_FLOAT R, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetension_gompertzmarmottant(APECSS_FLOAT R, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetensionderivative_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetensionderivative_gompertzmarmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetensionpressure_clean(APECSS_FLOAT R, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetensionpressure_marmottant(APECSS_FLOAT R, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetensionpressure_gompertzmarmottant(APECSS_FLOAT R, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);
APECSS_FLOAT apecss_interface_surfacetensionpressurederivative_gompertzmarmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Interface *Interface);

// ---------------------
// liquid.c

int apecss_liquid_setdefaultoptions(struct APECSS_Bubble *Bubble);
int apecss_liquid_processoptions(struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_density_fixed(APECSS_FLOAT p, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_density_nasg(APECSS_FLOAT p, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_soundspeed_fixed(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_soundspeed_nasg(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_enthalpy_quasiacoustic(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_enthalpy_nasg(APECSS_FLOAT p, APECSS_FLOAT rho, struct APECSS_Liquid *Liquid);
APECSS_FLOAT apecss_liquid_pressure_bubblewall(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_bubblewall_kelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_bubblewall_oldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_expl(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_explkelvinvoigt(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_zener(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_bubblewall_exploldroydb(APECSS_FLOAT *Sol, APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_infinity_noexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_infinity_sinexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_infinity_noexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_infinity_sinexcitation(APECSS_FLOAT t, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_viscous_clean(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressure_viscous_marmottant(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_viscous_cleanexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_viscous_marmottantexpl(APECSS_FLOAT R, APECSS_FLOAT U, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_viscous_nonimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_viscous_cleanimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_liquid_pressurederivative_viscous_marmottantimpl(APECSS_FLOAT R, struct APECSS_Bubble *Bubble);

// ---------------------
// odesolver.c

APECSS_FLOAT apecss_odesolver(struct APECSS_Bubble *Bubble);
int apecss_odesolver_settimestep(struct APECSS_NumericsODE *ODEs, APECSS_FLOAT err, APECSS_FLOAT timetoend, APECSS_FLOAT *dt);
int apecss_odesolver_rungekuttacoeffs(struct APECSS_NumericsODE *ODEs, int nODEs);

// ---------------------
// onscreen.c

int apecss_writeonscreen(char *str);
int apecss_erroronscreen(int num, char *message);
int apecss_infoscreen();

// ---------------------
// options.c

int apecss_options_setdefault(struct APECSS_Bubble *Bubble);
int apecss_options_process(struct APECSS_Bubble *Bubble);
int apecss_options_readfile(struct APECSS_Bubble *Bubble, char *OptionsDir);
int apecss_readoneoption(FILE *OptionsFile, char *stuk);
int apecss_lineget(char *ssring, FILE *fp);
int apecss_linegetsemi(char *ssring, FILE *fp);
int apecss_linegetskip(char *ssring, FILE *fp);

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
int apecss_results_rayleighplesset_write(struct APECSS_Bubble *Bubble);
int apecss_results_emissionstime_writenone(struct APECSS_Bubble *Bubble);
int apecss_results_emissionstime_writeall(struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_results_emissionstime_checknone(struct APECSS_Bubble *Bubble);
APECSS_FLOAT apecss_results_emissionstime_checktime(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsspace_storenone(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsspace_storeall(struct APECSS_Bubble *Bubble);
int apecss_results_emissionsspace_write(struct APECSS_Bubble *Bubble);
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