%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% SU2 configuration file                                                       %
% Case description: Optimization of convergent-divergent nozzle                %        
% Author: Abhishek Karn                                                        %
% Institution: IOE,Pulchowk Campus                                             %
% Date: 2021-07-23                                                             %
% File Version 7.1.1-linux64                                                   %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%
%
% Physical governing equations (EULER, NAVIER_STOKES,
%                               FEM_EULER, FEM_NAVIER_STOKES, FEM_RANS, FEM_LES,
%                               WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY,
%                               POISSON_EQUATION)
SOLVER= EULER
%
% Specify turbulence model (NONE, SA, SA_NEG, SST, SA_E, SA_COMP, SA_E_COMP)
KIND_TURB_MODEL= SST
%
% Mathematical problem (DIRECT, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT)
MATH_PROBLEM= CONTINUOUS_ADJOINT
%
% Restart solution (NO, YES)
RESTART_SOL= YES
%
% System of measurements (SI, US)
% International system of units (SI): ( meters, kilograms, Kelvins,
%                                       Newtons = kg m/s^2, Pascals = N/m^2,
%                                       Density = kg/m^3, Speed = m/s,
%                                       Equiv. Area = m^2 )
% United States customary units (US): ( inches, slug, Rankines, lbf = slug ft/s^2,
%                                       psf = lbf/ft^2, Density = slug/ft^3,
%                                       Speed = ft/s, Equiv. Area = ft^2 )
SYSTEM_MEASUREMENTS= SI
%
% -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------%
%
% Mach number (non-dimensional, based on the free-stream values)
MACH_NUMBER= 0.086
%
% Angle of attack (degrees, only for compressible flows)
AOA= 0.0
%
% Side-slip angle (degrees, only for compressible flows)
SIDESLIP_ANGLE= 0.0
%
% Init option to choose between Reynolds (default) or thermodynamics quantities
% for initializing the solution (REYNOLDS, TD_CONDITIONS)
INIT_OPTION= TD_CONDITIONS
%
% Free-stream option to choose between density and temperature (default) for
% initializing the solution (TEMPERATURE_FS, DENSITY_FS)
FREESTREAM_OPTION= TEMPERATURE_FS
%
% Free-stream pressure (101325.0 N/m^2, 2116.216 psf by default)
FREESTREAM_PRESSURE= 101325.0
%
% Free-stream temperature (288.15 K, 518.67 R by default)
FREESTREAM_TEMPERATURE= 288.15
%
% Compressible flow non-dimensionalization (DIMENSIONAL, FREESTREAM_PRESS_EQ_ONE,
%                              FREESTREAM_VEL_EQ_MACH, FREESTREAM_VEL_EQ_ONE)
REF_DIMENSIONALIZATION= DIMENSIONAL

% ---- IDEAL GAS, POLYTROPIC, VAN DER WAALS AND PENG ROBINSON CONSTANTS -------%
%
% Fluid model (STANDARD_AIR, IDEAL_GAS, VW_GAS, PR_GAS,
%              CONSTANT_DENSITY, INC_IDEAL_GAS, INC_IDEAL_GAS_POLY)
FLUID_MODEL= STANDARD_AIR
%
% Ratio of specific heats (1.4 default and the value is hardcoded
%                          for the model STANDARD_AIR, compressible only)
GAMMA_VALUE= 1.4
%
% Specific gas constant (287.058 J/kg*K default and this value is hardcoded
%                        for the model STANDARD_AIR, compressible only)
GAS_CONSTANT= 287.058
%
% Critical Temperature (131.00 K by default)
CRITICAL_TEMPERATURE= 131.00
%
% Critical Pressure (3588550.0 N/m^2 by default)
CRITICAL_PRESSURE= 3588550
%
% Acentric factor (0.035 (air))
ACENTRIC_FACTOR= 0.035

% --------------------------- VISCOSITY MODEL ---------------------------------%
%
% Viscosity model (SUTHERLAND, CONSTANT_VISCOSITY, POLYNOMIAL_VISCOSITY).
VISCOSITY_MODEL= CONSTANT_VISCOSITY
%
% Molecular Viscosity that would be constant (1.716E-5 by default)
MU_CONSTANT= 1.716E-05

% --------------------------- THERMAL CONDUCTIVITY MODEL ----------------------%
%
% Laminar Conductivity model (CONSTANT_CONDUCTIVITY, CONSTANT_PRANDTL,
% POLYNOMIAL_CONDUCTIVITY).
CONDUCTIVITY_MODEL= CONSTANT_PRANDTL
%
% Molecular Thermal Conductivity that would be constant (0.0257 by default)
KT_CONSTANT= 0.0257

% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%
%
% Navier-Stokes (no-slip), constant heat flux wall  marker(s) (NONE = no marker)
% Format: ( marker name, constant heat flux (J/m^2), ... )
MARKER_EULER= ( WALL)
%
% Symmetry boundary marker(s) (NONE = no marker)
MARKER_SYM= ( SYMMETRY)
%
% Riemann boundary marker(s) (NONE = no marker)
% Format: (marker, data kind flag, list of data)
MARKER_RIEMANN= ( INFLOW, TOTAL_CONDITIONS_PT, 101831.3, 286.75, 1.0, 0.0, 0.0, OUTFLOW, STATIC_PRESSURE, 500, 0.0, 0.0, 0.0, 0.0 )

% -------------------- TURBULENT NUMERICAL METHOD DEFINITION ------------------%
%
% Convective numerical method (SCALAR_UPWIND)
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
%
% Time discretization (EULER_IMPLICIT)
TIME_DISCRE_TURB= EULER_IMPLICIT
%
% Reduction factor of the CFL coefficient in the turbulence problem
CFL_REDUCTION_TURB= 0.5

% ------------------------ SURFACES IDENTIFICATION ----------------------------%
%
% Marker of the surface which is going to be plotted or designed
MARKER_PLOTTING= ( WALL)
%
% Marker of the surface where the functional (Cd, Cl, etc.) will be evaluated
MARKER_MONITORING= (OUTFLOW)
%



% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%
%
% Numerical method for spatial gradients (GREEN_GAUSS, LEAST_SQUARES, 
%                                         WEIGHTED_LEAST_SQUARES)
NUM_METHOD_GRAD= WEIGHTED_LEAST_SQUARES
%
% Courant-Friedrichs-Lewy condition of the finest grid
CFL_NUMBER= 0.5
%
% Adaptive CFL number (NO, YES)
CFL_ADAPT= NO
%
% Parameters of the adaptive CFL number (factor down, factor up, CFL min value,
%                                        CFL max value )
CFL_ADAPT_PARAM= ( 0.5, 1, 1.0, 100.0 )
%
% Runge-Kutta alpha coefficients
RK_ALPHA_COEFF= ( 0.66667, 0.66667, 1.000000 )
%
% Linear solver for the implicit formulation (BCGSTAB, FGMRES)
LINEAR_SOLVER= BCGSTAB
%
% Min error of the linear solver for the implicit formulation
LINEAR_SOLVER_ERROR= 1E-6
%
% Max number of iterations of the linear solver for the implicit formulation
LINEAR_SOLVER_ITER= 5
%
% Preconditioner of the Krylov linear solver (ILU, LU_SGS, LINELET, JACOBI)
LINEAR_SOLVER_PREC= ILU
%
% Linael solver ILU preconditioner fill-in level (1 by default)
LINEAR_SOLVER_ILU_FILL_IN= 0
% -------------------- FLOW NUMERICAL METHOD DEFINITION -----------------------%
%
% Convective numerical method (JST, LAX-FRIEDRICH, CUSP, ROE, AUSM, HLLC,
%                              TURKEL_PREC, MSW)
CONV_NUM_METHOD_FLOW= ROE
%
% Coefficient for the limiter (smooth regions)
VENKAT_LIMITER_COEFF= 0.006
%
% 2nd and 4th order artificial dissipation coefficients
JST_SENSOR_COEFF= ( 0.5, 0.02 )
%
% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT, EULER_EXPLICIT)
TIME_DISCRE_FLOW= EULER_IMPLICIT
% ----------- SLOPE LIMITER AND DISSIPATION SENSOR DEFINITION -----------------%
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the flow equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_FLOW= YES
%
% Slope limiter (NONE, VENKATAKRISHNAN, VENKATAKRISHNAN_WANG,
%                BARTH_JESPERSEN, VAN_ALBADA_EDGE)
SLOPE_LIMITER_FLOW= NONE
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the turbulence equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_TURB= NO


% -------------------------- MULTIGRID PARAMETERS -----------------------------%
%
% Multi-Grid Levels (0 = no multi-grid)
MGLEVEL= 2
%
% Multi-grid cycle (V_CYCLE, W_CYCLE, FULLMG_CYCLE)
MGCYCLE= V_CYCLE
%
% Multi-grid pre-smoothing level
MG_PRE_SMOOTH= ( 1, 2, 3, 3 )
%
% Multi-grid post-smoothing level
MG_POST_SMOOTH= ( 0, 0, 0, 0 )
%
% Jacobi implicit smoothing of the correction
MG_CORRECTION_SMOOTH= ( 0, 0, 0, 0 )
%
% Damping factor for the residual restriction
MG_DAMP_RESTRICTION= 0.8
%
% Damping factor for the correction prolongation
MG_DAMP_PROLONGATION= 0.8

% ---------------- ADJOINT-FLOW NUMERICAL METHOD DEFINITION -------------------%
%
% Convective numerical method (JST, LAX-FRIEDRICH, ROE)
CONV_NUM_METHOD_ADJFLOW= ROE
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the adjoint flow equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_ADJFLOW= YES
%
% Slope limiter (NONE, VENKATAKRISHNAN, BARTH_JESPERSEN, VAN_ALBADA_EDGE,
%                SHARP_EDGES, WALL_DISTANCE)
SLOPE_LIMITER_ADJFLOW= NONE
%
% 2nd, and 4th order artificial dissipation coefficients
ADJ_JST_SENSOR_COEFF= ( 0.5, 0.02 )
%
% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT)
TIME_DISCRE_ADJFLOW= EULER_IMPLICIT
%
% Relaxation coefficient
RELAXATION_FACTOR_ADJOINT= 1.0
%
% Reduction factor of the CFL coefficient in the adjoint problem
CFL_REDUCTION_ADJFLOW= 0.8
%
% Limit value for the adjoint variable
LIMIT_ADJFLOW= 1E15
%
% Multigrid adjoint problem (NO, YES)
MG_ADJFLOW= NO
%
% Objective function in gradient evaluation  (DRAG, LIFT, SIDEFORCE, MOMENT_X,
%                                             MOMENT_Y, MOMENT_Z, EFFICIENCY,
%                                             EQUIVALENT_AREA, NEARFIELD_PRESSURE,
%                                             FORCE_X, FORCE_Y, FORCE_Z, THRUST,
%                                             TORQUE, FREE_SURFACE, TOTAL_HEATFLUX,
%                                             MAXIMUM_HEATFLUX, INVERSE_DESIGN_PRESSURE,
%                                             INVERSE_DESIGN_HEATFLUX, SURFACE_TOTAL_PRESSURE, 
%                                             SURFACE_MASSFLOW)
% For a weighted sum of objectives: separate by commas, add OBJECTIVE_WEIGHT and MARKER_MONITORING in matching order.
OBJECTIVE_FUNCTION = SURFACE_TOTAL_PRESSURE
%
% List of weighting values when using more than one OBJECTIVE_FUNCTION. Separate by commas and match with MARKER_MONITORING.
OBJECTIVE_WEIGHT= 1
% Marker on which to track one-dimensionalized quantities
MARKER_ANALYZE  = (OUTFLOW)
%
% Method to compute the average value in MARKER_ANALYZE (AREA, MASSFLUX).
MARKER_ANALYZE_AVERAGE = AREA




% --------------------------- CONVERGENCE PARAMETERS --------------------------%
%
% Number of total iterations
ITER= 1000
%
% Convergence criteria (CAUCHY, RESIDUAL)
%
CONV_CRITERIA= RESIDUAL
%
%
% Min value of the residual (log10 of the residual)
CONV_RESIDUAL_MINVAL= -6
%
% Start convergence criteria at iteration number
CONV_STARTITER= 10

% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%
%
% Mesh input file
MESH_FILENAME= half_nozzle6.su2
%
% Mesh input file format (SU2, CGNS)
MESH_FORMAT= SU2
%
% Mesh output file
MESH_OUT_FILENAME= mesh_out.su2
%
% Restart flow input file
SOLUTION_FILENAME= solution_flow.dat
%
% Restart adjoint input file
SOLUTION_ADJ_FILENAME= solution_adj.dat
%
% Output file format (TECPLOT, TECPLOT_BINARY, PARAVIEW, PARAVIEW_BINARY,
%                     FIELDVIEW, FIELDVIEW_BINARY)
TABULAR_FORMAT= CSV
%
% Output file convergence history (w/o extension)
CONV_FILENAME= history
%
% Output file restart flow
RESTART_FILENAME= restart_flow.dat
%
% Output file flow (w/o extension) variables
VOLUME_FILENAME= flow
%
% Output file restart adjoint
RESTART_ADJ_FILENAME= restart_adj.dat
%
% Output file adjoint (w/o extension) variables
VOLUME_ADJ_FILENAME= adjoint
%
% Output Objective function gradient (using continuous adjoint)
GRAD_OBJFUNC_FILENAME= of_grad.dat
%
% Output file surface flow coefficient (w/o extension)
SURFACE_FILENAME= surface_flow
%
% Output file surface adjoint coefficient (w/o extension)
SURFACE_ADJ_FILENAME= surface_adjoint
%
% Writing solution file frequency
OUTPUT_WRT_FREQ= 1000
%
% Screen output
SCREEN_OUTPUT= (INNER_ITER, RMS_DENSITY, RMS_TKE, RMS_DISSIPATION, CL,CD)
% Output files
OUTPUT_FILES = (RESTART, PARAVIEW, SURFACE_PARAVIEW, SURFACE_CSV)
% -------------------- FREE-FORM DEFORMATION PARAMETERS -----------------------%
%
% Tolerance of the Free-Form Deformation point inversion
FFD_TOLERANCE= 1E-10
%
% Maximum number of iterations in the Free-Form Deformation point inversion
FFD_ITERATIONS= 500
%
% FFD box definition: 3D case (FFD_BoxTag, X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4,
%                              X5, Y5, Z5, X6, Y6, Z6, X7, Y7, Z7, X8, Y8, Z8)
%                     2D case (FFD_BoxTag, X1, Y1, 0.0, X2, Y2, 0.0, X3, Y3, 0.0, X4, Y4, 0.0,
%                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
FFD_DEFINITION= (MAIN_BOX, -0.04, 0.01, 0, 0.04, 0.01, 0, 0.04, 0.05,  0, -0.04, 0.05,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
%
% FFD box degree: 3D case (x_degree, y_degree, z_degree)
%                 2D case (x_degree, y_degree, 0)
FFD_DEGREE= (12, 1,0)
%
% Surface continuity at the intersection with the FFD (1ST_DERIVATIVE, 2ND_DERIVATIVE)
FFD_CONTINUITY= 2ND_DERIVATIVE
% ------------------------ GRID DEFORMATION PARAMETERS ------------------------%
%
% Kind of deformation (FFD_SETTING, HICKS_HENNE, HICKS_HENNE_NORMAL, PARABOLIC,
%                      HICKS_HENNE_SHOCK, NACA_4DIGITS, DISPLACEMENT, ROTATION, 
%                      FFD_CONTROL_POINT, FFD_DIHEDRAL_ANGLE, FFD_TWIST_ANGLE, 
%                      FFD_ROTATION)
% Marker of the surface in which we are going apply the shape deformation
DV_MARKER= ( WALL)
%
DV_KIND= FFD_CONTROL_POINT_2D
%
%
% Parameters of the shape deformation 
% 	- HICKS_HENNE_FAMILY ( Lower(0)/Upper(1) side, x_Loc )
% 	- NACA_4DIGITS ( 1st digit, 2nd digit, 3rd and 4th digit )
% 	- PARABOLIC ( 1st digit, 2nd and 3rd digit )
% 	- DISPLACEMENT ( x_Disp, y_Disp, z_Disp )
%       - FFD_CONTROL_POINT ( FFD_BoxTag, i_Ind, j_Ind, k_Ind, x_Disp, y_Disp, z_Disp )
% 	- ROTATION ( x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
DV_PARAM= (MAIN_BOX, 6,0, 0,0,1.0,0 )
%
% Value of the shape deformation deformation
DV_VALUE= 0.0
% Number of smoothing iterations for FEA mesh deformation
DEFORM_LINEAR_SOLVER_ITER= 500
%
% Number of nonlinear deformation iterations (surface deformation increments)
DEFORM_NONLINEAR_ITER= 1
%
% Print the residuals during mesh deformation to the console (YES, NO)
DEFORM_CONSOLE_OUTPUT= YES
%
% Minimum residual criteria for the linear solver convergence of grid deformation
DEFORM_LINEAR_SOLVER_ERROR= 1E-14
%
% Type of element stiffness imposed for FEA mesh deformation (INVERSE_VOLUME,
%                                          WALL_DISTANCE, CONSTANT_STIFFNESS)
DEFORM_STIFFNESS_TYPE= INVERSE_VOLUME


% --------------------- OPTIMAL SHAPE DESIGN DEFINITION -----------------------%
%
% Available flow based objective functions or constraint functions
%    DRAG, LIFT, SIDEFORCE, EFFICIENCY, BUFFET, 
%    FORCE_X, FORCE_Y, FORCE_Z,
%    MOMENT_X, MOMENT_Y, MOMENT_Z,
%    THRUST, TORQUE, FIGURE_OF_MERIT,
%    EQUIVALENT_AREA, NEARFIELD_PRESSURE,
%    TOTAL_HEATFLUX, MAXIMUM_HEATFLUX,
%    INVERSE_DESIGN_PRESSURE, INVERSE_DESIGN_HEATFLUX,
%    SURFACE_TOTAL_PRESSURE, SURFACE_MASSFLOW
%    SURFACE_STATIC_PRESSURE, SURFACE_MACH
%
% Available geometrical based objective functions or constraint functions
%    AIRFOIL_AREA, AIRFOIL_THICKNESS, AIRFOIL_CHORD, AIRFOIL_TOC, AIRFOIL_AOA,
%    WING_VOLUME, WING_MIN_THICKNESS, WING_MAX_THICKNESS, WING_MAX_CHORD, WING_MIN_TOC, WING_MAX_TWIST, WING_MAX_CURVATURE, WING_MAX_DIHEDRAL
%    STATION#_WIDTH, STATION#_AREA, STATION#_THICKNESS, STATION#_CHORD, STATION#_TOC,
%    STATION#_TWIST (where # is the index of the station defined in GEO_LOCATION_STATIONS)
%
% Available design variables
% 2D Design variables
%    FFD_CONTROL_POINT_2D   (  19, Scale | Mark. List | FFD_BoxTag, i_Ind, j_Ind, x_Mov, y_Mov )
%    FFD_CAMBER_2D          (  20, Scale | Mark. List | FFD_BoxTag, i_Ind )
%    FFD_THICKNESS_2D       (  21, Scale | Mark. List | FFD_BoxTag, i_Ind )
%    FFD_TWIST_2D           (  22, Scale | Mark. List | FFD_BoxTag, x_Orig, y_Orig )
%    HICKS_HENNE            (  30, Scale | Mark. List | Lower(0)/Upper(1) side, x_Loc )
%    ANGLE_OF_ATTACK        ( 101, Scale | Mark. List | 1.0 )
%
% 3D Design variables
%    FFD_CONTROL_POINT      (  11, Scale | Mark. List | FFD_BoxTag, i_Ind, j_Ind, k_Ind, x_Mov, y_Mov, z_Mov )
%    FFD_NACELLE            (  12, Scale | Mark. List | FFD_BoxTag, rho_Ind, theta_Ind, phi_Ind, rho_Mov, phi_Mov )
%    FFD_GULL               (  13, Scale | Mark. List | FFD_BoxTag, j_Ind )
%    FFD_CAMBER             (  14, Scale | Mark. List | FFD_BoxTag, i_Ind, j_Ind )
%    FFD_TWIST              (  15, Scale | Mark. List | FFD_BoxTag, j_Ind, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
%    FFD_THICKNESS          (  16, Scale | Mark. List | FFD_BoxTag, i_Ind, j_Ind )
%    FFD_ROTATION           (  18, Scale | Mark. List | FFD_BoxTag, x_Axis, y_Axis, z_Axis, x_Turn, y_Turn, z_Turn )
%    FFD_ANGLE_OF_ATTACK    (  24, Scale | Mark. List | FFD_BoxTag, 1.0 )
%
% Global design variables
%    TRANSLATION            (   1, Scale | Mark. List | x_Disp, y_Disp, z_Disp )
%    ROTATION               (   2, Scale | Mark. List | x_Axis, y_Axis, z_Axis, x_Turn, y_Turn, z_Turn )
%
% Optimization objective function with scaling factor, separated by semicolons.
% To include quadratic penalty function: use OPT_CONSTRAINT option syntax within the OPT_OBJECTIVE list.
% ex= Objective * Scale
OPT_OBJECTIVE= SURFACE_TOTAL_PRESSURE
%
% Optimization constraint functions with pushing factors (affects its value, not the gradient  
% in the python scripts), separated by semicolons
% ex= (Objective = Value ) * Scale, use '>','<','='
OPT_CONSTRAINT= NONE
%
% Factor to reduce the norm of the gradient (affects the objective function and gradient in the python scripts)
% In general, a norm of the gradient ~1E-6 is desired.
OPT_GRADIENT_FACTOR= 1E-6
%
% Factor to relax or accelerate the optimizer convergence (affects the line search in SU2_DEF)
% In general, surface deformations of 0.01'' or 0.0001m are desirable
OPT_RELAX_FACTOR= 1E3
%
% Maximum number of optimizer iterations
OPT_ITERATIONS= 100
%
% Requested accuracy
OPT_ACCURACY= 1E-10
%
% Upper bound for each design variable
OPT_BOUND_UPPER= 0.020
%
% Lower bound for each design variable
OPT_BOUND_LOWER= -0.020
%
% Optimization design variables, separated by semicolons
DEFINITION_DV= ( 19, 1.0| WALL| MAIN_BOX, 1,1,0,1.0);( 19, 1.0| WALL | MAIN_BOX, 2,1,0,1.0);( 19, 1.0| WALL| MAIN_BOX, 3,1,0,1.0);( 19, 1.0| WALL | MAIN_BOX, 4,1,0,1.0);( 19, 1.0| WALL | MAIN_BOX, 5,1,0,1.0);( 19, 1.0| WALL | MAIN_BOX, 6,1,0,1.0);( 19, 1.0| WALL| MAIN_BOX, 7,1,0,1.0);( 19, 1.0| WALL | MAIN_BOX, 8,1,0,1.0);( 19, 1.0| WALL | MAIN_BOX, 9,1,0,1.0);( 19, 1.0| WALL | MAIN_BOX, 10,1,0,1.0)
