#define  PHYSICS                 RHD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                CYLINDRICAL
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           PARABOLIC
#define  TIME_STEPPING           CHARACTERISTIC_TRACING
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     4
#define  USER_DEF_CONSTANTS      1

/* -- physics dependent declarations -- */

#define  EOS                     TAUB
#define  ENTROPY_SWITCH          NO
#define  USE_FOUR_VELOCITY       NO

/* -- user-defined parameters (labels) -- */

#define  BETA                    0
#define  RHO_IN                  1
#define  RHO_OUT                 2
#define  PRESS_IN                3

/* -- user-defined symbolic constants -- */

#define  EPS_PSHOCK_FLATTEN      2.5

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING      NO
#define  WARNING_MESSAGES       YES
#define  PRINT_TO_FILE          YES
#define  INTERNAL_BOUNDARY      NO
#define  SHOCK_FLATTENING       MULTID
#define  ARTIFICIAL_VISCOSITY   YES
#define  CHAR_LIMITING          YES
#define  LIMITER                MC_LIM
