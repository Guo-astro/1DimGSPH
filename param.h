#pragma once

namespace PARAM {
const double Dim = 3.0; //number of dimension
const double SMTH = 1.2;
//et
const double MU = 2.53;
const int OUTPUT_INTERVAL = 1;
const double C_CFL = 0.5; //CFL number
//const double Csmooth = 1.2; //Csmooth for calculation of smoothing length
const int SPACE_ORDER = 2; //space order for riemann solver
const double ACC = 1.0e-3; //infinitesimal constant
const double GAMMA = 5. / 3.; //specific heat ratio
const double G = .25 / M_PI;

//////////////////////////////////////////////////
///////////////Params CGS Units//////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
const double M_SUN_cgs = 1.988e33;
const double pc_cgs = 3.08567758e18;
const double accretion_alpha = .5;
const double G_cgs = 6.67384e-8;

const double PROTONMASS = 1.672621777e-24;
const double KBOLTZ_cgs = 1.3806488e-16;
const double yr = 3.1556952e7;
///////////////////////////////////////
const int NSIZE = 2;
const double poly_n = 1.5;
const double poly_gamma = 1.0 + 1.0 / poly_n;
const double poly_lambda = 1.0;
const double poly_K = 0.4;
const double poly_alpha = pow(pow(poly_lambda, -1.0 + 1.0 / poly_n) * (poly_n + 1.0) * poly_K / (4.0 * M_PI * G), .5);
const double xi = 3.65375374;
const double dxi = -0.2033013;
//const double xi = M_PI;
//const double dxi = -1./M_PI;
const double ddxi = -xi * xi * dxi;
const double poly_r = poly_alpha * xi;
const double poly_rhoc = poly_lambda;
const double poly_pc = poly_K * pow(poly_rhoc, 1.0 + 1. / poly_n);
const double poly_M = ddxi * 4.0 * M_PI * pow((poly_n + 1.0) * poly_K / (4.0 * M_PI * G), 1.5) * pow(poly_lambda, (3.0 - poly_n) / (2.0 * poly_n));
const double poly_E = (3. / (5. - poly_n)) * G * poly_M / poly_r;

//////////////////////////////////////////////////
///////////////Fundamental Simulation Units///////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

//const double SM = 30 * M_SUN_cgs;
const double SM = 3.00771 * M_SUN_cgs;
const double SL = 10 * pc_cgs;
const double ST = sqrt(pow(SL, 3.0) / SM / (4 * M_PI * G_cgs));
//////////////////////////////////////////////////
/////////////Combination SUnits////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
const double SPres = SM / (SL * ST * ST);
const double SVel = SL / ST;
const double SMassDens = SM / (SL * SL * SL);
const double SNumDens = 1. / (SL * SL * SL);
const double SEng = SM * SL * SL / (ST * ST);
const double SEng_per_Mass = SL * SL / (ST * ST);
const double SEng_dot_per_Mass = SL * SL / (ST * ST * ST);
const double SMyr = ST / yr / 1e6;

const double GSPH_SEng = SEng_per_Mass * SM;
const double GSPH_SEng_dot = SEng_dot_per_Mass * SMassDens;
//////////////////////////////////////////////////
/////////////Physical constants SUnits////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
const double SKB = KBOLTZ_cgs / SEng;
//////////////////////////////////////////////////
///////////////HVCC params CGSUnits///////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

/////////////Simulation Settings///////////////////////////////////////////
const double M_BH_cgs = 4.0e6 * M_SUN_cgs;
const double sM_BH = 1.e5 * M_SUN_cgs / SM;
const double D_GMC_BH_cgs = 10.0 * pc_cgs;
const double sD_GMC_BH = 10.0 * pc_cgs / SL;
const double sInit_Vel = 1e6 / SVel;
//const double Init_Vel_cgs = sqrt(2.0 * G_cgs * M_BH_cgs / D_GMC_BH_cgs);
//const double sInit_Vel = sqrt(2.0 * G_cgs * M_BH_cgs / D_GMC_BH_cgs) / SVel;
//const double sInit_Vel =10 ;
const double temprature_init = 20;

///////////////////////////////////////////////////
//Chemistry///////
const double tau = 3000.0;
const double OVERFLOW_LIMIT = 1.e-29;
const double PROTONMASS_CGS = 1.672621777e-24;
const double KBOLTZ_CGS = 1.3806488e-16;
const double PC_CGS = 3.08567758e18;
const double EV_CGS = 1.602176565e-12;
const double CR_ION_RATE = 1.8e-17;
const double G0 = 1.7;
const double CMB_T = 2.7;
const double YR_CGS = 3.1556952e7;
const double opratio = 2.4;
const double fortho = opratio / (1. + opratio);
const double fpara = 1. - fortho;
const double col_dens = 1.e20;

const int COMP = 11;

//3
const double i_abundance_OI = 4.6e-4;
//4
const double i_abundance_CI = 3.0e-4;
//2
const double i_abundance_HeI = 0.1;
//5
const double i_abundance_H2 = 1e-6;
//7
const double i_abundance_HII = 1e-5;
//9
const double i_abundance_Fe = 7.08e-7;
//10
const double i_abundance_Si = 3.55e-6;
//1
const double i_abundance_HI = fmax(1. - 2. * i_abundance_H2 - i_abundance_HII, 0.);
//0
const double i_abundance_e = 0.03;
//6
const double i_abundance_CO = 2.6e-4;
const double i_abundance_CII = i_abundance_CI - i_abundance_CO;

const double Grain_T = 8;
const double dust_to_gas_ratio = 1.0;
const double grain_size = 100.0;
}

