/*
 * heatingcooling.h
 *
 *  Created on: 2018/08/05
 *      Author: guo
 */

#ifndef APPS_KI2000_C10_HEATINGCOOLING_H_
#define APPS_KI2000_C10_HEATINGCOOLING_H_
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "param.h"
double fitLambda(double T);
double fitGamma(double T);
double fitTEMP(double numdens);
double fitAB_E(double numdens);
double fitAB_H2(double numdens);
double fitAB_CO(double numdens);
double fitAB_HI(double numdens);
// erg s-1///
double Gamma_PE(double T, double NUMDENS_CGS, double abundance_e);
// erg s-1///
double Lambda_PE(double T, double abundance_e, double NUMDENS_CGS);
// erg s-1///
double Gamma_CR(double NUMDENS_CGS, double abundance_e);
// erg s-1///
double Gamma_XR(double abundance_e, double col_dens_nutural);
// erg s-1 dissociating H2///
double Gamma_UV(double T, double NUMDENS_CGS, double x_2, double abundance_H, double col_den_H, double turb_vel);
// erg s-1 n n///
double Gamma_H2_grain_forming(double T, double T_gr, double NUMDENS_CGS, double x_2, double abundance_HI);

// s-1 n n///
double RATE_H2_UV_2H(double T, double col_den_H, double x_2, double turb_vel);
// cm3 s-1,  n_H,n///
double RATE_HHGrain_H2(double T, double T_gr);
// s-1///
double _D_CO();
//cm-3 s-1///
double _C_CO(double NUMDENS_CGS, double abundance_OI, double abundance_CII, double abundance_H2);
// s-1///
double RATE_XR_HII_HeII(double T, double abundance_e, double abundance_HI, double col_dens);
// s-1///
double RATE_CR_HII_HeII(double abundance_e);
// cm3 s-1///
double RATE_HHm_H2_e(double T, double abundance_HI, double abundance_HII);
// cm3 s-1///
double RATE_H_col_HII(double T);
// cm3 s-1///
double RATE_HIIe_HPhoton(double T);
// s-1///
double RATE_HH2_3H(double T);
// s-1///
double RATE_H2_CR_2H(double abundance_e);
// erg  cm-3 s-1///
double Lambda_Lya(double T, double abundance_e, double abundance_HI, double NUMDENS_CGS);
// erg cm-3 s-1///
double Lambda_OI(double T, double NUMDENS_CGS, double abundance_H, double abundance_H2, double abundance_e, double abundance_Hp, double abundance_OI, double tau,
		double col_dens_H2);
// erg cm-3 s-1 ///
double Lambda_CII(double T, double NUMDENS_CGS, double abundance_H, double abundance_H2, double abundance_e, double abundance_CII, double abundance_C, double tau_init,
		double col_dens_H2, double fpara, double fortho);
// erg cm-3 s-1 ///
double Lambda_SiII(double T, double NUMDENS_CGS, double abundance_H, double abundance_H2, double abundance_e, double abundance_Sip, double tau);
// erg cm-3 s-1 ///
double Lambda_FeII(double T, double NUMDENS_CGS, double abundance_HI, double abundance_H2, double abundance_e, double abundance_HII, double abundance_FeII);
// erg cm-3 s-1 ///
double Lambda_CO_rot(double T, double NUMDENS_CGS, double abundance_CO);
// erg cm-3 s-1 ///
double Lambda_Col_Dust(double T, double NUMDENS_CGS, double grain_size, double T_gr);

double rate_Lambda(double T, double NUMDENS_CGS, double abundance_HI, double abundance_e, double abundance_H2, double abundance_HII, double abundance_CII, double abundance_FeII,
		double abundance_SiII, double abundance_CI, double abundance_OI, double abundance_CO, double dust_to_gas_ratio, double T_gr, double grain_size);
double rate_Gamma(double T, double NUMDENS_CGS, double abundance_HI, double abundance_e, double abundance_H2, double abundance_HII, double abundance_CII, double dust_to_gas_ratio,
		double T_gr, double col_den_H2, double GAMMA, double mu);
#endif /* APPS_KI2000_C10_HEATINGCOOLING_H_ */
