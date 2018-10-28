/*
 * laneEmden.h
 *
 *  Created on: Oct 28, 2018
 *      Author: guo
 */

#ifndef LANEEMDEN_H_
#define LANEEMDEN_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double getPoly53Phi(double x) {
	double a = 9.999995857999977E-01;
	double b = -3.215355146247138E-04;
	double c = -1.667864474453594E-01;
	double d = 5.252845006191754E-04;
	double e = 1.155630978128144E-02;
	double f = 1.050795977407310E-03;
	double g = -1.389661348847463E-03;
	double h = 2.648078982351008E-04;
	double i = -1.692465157938468E-05;
	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7) + i * pow(x, 8);

	return phi;

}
double getPoly53Phi_dash(double x) {
	double a = -3.257098683322130E-04;
	double b = -3.335390631705378E-01;
	double c = 1.494801635423298E-03;
	double d = 4.629812975807299E-02;
	double e = 5.242083986884144E-03;
	double f = -8.358773124387748E-03;
	double g = 1.867615193618750E-03;
	double h = -1.388085203319444E-04;
	double i = 2.998565502661576E-07;
	double phi = a + b * pow(x, 1) + c * pow(x, 2) + d * pow(x, 3) + e * pow(x, 4) + f * pow(x, 5) + g * pow(x, 6) + h * pow(x, 7) + i * pow(x, 8);

	return phi;

}
double getPoly53Dens(double x) {
	double res = pow(getPoly53Phi(x), 1.5);
	return res;

}

#endif /* LANEEMDEN_H_ */
