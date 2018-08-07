#pragma once
#include "ps_defs.cpp"
struct kernel_t {
	kernel_t() {
	}
	//W
	F64 W(const F64vec dr, const F64 h) const {
		F64 value;
//		if ((dr * dr) / (h * h) < 9.0) {
		value = pow(h * sqrt(M_PI), -1.) * exp(-(dr * dr) / (h * h));
//		} else {
//					value = 0.0;
//				}
		return value;
	}
	//gradW
	F64vec gradW(const F64vec dr, const F64 h) const {
		F64vec vec_value;
		//gradient of gaussian kernel
		if (sqrt(dr * dr) < 3.0) {
		vec_value = dr
				* (-2.0*pow(h * sqrt(M_PI), -1.) * (1.0/ (h * h))
						* exp(-(dr * dr) / (h * h)));
		} else {
			vec_value = 0.0;
		}
		return vec_value;
	}
	F64 DW_h(const F64vec dr, const F64 h) const {
		F64 value;
		//gradient of gaussian kernel
		value = (- h * h + 2 * dr * dr/(h*h)) * exp(-(dr * dr) / (h * h))
				/ ( pow(M_PI, 0.5));
		return value;

	}
	static F64 supportRadius() {
		return 4.0;
	}
};
