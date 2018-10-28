#include "arraylist.cpp"
#include "relocarray.cpp"
#include "radix_sort.cpp"
#include "ps_defs.cpp"
#include "kernel.h"
#include "param.h"
#include "heatingcooling.h"
#include "laneEmden.h"
#include <cmath>
using namespace std;
double fRand_Double(double fMin, double fMax) {
	double f = (double) rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}
double GAMMA = 5. / 3.;
int TYPE_FLUID = 0;
int TYPE_GHOST = 1;
F64 CSMOOTH = 2.0;
struct Particle {
	F64 mu;
	F64 temp;
	F64 NUMDENS;
	F64vec pos;
	F64vec vel;
	F64vec acc;
	F64vec gradV;
	F64vec graddens;
	F64vec gradpres;
	F64vec gradvel_x;
	F64vec gradvel_y;
	F64vec gradvel_z;
	double dens;
	double snds;
	double pres;
	double smth;
	double mass;
	double eng;
	double dt;
	double OMEGA;
	F64vec vel_half;
	double eng_half;
	double eng_dot;
	int TYPE;
	int id;
	double cooling_timescale;
	double Gamma;
	double Lambda;

	const char* scalars[PARAM::COMP + 11] = { "mass", "pres", "dens", "vx", "vy", "vz", "T", "eng", "ab_e", "ab_HI", "ab_HeI", "ab_OI", "ab_CI", "ab_H2", "ab_CO", "ab_HII",
			"ab_CII", "ab_FeII", "ab_SiII", "cooling", "heating", "cooling_timescale" };
	const char* vectors[3] = { "vel", "acc", "pos" };

	double abundances[PARAM::COMP];
	double old_abundances[PARAM::COMP];
	U64 getKey() {
		return id;
	}
	bool operator==(const Particle& other) const {
		return id == other.id;
	}

//	Particle & operator =(const Particle& other) {
//
//		pos = other.pos     ;
//		vel = other.vel     ;
//		acc = other.acc     ;
//		dens = other.dens   ;
//		snds = other.snds   ;
//		pres = other.pres   ;
//		smth = other.smth   ;
//		mass = other.mass   ;
//		eng = other.eng;
//		dt = other.dt;
//		vel_half = other.vel_half;
//		eng_half = other.eng_half;
//		eng_dot = other.eng_dot;
//		TYPE = other.TYPE;
//		id = other.id;
//		return *this;
//	}
};
typedef struct {
	F64vec max;
	F64vec min;
	F64vec domain_len;
} Domain;
typedef RellocatableArray<Particle> Particles;
void evolve_abundance(Particle &hydro, const double oneDynTimeStep);

void find_nnbs(const Particle &pi, const Particles &ps, RellocatableArray<Particle> &nbrs, long nparts) {
	RellocatableArray<int> ptcl_id;

//	RadixSort<unsigned long int, 8> rs_;
//	rs_.lsdSort(ps.getPointer(),buf_.getPointer(),0,nparts);
	for (int i = 0; i < ps.size(); i++) {

		double rij = sqrt((pi.pos - ps[i].pos) * (pi.pos - ps[i].pos));
		if (rij < 6. * 2 * .5 * (pi.smth + ps[i].smth)) {

			nbrs.push_back(ps[i]);

		}

	}

}
void calc_density(Particles &ps, long nparts) {
	kernel_t ker;
	double hmin = 1e30;
	for (int i = 0; i < ps.size(); i++) {
		if (ps[i].id > 1000) {
			continue;
		}
		double tmp_dens = 0.0;
		ps[i].dens = 0.0;
		ps[i].OMEGA = 0.0;
		RellocatableArray<Particle> nbrs;
		find_nnbs(ps[i], ps, nbrs, nparts);

		for (int j = 0; j < nbrs.size(); j++) {
			F64vec rij = ps[i].pos - nbrs[j].pos;
			ps[i].dens += nbrs[j].mass * ker.W(rij, ps[i].smth);
			if (ps[i].dens != ps[i].dens || ps[i].smth != ps[i].smth || nbrs[j].smth != nbrs[j].smth) {
				cout << "Nan in dens " << nbrs[j].mass << " " << ps[i].dens << " " << rij.x << " " << ps[i].smth << " " << ps[i].id << " " << nbrs[j].dens << " " << nbrs[j].smth
						<< " " << ker.W(rij, ps[i].smth) << endl;
				while (true) {
				}
			}

//			if (ps[i].dens > 5) {
//				cout << "too small smnth in dens " << " mass " << nbrs[j].mass << " dens " << ps[i].dens << " rijx " << rij.x << " pos " << ps[i].pos.x << " smth " << ps[i].smth
//						<< " id " << ps[i].id << " densj " << nbrs[j].dens << " smthj " << nbrs[j].smth << " idj " << nbrs[j].id << " ker " << ker.W(rij, ps[i].smth) << endl;
////				ps[i].smth = 1e-4;
//
////				while (true) {
////				}
//			}
		}

		ps[i].smth = ps[i].mass / ps[i].dens;
	}
}
void calc_presure(Particles &ps, long nparts, double &glb_dt) {
	kernel_t ker;
//	nparts = ps.size();
	glb_dt = 1e30;
	for (int i = 0; i < ps.size(); i++) {
		if (ps[i].id > 1000) {
			continue;
		}
		ps[i].pres = .5 * pow(ps[i].dens, 2);
		ps[i].snds = sqrt(GAMMA * ps[i].pres / ps[i].dens);
		ps[i].dt = 0.4 * ps[i].smth / ps[i].snds;
		double vel_crit = 0.3 * fabs(ps[i].pos.x) / sqrt(ps[i].vel * ps[i].vel);
		glb_dt = fmin(glb_dt, ps[i].dt);
		glb_dt = fmin(glb_dt, vel_crit);
//
//		}

	}
}
//TODO higher dimension notice
void copyGhosts(const Domain &d, Particles &ps, long nparts) {
	double srch = 18;

	ps.remove_ghost();
	for (int i = 0; i < ps.size(); i++) {

		double partsx = ps[i].pos.x;
		double _minx = 0.0;
		double diff = partsx - _minx;
		if (diff < srch * ps[i].smth && diff > 0.0) {
			Particle ghost;
			ghost.pos.x = _minx - diff;

			if (ps[i].pos.x == 0 || ps[i].pos.x < 1e-5) {
				cout << "cp ghost ptcl in < boundary domain!!" << ps[i].pos.x << endl;

//				ghost.pos.x = d.min.x - diff - 1e-4;
				cout << ps[i].eng << endl;

			}

			ghost.acc = ps[i].acc;
			ghost.mass = ps[i].mass;
			ghost.dens = ps[i].dens;
			ghost.pres = ps[i].pres;
			ghost.eng = ps[i].eng;
			ghost.smth = ps[i].smth;
			ghost.vel = -ps[i].vel;
			ghost.id = ps[i].id + nparts + 1000;
			ghost.temp = ps[i].temp;
			ghost.NUMDENS = ps[i].NUMDENS;
//			cout << ghost.id << endl;

			ghost.TYPE = TYPE_GHOST;
			ps.push_back(ghost);

		}
	}
//	cout << ps.size() << endl;
}
void calc_gradients(Particles &ps, long nparts) {
	kernel_t kernel;

	for (S32 i = 0; i < nparts; ++i) {
		ps[i].gradV = 0.0;
		ps[i].graddens = 0.0;
		ps[i].gradpres = 0.0;
		ps[i].gradvel_x = 0.0;
		ps[i].gradvel_y = 0.0;
		ps[i].gradvel_z = 0.0;
//		nparts = ps.size();
		RellocatableArray<Particle> nbrs;
		find_nnbs(ps[i], ps, nbrs, nparts);
		const F64 rhoi_inv2 = 1.0 / (ps[i].dens * ps[i].dens);
		if (ps[i].pos.x < 0.4 && ps[i].pos.x > 0.39) {
			const F64 rhoi_inv2 = 1.0 / (ps[i].dens * ps[i].dens);
		}
		for (S32 j = 0; j < nbrs.size(); ++j) {
			if (ps[i].pos != nbrs[j].pos) {

				const F64vec dr = ps[i].pos - nbrs[j].pos;
				const F64vec eij = dr / sqrt(dr * dr);
				const F64vec dv = ps[i].vel - nbrs[j].vel;

				const F64vec gradW = kernel.gradW(dr, ps[i].smth);
				ps[i].gradV += -nbrs[j].mass * gradW;
				ps[i].graddens += nbrs[j].mass * gradW;
				ps[i].gradpres += (nbrs[j].mass / nbrs[j].dens) * (nbrs[j].pres - ps[i].pres) * gradW;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
				ps[i].gradvel_x += (nbrs[j].vel.x- ps[i].vel.x) * gradW * nbrs[j].mass / nbrs[j].dens;
				ps[i].gradvel_y += (nbrs[j].vel.y - ps[i].vel.y) * gradW * nbrs[j].mass / nbrs[j].dens;
#else
				ps[i].gradvel_x += (nbrs[j].vel.x - ps[i].vel.x) * gradW * nbrs[j].mass / nbrs[j].dens;
//				ps[i].gradvel_y += (nbrs[j].vel.y - ps[i].vel.y) * gradW
//						* nbrs[j].mass / nbrs[j].dens;
//				ps[i].gradvel_z += (nbrs[j].vel.z - ps[i].vel.z) * gradW
//						* nbrs[j].mass / nbrs[j].dens;
#endif
//					std::cout<<ps[i].gradvel_z<<std::endl;
			}

		}
		ps[i].gradV *= rhoi_inv2;
	}
}
void calc_force(Particles &ps, long nparts, double &glb_dt) {
	kernel_t kernel;

	for (int i = 0; i < nparts; i++) {
		if (ps[i].id > 1000) {
			continue;
		}
		ps[i].acc = 0.0;
		ps[i].eng_dot = 0.0;
		F64 v_sig_max = 0.0;
		RellocatableArray<Particle> nbrs;
		find_nnbs(ps[i], ps, nbrs, nparts);

		for (S32 j = 0; j < nbrs.size(); ++j) {
			const F64vec dr = ps[i].pos - nbrs[j].pos;
			const F64vec dv = ps[i].vel - nbrs[j].vel;
			const F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
			const F64 v_sig = ps[i].snds + nbrs[j].snds - 3.0 * w_ij;
			v_sig_max = std::max(v_sig_max, v_sig);
			const F64 AV = -0.5 * v_sig * w_ij / (0.5 * (ps[i].dens + nbrs[j].dens));

			const F64vec gradW = 0.5 * (kernel.gradW(dr, ps[i].smth) + kernel.gradW(dr, nbrs[j].smth));
			ps[i].acc -= nbrs[j].mass * (ps[i].pres / (ps[i].dens * ps[i].dens) + nbrs[j].pres / (nbrs[j].dens * nbrs[j].dens) + AV) * gradW;
			ps[i].eng_dot += nbrs[j].mass * (ps[i].pres / (ps[i].dens * ps[i].dens) + 0.5 * AV) * dv * gradW;
		}
		ps[i].dt = 0.5 * ps[i].smth / (ps[i].snds);
		glb_dt = fmin(glb_dt, ps[i].dt);
	}
}

void calc_Vij2_and_ss(const Particle ep_i, const Particle ep_j, F64&Vij2_hi, F64&Vij2_hj, F64& ssP, const F64 delta, const F64vec eij) {
	const F64 Vi = 1.0 / ep_i.dens;
	const F64 Vj = 1.0 / ep_j.dens;
	const F64 hi2 = ep_i.smth * ep_i.smth;
	const F64 hj2 = ep_j.smth * ep_j.smth;
//		const F64 rhoij = 0.5 * (ps.dens + ep_j.dens);
//		const F64 Vij2crit = 10.0 * (.0 / (rhoij * rhoij));
	const F64 Pij = 0.5 * (ep_i.pres + ep_j.pres);
	F64 Aij, Bij, Cij, Dij;

	Cij = (Vi - Vj) / delta;
	Dij = 0.5 * (Vi + Vj);
	Vij2_hi = 0.25 * hi2 * Cij * Cij + Dij * Dij;
	Vij2_hj = 0.25 * hj2 * Cij * Cij + Dij * Dij;
	ssP = (hi2 * Cij * Dij / (2 * Vij2_hi));
//		Aij = -2 * (Vi - Vj) * pow(delta, -3) + (dVi + dVj) * pow(delta, -2);
//		Bij = 0.5 * (dVi - dVj) / delta;
//		Cij = 1.5 * (Vi - Vj) / delta - 0.25 * (dVi + dVj);
//		Dij = 0.5 * (Vi + Vj) - 0.125 * (dVi - dVj) * delta;
//		Vij2_hi = ((0.234375 * hi2 * Aij * Aij + 0.1875 * (2 * Aij * Cij + Bij * Bij)) * hi2 + 0.25 * (2 * Bij * Dij + Cij * Cij)) * hi2 + Dij * Dij;
//		Vij2_hj = ((0.234375 * hj2 * Aij * Aij + 0.1875 * (2 * Aij * Cij + Bij * Bij)) * hj2 + 0.25 * (2 * Bij * Dij + Cij * Cij)) * hj2 + Dij * Dij;
//		ssP = 0.5 * (((0.46875 * hi2 * Aij * Bij + 0.375 * (Aij * Dij + Bij * Cij)) * hi2 + 0.5 * Cij * Dij) * hi2 / Vij2_hi + ((0.46875 * hj2 * Aij * Bij + 0.375 * (Aij * Dij + Bij * Cij)) * hj2 + 0.5 * Cij * Dij) * hj2 / Vij2_hj);
}

void calc_riemann_solver(const Particle &ep_i, const Particle &ep_j, const F64 &ss, const F64 delta, const F64vec &eij, const F64 dt, F64 &pstar, F64 &vstar) {

	F64 rho1, rho2, v1, v2, p1, p2;
	F64 vi = ep_i.vel * eij;
	F64 vj = ep_j.vel * eij;

	rho1 = ep_i.dens;
	p1 = ep_i.pres;
	v1 = vi;
	rho2 = ep_j.dens;
	p2 = ep_j.pres;
	v2 = vj;

	F64 ppre, p, v;
	F64 W1, W2;
	const F64 alpha = (2.0 * GAMMA) / (GAMMA - 1.0);
	p = .5 * (p1 + p2);
	double critif = 1.0 - 1.0e-6;
	if (p1 <= 0.0 || p2 <= 0.0) {
		pstar = .5 * (p1 + p2);
		vstar = .5 * (v1 + v2);
	} else {
		for (U32 loop = 0; loop < 50; loop++) {
			ppre = p;
			if ((p / p1 < critif))
				W1 = sqrt(p1 * rho1) * ((GAMMA - 1.0) / (2.0 * sqrt(GAMMA))) * (1.0 - p / p1) / (1.0 - pow(p / p1, 1.0 / alpha));
			else
				W1 = sqrt(p1 * rho1) * sqrt(0.5 * (GAMMA + 1.0) * p / p1 + 0.5 * (GAMMA - 1.0));
			if ((p / p2 < critif))
				W2 = sqrt(p2 * rho2) * ((GAMMA - 1.0) / (2.0 * sqrt(GAMMA))) * (1.0 - p / p2) / (1.0 - pow(p / p2, 1.0 / alpha));
			else
				W2 = sqrt(p2 * rho2) * sqrt(0.5 * (GAMMA + 1.0) * p / p2 + 0.5 * (GAMMA - 1.0));
			p = ((p2 / W2 + p1 / W1) + v2 - v1) / (1.0 / W2 + 1.0 / W1);
			if (p < 0.0)
				p = 0.5 * ppre;
			if (fabs(p - ppre) < 1e-3) {
//				std::cout <<W1 << "==" << W2<< "=="<< v2 <<" "<<v1<<" "<<p2<<" "<<p1<<" "<< p << std::endl;

				break;
			}
		}
		vstar = ((W1 * v1 + W2 * v2) + p2 - p1) / (W1 + W2);
		pstar = p;
		if ((pstar != pstar) || vstar != vstar || pstar <= 0.0) {
			vstar = .5 * (v1 + v2);
			pstar = .5 * (p1 + p2);
			std::cout << p1 << "==" << p2 << "==" << pstar << std::endl;
		}

	}
}
struct interP {
	F64 VIJ_I;
	F64 VIJ_J;
	F64 sstar;
	F64 PSTAR;
	F64 VSTAR;

};
void calc_force_G(Particles &ps, long nparts, const double &glb_dt) {
	kernel_t kernel;

//	calc_gradients(ps, nparts);
	for (S32 i = 0; i < nparts; ++i) {
		if (ps[i].id > 1000) {
			continue;
		}
		ps[i].acc = 0.0;
		ps[i].eng_dot = 0.0;
//		nparts = ps.size();
		RellocatableArray<Particle> nbrs;
		find_nnbs(ps[i], ps, nbrs, nparts);
		RellocatableArray<interP> intervals;

		for (S32 j = 0; j < nbrs.size(); ++j) {
			F64 ss, Vij2_hi, Vij2_hj;

			if (ps[i].pos != nbrs[j].pos) {

				F64 VIJ_I, VIJ_J, PSTAR, VSTAR;
				F64 sstar = 0.0;
				F64vec interpo = 0.0;
				F64vec eij;
				F64vec gradW_hsi;
				F64vec gradW_hsj;
				const F64vec dv = ps[i].vel - nbrs[j].vel;
				const F64vec dr = ps[i].pos - nbrs[j].pos;
				const F64 delta = sqrt(dr * dr);

				eij = dr / delta;
				const F64 si = ps[i].pos * eij;
				const F64 sj = nbrs[j].pos * eij;
				const F64 dsij = si - sj;
				F64vec vij = 0.0;
				gradW_hsi = kernel.gradW(dr, ps[i].smth);

				ps[i].acc -= nbrs[j].mass * gradW_hsi;

			}
		}
		ps[i].acc +=   getPoly53Phi_dash(sqrt(ps[i].pos.x * ps[i].pos.x)) * ps[i].pos.x / sqrt(ps[i].pos.x * ps[i].pos.x);
	}

}

void relaxation(Particles &ps, const vector<Particles> &res, long nparts, double &glb_dt) {
	kernel_t kernel;

	for (int i = 0; i < nparts; i++) {
		ps[i].acc = 0.0;
		ps[i].eng_dot = 0.0;
		F64 v_sig_max = 0.0;
//			cout<<res[i].size()<<endl;
		for (S32 j = 0; j < res[i].size(); ++j) {
			const F64vec dr = ps[i].pos - ps[j].pos;
			const F64vec dv = ps[i].vel - ps[j].vel;
			const F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;
			const F64 v_sig = ps[i].snds + ps[j].snds - 3.0 * w_ij;
			v_sig_max = std::max(v_sig_max, v_sig);
			const F64 AV = -0.5 * v_sig * w_ij / (0.5 * (ps[i].dens + ps[j].dens));

			const F64vec gradW = 0.5 * (kernel.gradW(dr, ps[i].smth) + kernel.gradW(dr, ps[j].smth));
			ps[i].acc -= ps[j].mass * gradW;

		}
		ps[i].acc += -sin(ps[i].pos.x);
		ps[i].dt = .5 * ps[i].smth / ps[i].snds;
		glb_dt = fmin(glb_dt, ps[i].dt);
	}
}
void leapfrog_step(Particles &ps, double dt, long nparts) {
	for (int i = 0; i < nparts; i++) {
		ps[i].vel_half += .5 * ps[i].acc * dt;
		ps[i].pos += ps[i].vel_half * dt;

	}

}
void leapfrog_first_step(Particles &ps, double dt, long nparts) {
	for (int i = 0; i < nparts; i++) {
		ps[i].vel_half += .5 * ps[i].acc * dt;
		ps[i].pos += ps[i].vel_half * dt;

	}

}
void relax_step(Particles &ps, double dt, long nparts) {
	for (int i = 0; i < nparts; i++) {
		ps[i].vel_half = ps[i].vel = 0.0;
		ps[i].pos += .5 * ps[i].acc * dt * dt;

	}

}
void relax_first_step(Particles &ps, double dt, long nparts) {
	for (int i = 0; i < nparts; i++) {
		ps[i].vel_half = ps[i].vel = 0.0;
		ps[i].pos += .5 * ps[i].acc * dt * dt;

	}

}

void InitialKick(Particles &ps, double dt, long nparts) {
	for (S32 i = 0; i < nparts; ++i) {
		ps[i].vel_half = ps[i].vel + 0.5 * dt * ps[i].acc;
		ps[i].eng_half = ps[i].eng + 0.5 * dt * ps[i].eng_dot;
	}
}

void FullDrift(Particles &ps, double dt, long nparts) {
//time becomes t + dt;
	for (S32 i = 0; i < nparts; ++i) {
		ps[i].pos += dt * ps[i].vel_half;
	}
}

void Predict(Particles &ps, double dt, long nparts) {
	for (S32 i = 0; i < nparts; ++i) {
		ps[i].vel += dt * ps[i].acc;
		ps[i].eng += dt * ps[i].eng_dot;
	}
}

void FinalKick(Particles &ps, double dt, long nparts) {
//	nparts = ps.size();
	for (S32 i = 0; i < ps.size(); ++i) {
		if (ps[i].id > 1000) {
			continue;
		}
		ps[i].pos += .5 * ps[i].acc * dt * dt;
	}

}

void evolve_abundance(Particle &hydro, const double oneDynTimeStep) {

	double abundance_e = hydro.abundances[0];
	double abundance_HI = hydro.abundances[1];
	double abundance_HeI = hydro.abundances[2];
	double abundance_OI = hydro.abundances[3];
	double abundance_CI = hydro.abundances[4];
	double abundance_H2 = hydro.abundances[5];
	double abundance_CO = hydro.abundances[6];
	double abundance_HII = hydro.abundances[7];
	double abundance_CII = hydro.abundances[8];
	double abundance_FeII = hydro.abundances[9];
	double abundance_SiII = hydro.abundances[10];
	double abundance_e_old = abundance_e;
	double abundance_HI_old = abundance_HI;
	double abundance_HeI_old = abundance_HeI;
	double abundance_OI_old = abundance_OI;
	double abundance_CI_old = abundance_CI;
	double abundance_H2_old = abundance_H2;
	double abundance_CO_old = abundance_CO;
	double abundance_HII_old = abundance_HII;
	double abundance_CII_old = abundance_CII;
	double abundance_FeII_old = abundance_FeII;
	double abundance_SiII_old = abundance_SiII;
	double mu = hydro.mu;
	double dens = hydro.dens;

	double NUMDENS_CGS = hydro.NUMDENS;
	double MASSDENS_CGS = hydro.dens * PARAM::SMassDens;

	double tdust = PARAM::Grain_T;
	double dust_to_gas_ratio = PARAM::dust_to_gas_ratio;

	double energy_old = hydro.eng * PARAM::SEng_per_Mass;

	double energy_init = hydro.eng * PARAM::SEng_per_Mass;

	double abundance_e_init = abundance_e;
	double abundance_HI_init = abundance_HI;
	double abundance_HeI_init = abundance_HeI;
	double abundance_OI_init = abundance_OI;
	double abundance_CI_init = abundance_CI;
	double abundance_H2_init = abundance_H2;
	double abundance_CO_init = abundance_CO;
	double abundance_HII_init = abundance_HII;
	double abundance_CII_init = abundance_CII;
	double abundance_FeII_init = abundance_FeII;
	double abundance_SiII_init = abundance_SiII;
	double GAMMA_GAS = 5. / 3.;
	double temprature_init = hydro.mu * PARAM::PROTONMASS_CGS * (hydro.pres) * PARAM::SEng_per_Mass / (hydro.dens * PARAM::KBOLTZ_CGS * (1.1 + abundance_e - abundance_H2));

	double ydot_H2 = 0.0;
	double ydot_Hp = 0.0;
	double ydot_CO = 0.0;

	double Gamma, Lambda;
	double ylam = 0.0;

	double energy = 0.0;
	int main_try_steps = 1e8, sub_try_steps = 1000, chem_try_steps = 100;
	double t_start = 0.;

	double fac = 0.1;
	double time_passed = 0.0, tleft = 0.0;
	double dt = oneDynTimeStep * PARAM::ST;
	double dt_substep = dt;
	int num_converged = 0;
	double abs_tol_H = 1e-7, rel_tol_H = fac, rel_diff_H;
	double abs_tol_HII = 1e-7, rel_tol_HII = fac, rel_diff_HII;
	double abs_tol_H2 = 1e-7, rel_tol_H2 = fac, rel_diff_H2;
	double abs_tol_CII = 1e-7, rel_tol_CII = fac, rel_diff_CII;
	double abs_tol_CO = 1e-7, rel_tol_CO = fac, rel_diff_CO;
	double abs_tol_e = 1e-7, rel_tol_e = fac, rel_diff_e;
	double abs_tol_eng, rel_tol_eng = fac, rel_diff_eng;
	double loop = 1000;
	double T = 0.0;
	for (int step = 0; step < main_try_steps; step++) {

//		std::cout << dt << " " << time_passed << " " << tleft << "out loop" << std::endl;
		int iter_success = 0;
		tleft = oneDynTimeStep * PARAM::ST - time_passed;
		dt = (tleft < dt) ? tleft : dt;
		for (int sub_step = 0; sub_step < sub_try_steps; sub_step++) {
			abundance_e = abundance_e_init;
			abundance_HI = abundance_HI_init;
			abundance_OI = abundance_OI_init;
			abundance_CI = abundance_CI_init;
			abundance_H2 = abundance_H2_init;
			abundance_CO = abundance_CO_init;
			abundance_HII = abundance_HII_init;
			abundance_CII = abundance_CII_init;
			abundance_FeII = abundance_FeII_init;
			abundance_SiII = abundance_SiII_init;
			T = temprature_init;
			energy = energy_init;
//
			double K1 = RATE_CR_HII_HeII(abundance_e);
			double K2 = RATE_XR_HII_HeII(T, abundance_e, abundance_HI, PARAM::col_dens);
			double K3 = RATE_H_col_HII(T);
			double K4 = RATE_HIIe_HPhoton(T);
			double K5 = RATE_HHGrain_H2(T, tdust);
			double K6 = RATE_HHm_H2_e(T, abundance_HI, abundance_HII);
			double K7 = RATE_H2_UV_2H(T, PARAM::col_dens, mu, abundance_H2);
			double K8 = RATE_HH2_3H(T);
			double K9 = RATE_H2_CR_2H(abundance_e);

			double C_HII = K1 * abundance_HI * NUMDENS_CGS + K2 * abundance_HI * NUMDENS_CGS + K3 * NUMDENS_CGS * abundance_e * NUMDENS_CGS * abundance_HI;
			double D_HII = K4 * abundance_e * NUMDENS_CGS;
			double C_HI = K4 * abundance_e * abundance_HII * NUMDENS_CGS * NUMDENS_CGS + K7 * abundance_H2 * NUMDENS_CGS + K8 * abundance_H2 * abundance_HI * NUMDENS_CGS
					+ K9 * abundance_H2 * NUMDENS_CGS;
			double D_HI = K1 + K2 + K3 * abundance_e * NUMDENS_CGS + K5 * abundance_HI * NUMDENS_CGS;
			double C_H2 = K5 * abundance_HI * abundance_HI * NUMDENS_CGS * NUMDENS_CGS + K6 * abundance_e * abundance_HI * NUMDENS_CGS * abundance_HI * NUMDENS_CGS;
			double D_H2 = K7 + K8 * NUMDENS_CGS * abundance_HI + K9;
			double C_CO = _C_CO(NUMDENS_CGS, abundance_OI, abundance_CII, abundance_H2);
			double D_CO = _D_CO();

			abundance_HII = C_HII / (NUMDENS_CGS) / D_HII + (abundance_HII_init - C_HII / (NUMDENS_CGS) / D_HII) * exp(-D_HII * .5 * dt);
			abundance_HI = C_HI / (NUMDENS_CGS) / D_HI + (abundance_HI_init - C_HI / (NUMDENS_CGS) / D_HI) * exp(-D_HI * .5 * dt);
			abundance_H2 = C_H2 / (NUMDENS_CGS) / D_H2 + (abundance_H2_init - C_H2 / (NUMDENS_CGS) / D_H2) * exp(-D_H2 * .5 * dt);
			abundance_CO = C_CO / (NUMDENS_CGS) / D_CO + (abundance_CO_init - C_CO / (NUMDENS_CGS) / D_CO) * exp(-D_CO * .5 * dt);
			abundance_HII = abundance_HII / (abundance_HII + abundance_HI + 2. * abundance_H2);
			abundance_HI = abundance_HI / (abundance_HII + abundance_HI + 2. * abundance_H2);
			abundance_H2 = abundance_H2 / (abundance_HII + abundance_HI + 2. * abundance_H2);
			abundance_CII = fmax(abundance_CI - abundance_CO, 0.0);
			abundance_e = fmin(abundance_HII + abundance_CII + abundance_SiII, 1.0);
//			std::cout << "temp: " << abundance_HII_init << " xe: " << K1<<" "<<K2<<" "<<K3 << " "<<K4<<" CHII "<<C_HII<<" DHII "<<D_HII<<" dt "<< dt <<std::endl;
//						std::cout << "xe: " << abundance_e  <<std::endl;

			GAMMA_GAS = PARAM::GAMMA;
			Gamma = rate_Gamma(T, NUMDENS_CGS, abundance_HI, abundance_e, abundance_H2, abundance_HII, abundance_CII, dust_to_gas_ratio, PARAM::Grain_T, PARAM::col_dens, GAMMA_GAS,
					mu);
			Lambda = rate_Lambda(T, NUMDENS_CGS, abundance_HI, abundance_e, abundance_H2, abundance_HII, abundance_CII, abundance_FeII, abundance_SiII, abundance_CI, abundance_OI,
					abundance_CO, dust_to_gas_ratio, PARAM::Grain_T, PARAM::grain_size);
			rel_diff_eng = (fabs(Gamma - Lambda) * dt) / (energy * MASSDENS_CGS);

			if (rel_diff_eng + 1.e-6 > rel_tol_eng) {
				iter_success = -1;
				break;
			}
			energy = (energy_init + (Gamma - Lambda) * dt / (MASSDENS_CGS));
			double mu = abundance_HI + abundance_HII + 2.0 * abundance_H2 + 4.0 * abundance_HeI;
			mu /= (abundance_HI + abundance_HII + abundance_H2 + abundance_HeI + abundance_e);

//			 mu * PARAM::PROTONMASS_CGS * (pres) * PARAM::SEng_per_Mass
//						/ (dens * PARAM::KBOLTZ_CGS * (1.1 + abundance_e - abundance_H2))
//			 mu * PARAM::PROTONMASS_CGS * (energy * (GAMMA_GAS - 1.0) ) * PARAM::SEng_per_Mass
//						/ ( PARAM::KBOLTZ_CGS * (1.1 + abundance_e - abundance_H2))
			T = mu * PARAM::PROTONMASS_CGS * energy * (GAMMA_GAS - 1.0) / (PARAM::KBOLTZ_CGS * (1.1 + abundance_e - abundance_H2));

			num_converged = 2;

			rel_diff_eng = (fabs(Gamma - Lambda) * dt) / (energy * MASSDENS_CGS);

//			std::cout << T << std::endl;
			if (rel_diff_eng + 1.e-6 < rel_tol_eng) {
				iter_success = num_converged;
				break;
			}

		}

		if (iter_success > 0) {
			double eng_eq = 0.0;
//						std::cout << T << "  "<<temprature_init<<std::endl;

			if (fabs(T - temprature_init) + 1e-6 < 1.e-3) {

				hydro.NUMDENS = NUMDENS_CGS;
				hydro.temp = T;
				hydro.mu = mu;
				hydro.Lambda = Lambda * 1e26;
				hydro.Gamma = Gamma * 1e26;
				hydro.cooling_timescale = ((energy / fabs(Lambda)) * MASSDENS_CGS) / PARAM::yr;
				hydro.eng = eng_eq = energy / PARAM::SEng_per_Mass;
				hydro.abundances[0] = abundance_e;
				hydro.abundances[1] = abundance_HI;
				hydro.abundances[2] = abundance_HeI;
				hydro.abundances[3] = abundance_OI;
				hydro.abundances[4] = abundance_CI;
				hydro.abundances[5] = abundance_H2;
				hydro.abundances[6] = abundance_CO;
				hydro.abundances[7] = abundance_HII;
				hydro.abundances[8] = abundance_CII;
				hydro.abundances[9] = abundance_FeII;
				hydro.abundances[10] = abundance_SiII;
//				std::cout << time_passed << " " << oneDynTimeStep * PARAM::ST
//									<< " " << hydro.cooling_timescale << std::endl;
				break;

			}
			abundance_e_init = abundance_e;
			abundance_HI_init = abundance_HI;
			abundance_OI_init = abundance_OI;
			abundance_CI_init = abundance_CI;
			abundance_H2_init = abundance_H2;
			abundance_CO_init = abundance_CO;
			abundance_HII_init = abundance_HII;
			abundance_CII_init = abundance_CII;
			abundance_FeII_init = abundance_FeII;
			abundance_SiII_init = abundance_SiII;
			temprature_init = T;
			energy_init = energy;
			time_passed += dt;

//			if (time_passed >= oneDynTimeStep * PARAM::ST) {
//				hydro.Lambda = Lambda * 1e26;
//				hydro.NUMDENS = NUMDENS_CGS;
//				hydro.temp = T;
//				hydro.mu = mu;
//				hydro.Gamma = Gamma * 1e26;
//				hydro.cooling_timescale = ((energy / fabs(Lambda)) * MASSDENS_CGS) / PARAM::yr;
//				hydro.eng = eng_eq = energy / PARAM::SEng_per_Mass;
//				hydro.abundances[0] = abundance_e;
//				hydro.abundances[1] = abundance_HI;
//				hydro.abundances[2] = abundance_HeI;
//				hydro.abundances[3] = abundance_OI;
//				hydro.abundances[4] = abundance_CI;
//				hydro.abundances[5] = abundance_H2;
//				hydro.abundances[6] = abundance_CO;
//				hydro.abundances[7] = abundance_HII;
//				hydro.abundances[8] = abundance_CII;
//				hydro.abundances[9] = abundance_FeII;
//				hydro.abundances[10] = abundance_SiII;
//				hydro.old_abundances[0] = abundance_e;
//				hydro.old_abundances[1] = abundance_HI;
//				hydro.old_abundances[2] = abundance_HeI;
//				hydro.old_abundances[3] = abundance_OI;
//				hydro.old_abundances[4] = abundance_CI;
//				hydro.old_abundances[5] = abundance_H2;
//				hydro.old_abundances[6] = abundance_CO;
//				hydro.old_abundances[7] = abundance_HII;
//				hydro.old_abundances[8] = abundance_CII;
//				hydro.old_abundances[9] = abundance_FeII;
//				hydro.old_abundances[10] = abundance_SiII;
//
////				std::cout << time_passed << " " << oneDynTimeStep * PARAM::ST
////									<< " " << hydro.cooling_timescale << std::endl;
//				break;
//			}
		} else {
			abundance_e = abundance_e_init;
			abundance_HI = abundance_HI_init;
			abundance_OI = abundance_OI_init;
			abundance_CI = abundance_CI_init;
			abundance_H2 = abundance_H2_init;
			abundance_CO = abundance_CO_init;
			abundance_HII = abundance_HII_init;
			abundance_CII = abundance_CII_init;
			abundance_FeII = abundance_FeII_init;
			abundance_SiII = abundance_SiII_init;
			T = temprature_init;
			energy = energy_init;

			//

			dt *= 0.8;
		}

	}
}
int main() {
	Particles ps;
#ifndef RESTART
	double xmin = -PARAM::xi + 0.1;
	double xmax = PARAM::xi - 0.1;
	double domain_len = xmax - xmin;
	int nparts = 0;
	int STEP_LIMIT = 0;
	int id = 0;
	STEP_LIMIT = 90;
	nparts = 200;
	double dx = (xmax - xmin) / nparts;
	double totMass = 0;
	std::cout << xmax << std::endl;

	for (F64 x = xmin; x < xmax; x += dx) {
		id++;
		Particle pi;
		pi.pos = x;
		pi.id = id;
		pi.vel = 0.0;
		pi.dens = getPoly53Dens(sqrt(x * x));
		totMass += pi.dens * dx;
		pi.pres = pow(pi.dens, PARAM::GAMMA);

		pi.eng = pi.pres / (pi.dens * (GAMMA - 1));
		pi.TYPE = TYPE_FLUID;
		ps.push_back(pi);

	}

	for (int i = 0; i < ps.size(); i++) {
		ps[i].mass = totMass / ps.size();
		ps[i].smth = ps[i].mass / ps[i].dens;

	}
#endif



#ifdef RESTART
	int nparts = 200;
		double xmin = -PARAM::xi;
		double xmax = PARAM::xi;
		char str[80];
		float f;
		FILE * pFile;
		pFile = fopen("result/20000.dat", "r");

		for (int i = 0; i < nparts; i++) {
			Particle pi;

			fscanf(pFile, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &pi.id, &pi.mass, &pi.pos.x, &pi.dens, &pi.eng, &pi.pres, &pi.acc.x,
					&pi.eng_dot, &pi.vel.x, &pi.smth, &pi.mu, &pi.temp, &pi.NUMDENS, &pi.abundances[0], &pi.abundances[5]);
			pi.vel.x = 0.0;
			pi.TYPE = TYPE_FLUID;
			ps.push_back(pi);
	//		cout << pi.eng << endl;
		}
#endif


	char filename[256];
	sprintf(filename, "result/init.dat", 0);
	FILE* fp;
	fp = fopen(filename, "w");

	for (int i = 0; i < ps.size(); i++) {
		//				std::cout << ps[i].pos.x << " "<< ps[i].dens<<std::endl;

		fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", ps[i].id, ps[i].mass, ps[i].pos.x, ps[i].dens, ps[i].eng, ps[i].pres,
				ps[i].acc.x, ps[i].eng_dot, ps[i].vel.x, ps[i].smth, ps[i].mu, ps[i].temp, ps[i].NUMDENS, ps[i].abundances[0], ps[i].abundances[5]);
	}
	Domain domain;
	domain.min.x = xmin;
	domain.max.x = xmax;
	domain.domain_len.x = xmax - xmin;

	double glb_dt = 1e30;

//	copyGhosts(domain, ps, nparts);
	for (int i = 0; i < 10; i++) {
		calc_density(ps, nparts);
	}
	calc_presure(ps, nparts, glb_dt);
	calc_force_G(ps, nparts, glb_dt);
//	copyGhosts(domain, ps, nparts);

	double passtime = 0.0;
	for (int step = 1; step < 20001; step++) {

		passtime += glb_dt;
		cout << "in main time passed:  " << passtime * PARAM::ST / PARAM::yr << " step: " << step << endl;

//		nparts = ps.size();
//		InitialKick(ps, glb_dt, nparts);
//		FullDrift(ps, glb_dt, nparts);
//		Predict(ps, glb_dt, nparts);
//		nparts = ps.size();
		for (int i = 0; i < 10; i++) {

			calc_density(ps, nparts);
		}
		calc_presure(ps, nparts, glb_dt);
		calc_force_G(ps, nparts, glb_dt);
		FinalKick(ps, glb_dt, nparts);
		for (int i = 0; i < 10; i++) {

			calc_density(ps, nparts);
		}
//		copyGhosts(domain, ps, nparts);
		nparts = ps.size();

		if (step % 1000 == 0) {
			char filename[256];
			sprintf(filename, "result/%04d.dat", step);
			FILE* fp;
			fp = fopen(filename, "w");
			for (int i = 0; i < ps.size(); i++) {
//				std::cout << ps[i].pos.x << " "<< ps[i].dens<<std::endl;

				fprintf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", ps[i].id, ps[i].mass, ps[i].pos.x, ps[i].dens, ps[i].eng, ps[i].pres,
						ps[i].acc.x, ps[i].eng_dot, ps[i].vel.x, ps[i].smth, ps[i].mu, ps[i].temp, ps[i].NUMDENS, ps[i].abundances[0], ps[i].abundances[5]);
			}
		}
	}

//	cout << res[0][3].pos << endl;
	return 0;
}
