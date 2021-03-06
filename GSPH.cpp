#include "arraylist.cpp"
#include "relocarray.cpp"
#include "radix_sort.cpp"
#include "ps_defs.cpp"
#include "kernel.h"
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
//	nparts = ps.size();
	double hmin = 1e30;
	for (int i = 0; i < nparts; i++) {
		double tmp_dens = 0.0;
		ps[i].dens = 0.0;
		ps[i].OMEGA = 0.0;
//		nparts = ps.size();
		RellocatableArray<Particle> nbrs;
		find_nnbs(ps[i], ps, nbrs, nparts);
#ifdef I2000_VARIABLE_SMTH
		for (int j = 0; j < nbrs.size(); j++) {
			F64vec rij = ps[i].pos - nbrs[j].pos;
//			ps[i].dens += nbrs[j].mass * ker.W(rij, ps[i].smth);
			tmp_dens += nbrs[j].mass * ker.W(rij, CSMOOTH * ps[i].smth);
			ps[i].OMEGA += nbrs[j].mass * ker.DW_h(rij, ps[i].smth);
		}

		ps[i].smth = ps[i].mass / tmp_dens;
		hmin = fmin(ps[i].smth,hmin);
		for (int j = 0; j < nbrs.size(); j++) {
			F64vec rij = ps[i].pos - nbrs[j].pos;
			ps[i].dens += nbrs[j].mass * ker.W(rij, ps[i].smth);
		}
		ps[i].OMEGA *= ps[i].smth / (1. * ps[i].dens);
		ps[i].OMEGA = 1.0;

#elif NAIVE_VARIABLE_SMTH
		for (int j = 0; j < nbrs.size(); j++) {
			F64vec rij = ps[i].pos - nbrs[j].pos;
			ps[i].dens += nbrs[j].mass * ker.W(rij, ps[i].smth);
		}
		ps[i].smth = ps[i].mass / ps[i].dens;
#else
		for (int j = 0; j < nbrs.size(); j++) {
			F64vec rij = ps[i].pos - nbrs[j].pos;
			ps[i].dens += nbrs[j].mass * ker.W(rij, ps[i].smth);
		}
		ps[i].smth = 0.01;
#endif
	}
//	for (int i = 0; i < nparts; i++) {
//		ps[i].smth = (ps[i].smth > 2. * hmin) ? 2. * hmin : ps[i].smth;
//	}
}
void calc_presure(Particles &ps, long nparts, double &glb_dt) {
	kernel_t ker;
//	nparts = ps.size();
	for (int i = 0; i < nparts; i++) {
		ps[i].pres = ps[i].dens * ps[i].eng * (GAMMA - 1.0);
		ps[i].snds = sqrt(GAMMA * ps[i].pres / ps[i].dens);
		ps[i].dt = 0.5 * ps[i].smth / ps[i].snds;
		glb_dt = fmin(glb_dt, ps[i].dt);
	}
}
//TODO higher dimension notice
void copyGhosts(const Domain &d, Particles &ps, long nparts) {
	double srch = 6.;
//		if (ps[i].TYPE == TYPE_GHOST) {
	ps.remove_ghost(nparts);

//		}
//	nparts = ps.size();
	for (int i = 0; i < nparts; i++) {
		if (ps[i].pos.x < 0.0 || ps[i].pos.x >= 0.8) {
			ps.remove(ps[i].id);
		}
		double partsx = ps[i].pos.x;
		if (fabs(partsx - d.min.x) < srch * ps[i].smth) {
//			ps[i].pos.x = 2.0 * d.min.x - partsx;
//			ps[i].vel *= -1.;
			Particle ghost;
			ghost.pos.x = 2.0 * d.min.x - partsx;
			ghost.acc *= -1;
			ghost.mass = ps[i].mass;
			ghost.dens = ps[i].dens;
			ghost.pres = ps[i].pres;
			ghost.eng = ps[i].eng;
			ghost.smth = ps[i].smth;
			ghost.vel *= -1.0;
			ghost.id = ps[i].id + nparts;
			ghost.TYPE = TYPE_GHOST;
			ps.push_back(ghost);
		} else if (fabs(partsx - d.max.x) < srch * ps[i].smth) {
//			ps[i].pos.x = 2.0 * d.max.x - partsx;
//			ps[i].vel *= -1.;
			Particle ghost;
			ghost.pos.x = 2.0 * d.max.x - partsx;
			ghost.acc *= -1;
			ghost.mass = ps[i].mass;
			ghost.dens = ps[i].dens;
			ghost.pres = ps[i].pres;
			ghost.smth = ps[i].smth;
			ghost.eng = ps[i].eng;
			ghost.vel *= -1.0;
			ghost.id = ps[i].id + nparts;
			ghost.TYPE = TYPE_GHOST;
			ps.push_back(ghost);
		}
	}
//	for (int i = 0; i < nparts; i++) {
//		double partsx = ps[i].pos.x;
//		if(fabs(ps[i].pos.x -0.8)<srch * ps[i].smth){
//			ps[i].pos.x -= 0.80;
//		}
//		if(fabs(ps[i].pos.x -0.)<srch * ps[i].smth){
//					ps[i].pos.x += 0.80;
//				}
//	}
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
		if(ps[i].pos.x<0.4&&ps[i].pos.x>0.39){
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
#ifdef I2000_SECOND_ORDER_RIEMMAN_SOLVER
	F64 drdsi, drdsj, dpdsi, dpdsj, dvdsi, dvdsj;

	F64 si = .5 * delta;
	F64 sj = .5 * delta;
	drdsi = ep_i.graddens * eij;
	drdsj = ep_j.graddens * eij;
	dpdsi = ep_i.gradpres * eij;
	dpdsj = ep_j.gradpres * eij;
	dvdsi = ep_i.gradvel_x * eij*eij.x;
	dvdsj = ep_j.gradvel_x * eij *eij.x;
//	std::cout<<drdsi<<std::endl;
	if (dvdsi * dvdsj <0.0||drdsi *drdsj<0.0||dpdsi *dpdsj<0.0) {
		dvdsi = dvdsj = 0.0;
		drdsi = drdsj = 0.0;
		dpdsi = dpdsj = 0.0;
	}
	if (3.0 * fabs(vj - vi) > ((ep_i.snds < ep_j.snds) ? ep_i.snds : ep_j.snds)) {
		//		//monotonisity constraint

		drdsi = drdsj = 0.0;
		dpdsi = dpdsj = 0.0;
		dvdsi = dvdsj = 0.0;
	}
	double ddensi = drdsi * (ss + ep_i.snds * 0.5 * dt - si);
	double dpresi = dpdsi * (ss + ep_i.snds * 0.5 * dt - si);
	double dveli = dvdsi * (ss + ep_i.snds * 0.5 * dt - si);
	double ddensj = drdsj * (ss - ep_j.snds * 0.5 * dt + sj);
	double dpresj = dpdsj * (ss - ep_j.snds * 0.5 * dt + sj);
	double dvelj = dvdsj * (ss - ep_j.snds * 0.5 * dt +sj);
//	double ddensi = fabs(drdsi) * (ss + ep_i.snds * 0.5 * dt - si);
//	double dpresi = fabs(dpdsi) * (ss + ep_i.snds * 0.5 * dt - si);
//	double dveli = fabs(dvdsi) * (ss + ep_i.snds * 0.5 * dt - si);
//	double ddensj = fabs(drdsj) * (ss - ep_j.snds * 0.5 * dt + sj);
//	double dpresj = fabs(dpdsj) * (ss - ep_j.snds * 0.5 * dt + sj);
//	double dvelj = fabs(dvdsj) * (ss - ep_j.snds * 0.5 * dt +sj);
// ddensi =  0.0;
// dpresi =  0.0;
// ddensj =  0.0;
// dpresj =  0.0;
// dveli =   0.0;
// dvelj =   0.0;

	rho1 = ep_i.dens + ddensi;
	p1 = ep_i.pres + dpresi;
	v1 = vi + dveli;
	rho2 = ep_j.dens + ddensj;
	p2 = ep_j.pres + dpresj;
	v2 = vj + dvelj;
	if (rho1 < 0.0||p1 < 0.0||rho2 < 0.0||p2 < 0.0) {
		rho1 = ep_i.dens;
		p1 = ep_i.pres;
		v1 = vi;
		rho2 = ep_j.dens;
		p2 = ep_j.pres;
		v2 = vj;
	}

#elif VANLEER__SECOND_ORDER_RIEMMAN_SOLVER
	F64vec dr = ep_i.pos - ep_j.pos;
	F64 dsij = sqrt(dr*dr);
	F64 drdsi, drdsj, dpdsi, dpdsj, dvdsi, dvdsj;

	F64 si = .5 * delta;
	F64 sj = .5 * delta;
	drdsi = ep_i.graddens * eij;
	drdsj = ep_j.graddens * eij;
	dpdsi = ep_i.gradpres * eij;
	dpdsj = ep_j.gradpres * eij;
	dvdsi = ep_i.gradvel_x * eij * eij.x;
	dvdsj = ep_j.gradvel_x * eij * eij.x;
	F64 dri = dsij*drdsi;
	F64 dpi = dsij*dpdsi;
	F64 dvi = dsij*dvdsi;

	F64 drj = dsij*drdsj;
	F64 dpj = dsij*dpdsj;
	F64 dvj = dsij*dvdsj;

	F64 dri_dash =2*dri-(ep_i.dens-ep_j.dens);
	F64 dpi_dash = 2*dpi-(ep_i.pres-ep_j.pres);
	F64 dvi_dash =2*dvi-(vi-vj);

	F64 drj_dash = 2*drj-(ep_j.dens-ep_i.dens);
	F64 dpj_dash = 2*dpj-(ep_j.pres-ep_i.pres);
	F64 dvj_dash = 2*dvj-(vj-vi);
	bool same_densi_sign = (std::signbit(ep_i.dens-ep_j.dens)+
			std::signbit(dri)+
			std::signbit(dri_dash)==3)| (std::signbit(ep_i.dens-ep_j.dens)+
			std::signbit(dri)+
			std::signbit(dri_dash)==0);
	bool same_presi_sign = (std::signbit(ep_i.pres-ep_j.pres)+
			std::signbit(dpi)+
			std::signbit(dpi_dash)==3)| (std::signbit(ep_i.pres-ep_j.dens)+
			std::signbit(dpi)+
			std::signbit(dpi_dash)==0);
	bool same_veli_sign = (std::signbit(vj-vi)+
			std::signbit(dvi)+
			std::signbit(dvi_dash)==3)| (std::signbit(vj-vi)+
			std::signbit(dvi)+
			std::signbit(dvi_dash)==0);
	bool same_densj_sign = (std::signbit(ep_j.dens-ep_i.dens)+
			std::signbit(drj)+
			std::signbit(drj_dash)==3)| (std::signbit(ep_j.dens-ep_i.dens)+
			std::signbit(drj)+
			std::signbit(drj_dash)==0);
	bool same_presj_sign = (std::signbit(ep_j.pres-ep_i.pres)+
			std::signbit(dpj)+
			std::signbit(dpj_dash)==3)| (std::signbit(ep_j.pres-ep_i.dens)+
			std::signbit(dpj)+
			std::signbit(dpj_dash)==0);
	bool same_velj_sign = (std::signbit(vi-vj)+
			std::signbit(dvj)+
			std::signbit(dvj_dash)==3)| (std::signbit(vi-vj)+
			std::signbit(dvj)+
			std::signbit(dvj_dash)==0);
	F64 drmonoi =(same_densi_sign)? (fmin(fmin(2.0*fabs(ep_i.dens-ep_j.dens),fabs(dri)) ,2.0*fabs(dri_dash))*(dri>0.0?1.0:-1.0)):0.0;
	F64 dpmonoi =(same_presi_sign)?(fmin(fmin(2.0*fabs(ep_i.pres-ep_j.pres),fabs(dpi)) ,2.0*fabs(dpi_dash))*(dpi>0.0?1.0:-1.0)):0.0;
	F64 dvmonoi =(same_veli_sign)?(fmin(fmin(2.0*fabs(vi-vj),fabs(dvi)) ,2.0*fabs(dvi_dash))*(dvi>0.0?1.0:-1.0)):0.0;

	F64 drmonoj =(same_densj_sign)? (fmin(fmin(2.0*fabs(ep_j.dens-ep_i.dens),fabs(drj)) ,2.0*fabs(drj_dash))*(drj>0.0?1.0:-1.0)):0.0;
	F64 dpmonoj =(same_presj_sign)?(fmin(fmin(2.0*fabs(ep_j.pres-ep_i.pres),fabs(dpj)) ,2.0*fabs(dpj_dash))*(dpj>0.0?1.0:-1.0)):0.0;
	F64 dvmonoj = (same_velj_sign)?(fmin(fmin(2.0*fabs(vj-vi),fabs(dvj)) ,2.0*fabs(dvj_dash))*(dvj>0.0?1.0:-1.0)):0.0;
	if ((3. * fabs(vj - vi) > ((ep_i.snds < ep_j.snds) ? ep_i.snds : ep_j.snds))) {
		//						std::cout<< ps.graddens<<" "<<ps.grav<<" "<<ps.snds<<" "<<ps.pres<<std::endl;
		drmonoi = 0.0;
		drmonoj= 0.0;
		dpmonoi = 0.0;
		dpmonoj= 0.0;
		dvmonoi = 0.0;
		dvmonoj = 0.0;
	}

	//set initial condition of riemann problem

	double ddensi = .5*drmonoi* ( 1.0- ep_i.snds *dt/dsij );
	double dpresi = .5*dpmonoi* ( 1.0- ep_i.snds *dt/dsij );
	double ddensj = .5*drmonoj*( 1.0- ep_j.snds *dt/dsij );
	double dpresj = .5*dpmonoj*( 1.0- ep_j.snds *dt/dsij );
	double dveli = .5*dvmonoi * ( 1.0- ep_i.snds* dt/dsij);
	double dvelj = .5*dvmonoj*( 1.0- ep_j.snds* dt/dsij);

	rho1 = ep_i.dens + ddensi;
	p1 = ep_i.pres + dpresi;
	v1 = vi + dveli;
	rho2 = ep_j.dens + ddensj;
	p2 = ep_j.pres + dpresj;
	v2 = vj + dvelj;
#else
	rho1 = ep_i.dens;
	p1 = ep_i.pres;
	v1 = vi;
	rho2 = ep_j.dens;
	p2 = ep_j.pres;
	v2 = vj;
#endif
	F64 ppre, p, v;
	F64 W1, W2;
	const F64 alpha = (2.0 * GAMMA) / (GAMMA - 1.0);
	p = .5 * (p1 + p2);
//	p=ep_j.snds*p1+ep_i.snds*p2-ep_i.snds*ep_j.snds*(v1-v2)/(ep_i.snds+ep_j.snds);
	double critif = 1.0 - 1.0e-6;
	p1 += 1e-8;
	p2 += 1e-8;
	if (p1 <= 0.0 || p2 <= 0.0) {
		pstar = .5 * (p1 + p2);
		vstar = .5 * (v1 + v2);
	} else {
		for (U32 loop = 0; loop < 400; loop++) {
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
			if (fabs(p - ppre) < 1e-18)
				break;
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

	calc_gradients(ps, nparts);
	for (S32 i = 0; i < nparts; ++i) {
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
//				const F64vec gradW = 0.5
//									* (kernel.gradW(dr, ps[i].smth)
//											+ kernel.gradW(dr, nbrs[j].smth));
				gradW_hsi = kernel.gradW(dr, sqrt(2.0) * ps[i].smth);
				gradW_hsj = kernel.gradW(dr, sqrt(2.0) * nbrs[j].smth);

				calc_Vij2_and_ss(ps[i], nbrs[j], VIJ_I, VIJ_J, sstar, delta, eij);
				calc_riemann_solver(ps[i], nbrs[j], sstar, delta, eij, glb_dt, PSTAR, VSTAR);
//				interP ips;
//				ips.PSTAR = PSTAR;
//				ips.VIJ_I = VIJ_I;
//				ips.VIJ_J = VIJ_J;
//				ips.VSTAR = VSTAR;
//				ips.sstar = sstar;
//				intervals.push_back(ips);
//				cout << VSTAR << " " << i << " " << j << endl;
//				vij = .5 * (dv);
//				PSTAR = .5 * (ps[i].pres + ps[j].pres);
#ifdef VARIABLE_SMTH
				interpo = (PSTAR*eij) * (gradW_hsi * VIJ_I + gradW_hsj * VIJ_J );
#else
//				VIJ_I =VIJ_J=pow(PSTAR/(ps[i].eng*2/3.),-2.);
				interpo = (PSTAR) * (gradW_hsi * VIJ_I / ps[i].OMEGA + gradW_hsj * VIJ_J / ps[i].OMEGA);
#endif
				ps[i].acc -= nbrs[j].mass * interpo;

			}
		}
		ps[i].vel_half = ps[i].vel + 0.5 * ps[i].acc * glb_dt;
		for (S32 j = 0; j < nbrs.size(); ++j) {
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
//				const F64 si = ps[i].pos * eij;
//				const F64 sj = nbrs[j].pos * eij;
//				const F64 dsij = si - sj;
				F64vec vij = 0.0;
				//				const F64vec gradW = 0.5
				//									* (kernel.gradW(dr, ps[i].smth)
				//											+ kernel.gradW(dr, nbrs[j].smth));
				gradW_hsi = kernel.gradW(dr, sqrt(2.0) * ps[i].smth);
				gradW_hsj = kernel.gradW(dr, sqrt(2.0) * nbrs[j].smth);
				calc_Vij2_and_ss(ps[i], nbrs[j], VIJ_I, VIJ_J, sstar, delta, eij);
				calc_riemann_solver(ps[i], nbrs[j], sstar, delta, eij, glb_dt, PSTAR, VSTAR);
//				PSTAR = intervals[j].PSTAR;
//				VIJ_I = intervals[j].VIJ_I;
//				VIJ_J = intervals[j].VIJ_J;
//				VSTAR = intervals[j].VSTAR;
//				sstar = intervals[j].sstar;
//				calc_Vij2_and_ss(ps[i], nbrs[j], VIJ_I, VIJ_J, sstar, delta,
//						eij);
//				calc_riemann_solver(ps[i], nbrs[j], sstar, delta, eij, glb_dt,
//						PSTAR, VSTAR);
//				double fac1, fac2;
//					if (fabs(dr.x) > 0.0) {
////						fac1 = 0;
////						fac2 = 0;
////						break;
//					}
//					cout<<(fabs(dr.x)>0)<<endl;
//					double vi = VSTAR * fac1;
//					double vj = VSTAR * fac2;
				vij = VSTAR * eij;
				ps[i].eng_dot -= nbrs[j].mass * PSTAR * (vij - ps[i].vel_half) * (gradW_hsi * VIJ_I / ps[i].OMEGA + gradW_hsj * VIJ_J / ps[i].OMEGA);
//					ps[i].eng_dot -= nbrs[j].mass * PSTAR
//												* (vi*gradW_hsi.x * VIJ_I- vhalf * gradW_hsi * VIJ_I
//														+ vj*gradW_hsj.x * VIJ_J- vhalf * gradW_hsj * VIJ_J);
			}
		}

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
	for (S32 i = 0; i < nparts; ++i) {
//		ps[i].vel=0.0;
//		ps[i].eng=0.0;
//		ps[i].vel = ps[i].vel_half + 0.5 * dt * ps[i].acc;
//		ps[i].eng = ps[i].eng_half + 0.5 * dt * ps[i].eng_dot;
//		ps[i].vel += dt * ps[i].acc;
//		ps[i].pos += ps[i].vel * dt;
//		ps[i].eng += dt * ps[i].eng_dot;
		ps[i].pos += ps[i].vel_half * dt;
			ps[i].vel += dt * ps[i].acc;
			ps[i].eng += dt * ps[i].eng_dot;
#ifdef ENERGY_GUARD
		if (ps[i].eng <= 0.0) {
			ps[i].eng -= dt * ps[i].eng_dot;
			ps[i].eng = .5*ps[i].eng;
		}
#endif
	}

}

int main() {
	Particles ps;

	double xmin = 0;
	double xmax = 0.8;
	double domain_len = xmax - xmin;
	int nparts = 0;
	int STEP_LIMIT = 0;
	int id = 0;
#ifdef SodST
	STEP_LIMIT = 80;
	double xmid = .5 * domain_len;
	int npartsl = 80;
	int npartsr = 40;
	nparts = npartsl + npartsr;
	double dxl = xmid / npartsl;
	double dxr = xmid / npartsr;
	for (F64 x = xmin; x < xmin + xmid; x += dxl) {
		id++;
		Particle pi;
		pi.pos = x+dxl*.5;
		pi.id = id;
		pi.vel = 0;
		pi.dens = 1.0;
		pi.mass = pi.dens * dxl;
		pi.pres = 1.0;
		pi.smth = pi.mass / pi.dens;
		pi.eng = pi.pres / (pi.dens * (GAMMA - 1));
		pi.TYPE = TYPE_FLUID;
		ps.push_back(pi);

	}
	for (F64 x = xmin + xmid; x < xmax; x += dxr) {
		id++;
		Particle pi;
		pi.pos = x+dxr*.5;
		pi.id = id;
		pi.vel = 0;
		pi.dens = 0.5;
		pi.mass = pi.dens * dxr;
		pi.pres = 0.2;
		pi.smth = pi.mass / pi.dens;
		pi.eng = pi.pres / (pi.dens * (GAMMA - 1));
		pi.TYPE = TYPE_FLUID;
		ps.push_back(pi);

	}
#elif BLASTWAVE
#ifdef I2000_SECOND_ORDER_RIEMMAN_SOLVER
	STEP_LIMIT = 20;
#else
	STEP_LIMIT = 250;
#endif
	double xmid = .5 * domain_len;
	int npartsl = 200;
	int npartsr = 200;
	nparts = npartsl + npartsr;
	double dxl = xmid / npartsl;
	double dxr = xmid / npartsr;
	for (F64 x = xmin; x < xmin + xmid; x += dxl) {
		id++;
		Particle pi;
		pi.pos = x+dxl*.5;
		pi.id = id;
		pi.vel = 0;
		pi.dens = 1.0;
		pi.mass = pi.dens * dxl;
		pi.pres = 20.0;
		pi.smth = pi.mass / pi.dens;
		pi.eng = pi.pres / (pi.dens * (GAMMA - 1));
		pi.TYPE = TYPE_FLUID;
		ps.push_back(pi);

	}
	for (F64 x = xmin + xmid; x < xmax; x += dxr) {
		id++;
		Particle pi;
		pi.pos = x+dxr*.5;
		pi.id = id;
		pi.vel = 0;
		pi.dens =1.0;
		pi.mass = pi.dens * dxr;
		pi.pres =1.0;
		pi.smth = pi.mass / pi.dens;
		pi.eng = pi.pres / (pi.dens * (GAMMA - 1));
		pi.TYPE = TYPE_FLUID;
		ps.push_back(pi);

	}

#elif SJOGREEN
#ifdef I2000_SECOND_ORDER_RIEMMAN_SOLVER
	STEP_LIMIT = 20;
#else
	STEP_LIMIT = 250;
#endif
	double xmid = .5 * domain_len;
	int npartsl = 200;
	int npartsr = 200;
	nparts = npartsl + npartsr;
	double dxl = xmid / npartsl;
	double dxr = xmid / npartsr;
	for (F64 x = xmin; x < xmin + xmid; x += dxl) {
		id++;
		Particle pi;
		pi.pos = x;
		pi.id = id;
		pi.vel = -2.0;
		pi.dens = 1.0;
		pi.mass = pi.dens * dxl;
		pi.pres = 0.4;
		pi.smth = pi.mass / pi.dens;
		pi.eng = pi.pres / (pi.dens * (GAMMA - 1));
		pi.TYPE = TYPE_FLUID;
		ps.push_back(pi);

	}
	for (F64 x = xmin + xmid; x < xmax; x += dxr) {
		id++;
		Particle pi;
		pi.pos = x;
		pi.id = id;
		pi.vel = 2.0;
		pi.dens =1.0;
		pi.mass = pi.dens * dxr;
		pi.pres =0.4;
		pi.smth = pi.mass / pi.dens;
		pi.eng = pi.pres / (pi.dens * (GAMMA - 1));
		pi.TYPE = TYPE_FLUID;
		ps.push_back(pi);

	}
#endif
#ifdef TRANSPORT
	STEP_LIMIT =90;
	xmin = -M_PI+0.2;
	xmax =M_PI-0.2;
	nparts = 100;
	double dx =( xmax-xmin)/nparts;
	for (F64 x = xmin; x < xmax; x += dx) {
		id++;
		Particle pi;
		pi.pos = x;
		pi.id = id;
		pi.vel = 0.0;
		pi.dens = 1.0 + cos(x);
		pi.mass = pi.dens * dx;
		pi.pres = 1.0;
		pi.smth = pi.mass / pi.dens;
		pi.eng = pi.pres / (pi.dens * (GAMMA - 1));
		pi.TYPE = TYPE_FLUID;
		ps.push_back(pi);

	}
#endif
	FILE* fp;
	fp = fopen("test1.txt", "w");

	Domain domain;
	domain.min.x = xmin;
	domain.max.x = xmax;
	domain.domain_len.x = xmax - xmin;
	double glb_dt = 1e30;
//	copyGhosts(domain, ps, nparts);
	calc_density(ps, nparts);
	calc_presure(ps, nparts, glb_dt);
	calc_force_G(ps, nparts, glb_dt);
	copyGhosts(domain, ps, nparts);

	double passtime = 0.0;
	for (int step = 0; step < 250; step++) {
		passtime += glb_dt;
		cout <<"in main time passed: " <<passtime << endl;
//		nparts = ps.size();
//		InitialKick(ps, glb_dt, nparts);
//		FullDrift(ps, glb_dt, nparts);
//		Predict(ps, glb_dt, nparts);
//		nparts = ps.size();
		calc_density(ps, nparts);
		calc_presure(ps, nparts, glb_dt);
		calc_force_G(ps, nparts, glb_dt);
//		glb_dt =0.001;S
//		cout << glb_dt << endl;
		FinalKick(ps, glb_dt, nparts);
		copyGhosts(domain, ps, nparts);
	}
//	nparts = ps.size();
	for (int i = 0; i < nparts; i++) {
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", ps[i].pos.x, ps[i].dens, ps[i].eng, ps[i].pres, ps[i].acc.x, ps[i].eng_dot,ps[i].vel.x);
	}
//	cout << res[0][3].pos << endl;
	return 0;
}
