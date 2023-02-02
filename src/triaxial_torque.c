/**
 * @file    triaxial_torque.c
 * @brief   Torque on triaxial bodies
 * @author  Henry Yuan
 * 
 * 
 * *** COMMENT SECTION BELOW THIS HAS NOT BEEN CHANGED FROM COPIED GRAVITATIONAL_HARMONICS.C
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Gravity Fields$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    `Tamayo, Rein, Shi and Hernandez, 2019 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T/abstract>`_. 
 * Based on                None
 * C Example               :ref:`c_example_J2`
 * Python Example          `J2.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/J2.ipynb>`_.
 * ======================= ===============================================
 * 
 * * *** COMMENT SECTION ABOVE THIS HAS NOT BEEN CHANGED FROM COPIED GRAVITATIONAL_HARMONICS.C
 * 
 * Adds the effects of a particle having 3 differing moments of inertia (triaxial particle) on the spin vector 
 * of the particle. Assumes that changes to the spin angular momentum are negligible compared to the orbital 
 * angular momentum; i.e., the particle's orbit is not affected by the changes to the spin vector.
 *
 * 
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 * 
 * x,y,z: refers to the x,y,z coordinate system underlying rebound
 * 
 * I define a new coordinate system in addition to the one used by rebound (described above):
 * -> i,j,k: principal axes corresponding to the particle's principal moments of inertia. axis i has the
 *    lowest moment of inertia and axis k has the highest moment of inertia
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * tt_Ii (double)                  Yes         Moment for axis i (<= Ij)
 * tt_Ij (double)                  Yes         Moment for axis j (<= Ik, >= Ii)
 * tt_Ik (double)                  Yes         Moment for axis k (>= Ij)
 * tt_omega (double)               Yes         spin rate
 * tt_ix (double)                  Yes         x-component of unit vector for lowest moment axis
 * tt_iy (double)                  Yes         y-component of unit vector for lowest moment axis
 * tt_iz (double)                  Yes         z-component of unit vector for lowest moment axis
 * tt_jx (double)                  Yes         x-component of unit vector for middle moment axis
 * tt_jy (double)                  Yes         y-component of unit vector for middle moment axis
 * tt_jz (double)                  Yes         z-component of unit vector for middle moment axis
 * tt_kx (double)                  Yes         x-component of unit vector for highest moment axis
 * tt_ky (double)                  Yes         y-component of unit vector for highest moment axis
 * tt_kz (double)                  Yes         z-component of unit vector for highest moment axis
 * tt_si (double)                  Yes         i component of spin unit vector
 * tt_sj (double)                  Yes         j component of spin unit vector
 * tt_sk (double)                  Yes         k component of spin unit vector
 * tt_tidal_dt  (double)           Yes         constant tidal timelag
 * tt_k2 (double)                  Yes         2nd degree Love number
 * tt_R  (double)                  Yes         mean radius
 * ============================ =========== ==================================================================
 * 
 * Parameter requirements:
 * ============================
 * i, j, k, and s must be unit vectors (e.g. such that ix^2 + iy^2 + iz^2 = 1)
 * i cross j = k
 */

/* DEBUGGING CHANGES:
- commented out torque calc
- domega_dts_ijk all 0 
- print statements
- only using first derivative calc
*/

/* CHANGES:
- calling dijk_dt_acc instead of old one
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

// global variables
double PI = 3.14159265358979323846;

// computes and returns the dot product between vectors u and v
static double rebx_dot_prod(double u[3], double v[3]) {
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

// computes cross product between vectors u and v, and puts result in w
static void rebx_cross_prod(double u[3], double v[3], double w[3]) {
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
}

// convert vector ijk to xyz (same vector in xyz basis), given ijk_xyz
static void rebx_ijk_to_xyz(double ijk[3], double xyz[3], double ijk_xyz[3][3]){
    for (int i=0; i<3; i++) {
        xyz[i] = ijk[0]*ijk_xyz[0][i] + ijk[1]*ijk_xyz[1][i] + ijk[2]*ijk_xyz[2][i];
    }
}

// convert vector xyz to ijk (same vector in ijk basis), given ijk_xyz
static void rebx_xyz_to_ijk(double xyz[3], double ijk[3], double ijk_xyz[3][3]){
    for (int i=0; i<3; i++) {
        ijk[i] = rebx_dot_prod(xyz,ijk_xyz[i]);
    }
}

// linearly interpolates position of particle p at time dt
static void rebx_interpolate_xyz(struct reb_particle* p, double xyz[3], double dt){
    xyz[0] = p->x + (dt*p->vx);
    xyz[1] = p->y + (dt*p->vy);
    xyz[2] = p->z + (dt*p->vz);
}

// interpolate position assuming orbit around stationary primary
// orbit assumed to be near circular??
        // *** when calling with p = primary, dtheta not used, make sure to incorporate additional time offset into dt ***
static void rebx_interpolate_xyz_acc(double sim_G, struct reb_particle* p, struct reb_particle* primary, double xyz[3], double dt){
    
    // dtheta not used if p == primary
    if (p == primary) {
        xyz[0] = p->x + (dt*p->vx);
        xyz[1] = p->y + (dt*p->vy);
        xyz[2] = p->z + (dt*p->vz);
        return;
    }

    struct reb_orbit o = reb_tools_particle_to_orbit(sim_G, *p, *primary);

    double r_xyz[3] = {p->x - primary->x,p->y - primary->y,p->z - primary->z};
    double v_xyz[3] = {p->vx - primary->vx, p->vy - primary->vy, p->vz - primary->vz};
    double r = sqrt(r_xyz[0]*r_xyz[0] + r_xyz[1]*r_xyz[1] + r_xyz[2]*r_xyz[2]);
    double v = sqrt(v_xyz[0]*v_xyz[0] + v_xyz[1]*v_xyz[1] + v_xyz[2]*v_xyz[2]);

    double dtheta = (dt / o.P) * 2*PI;
    double cos_dtheta = cos(dtheta);
    double sin_dtheta = sin(dtheta);

    xyz[0] = p->x + (dt*primary->vx) + r*sin_dtheta*v_xyz[0]/v - (1-cos_dtheta)*r_xyz[0];
    xyz[1] = p->y + (dt*primary->vy) + r*sin_dtheta*v_xyz[1]/v - (1-cos_dtheta)*r_xyz[1];
    xyz[2] = p->z + (dt*primary->vz) + r*sin_dtheta*v_xyz[2]/v - (1-cos_dtheta)*r_xyz[2];
    // printf("a=%f e=%f\n", o.a, o.e);
}

// computes time-derivative of spin vector omega_i,omega_j,omega_k using torque vector Mi, Mj, Mk according to Euler's equations
static void rebx_domega_dt(double omega_ijk[3], double M_ijk[3], const double I_ijk[3], double domega_dts_ijk[3]){
    domega_dts_ijk[0] = (M_ijk[0] + (I_ijk[1]-I_ijk[2])*omega_ijk[1]*omega_ijk[2]) / I_ijk[0];
    domega_dts_ijk[1] = (M_ijk[1] + (I_ijk[2]-I_ijk[0])*omega_ijk[2]*omega_ijk[0]) / I_ijk[1];
    domega_dts_ijk[2] = (M_ijk[2] + (I_ijk[0]-I_ijk[1])*omega_ijk[0]*omega_ijk[1]) / I_ijk[2];
    // domega_dts_ijk[0] = 0.0; // [DEBUG]
    // domega_dts_ijk[1] = 0.0; // [DEBUG]
    // domega_dts_ijk[2] = 0.0; // [DEBUG]
}

// computes time-derivative of vectors i,j,k (components in old ijk basis)
static void rebx_dijk_dt(double ijk_ijk[3][3], double omega_ijk[3], double dijk_dts_ijk[3][3]){
    rebx_cross_prod(omega_ijk,ijk_ijk[0],dijk_dts_ijk[0]); // di/dt
    rebx_cross_prod(omega_ijk,ijk_ijk[1],dijk_dts_ijk[1]); // dj/dt
    rebx_cross_prod(omega_ijk,ijk_ijk[2],dijk_dts_ijk[2]); // dk/dt
}

// calculates the triaxial torque from all other bodies on the 'index'th particle
static void rebx_calc_triax_torque(struct reb_simulation* const sim, int index, double M_ijk[3], const double I_ijk[3], double ijk_xyz[3][3], double dt, const double sim_dt){
    
    /*************BELOW FOR JUST TORQUE FROM HOST STAR**********/
    // if primary, ignore
    if (index == 0) {
        return;
    }
    /*************ABOVE FOR JUST TORQUE FROM HOST STAR**********/

    struct reb_particle* p = &sim->particles[index];
    struct reb_particle* torquer;
    double p_xyz[3];
    double torquer_xyz[3];
    double r_xyz[3];
    double r;
    double r_dot_i;
    double r_dot_j;
    double r_dot_k;
    double prefac;

    // rebx_interpolate_xyz(p,p_xyz,dt-sim_dt);
    rebx_interpolate_xyz_acc(sim->G, p, &sim->particles[0],p_xyz,dt-sim_dt);

    /*************BELOW FOR TORQUE FROM ALL OTHER BODIES**********/
    // const int _N_real = sim->N - sim->N_var;
	// for(int i=0; i<_N_real; i++){
    //     if (i == index) {
    //         continue;
    //     }
    //    torquer = &sim->particles[i];
    /*************ABOVE FOR TORQUE FROM ALL OTHER BODIES**********/
    torquer = &sim->particles[0];
    // rebx_interpolate_xyz(torquer,torquer_xyz,dt-sim_dt);
    rebx_interpolate_xyz_acc(sim->G, torquer, &sim->particles[0],torquer_xyz,dt-sim_dt);
    r_xyz[0] = torquer_xyz[0] - p_xyz[0];
    r_xyz[1] = torquer_xyz[1] - p_xyz[1];
    r_xyz[2] = torquer_xyz[2] - p_xyz[2];
    r = sqrt(rebx_dot_prod(r_xyz, r_xyz));
    prefac = 3 * sim->G * torquer->m / pow(r,5);

    r_dot_i = rebx_dot_prod(r_xyz,ijk_xyz[0]);
    r_dot_j = rebx_dot_prod(r_xyz,ijk_xyz[1]);
    r_dot_k = rebx_dot_prod(r_xyz,ijk_xyz[2]);

    M_ijk[0] += prefac*(I_ijk[2]-I_ijk[1])*r_dot_j*r_dot_k;
    M_ijk[1] += prefac*(I_ijk[0]-I_ijk[2])*r_dot_k*r_dot_i;
    M_ijk[2] += prefac*(I_ijk[1]-I_ijk[0])*r_dot_i*r_dot_j;

        // printf("triax: %.10e\n", M_ijk[2]);  // DEBUG
    // }
}

// calculates the tidal torque from host star on the 'index'th particle
// assumes host star is at index 0
static void rebx_calc_tidal_torque(struct reb_simulation* const sim, int index, double M_ijk[3], double omega_ijk[3], double ijk_xyz[3][3], const double tidal_dt, const double k2, 
    const double R, double dt, const double sim_dt){

    // if primary, ignore
    if (index == 0) {
        return;
    }

    struct reb_particle* p = &sim->particles[index];
    struct reb_particle* primary = &sim->particles[0];
    struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, *p, *primary);
 
    double orbit_normal_ijk[3];
    double p_xyz[3];
    double pv_xyz[3] = {p->vx,p->vy,p->vz};
    double pv_ijk[3];
    double primary_xyz[3];
    double r_xyz[3];
    double r_ijk[3];
    double r;
    double rho_ijk[3];
    double rho;
    double prefac;
    double rho_cross_r_ijk[3];
    double norm_cross_r_ijk[3];
    double omega_cross_r_ijk[3];

    rebx_interpolate_xyz_acc(sim->G,p,primary,p_xyz,dt-sim_dt);
    rebx_interpolate_xyz_acc(sim->G,primary,primary,primary_xyz,dt-sim_dt);

    // calculate r vector
    r_xyz[0] = primary_xyz[0] - p_xyz[0];
    r_xyz[1] = primary_xyz[1] - p_xyz[1];
    r_xyz[2] = primary_xyz[2] - p_xyz[2];
    r = sqrt(rebx_dot_prod(r_xyz,r_xyz));
    rebx_xyz_to_ijk(r_xyz,r_ijk,ijk_xyz);
    
    // calculate orbit normal vector by v cross r
    rebx_xyz_to_ijk(pv_xyz,pv_ijk,ijk_xyz);
    rebx_cross_prod(pv_ijk,r_ijk,orbit_normal_ijk);
    double orbit_normal_norm = sqrt(rebx_dot_prod(orbit_normal_ijk,orbit_normal_ijk));
    orbit_normal_ijk[0] /= orbit_normal_norm;
    orbit_normal_ijk[1] /= orbit_normal_norm;
    orbit_normal_ijk[2] /= orbit_normal_norm;

    // calculate rho_ijk
    rebx_cross_prod(orbit_normal_ijk,r_ijk,norm_cross_r_ijk);
    rebx_cross_prod(omega_ijk,r_ijk,omega_cross_r_ijk);
    rho_ijk[0] = r_ijk[0] + tidal_dt*(omega_cross_r_ijk[0]-o.n*norm_cross_r_ijk[0]);
    rho_ijk[1] = r_ijk[1] + tidal_dt*(omega_cross_r_ijk[1]-o.n*norm_cross_r_ijk[1]);
    rho_ijk[2] = r_ijk[2] + tidal_dt*(omega_cross_r_ijk[2]-o.n*norm_cross_r_ijk[2]);

    rho = sqrt(rebx_dot_prod(rho_ijk,rho_ijk));

    // print components of r and rho
    // printf("xyz: %.5e, %.5e, %.5e\n", p->x,p->y,p->z);
    // printf("r: %.5e, %.5e, %.5e\n", r_ijk[0], r_ijk[1], r_ijk[2]);
    // printf("rho: %.5e, %.5e, %.5e\n", rho_ijk[0], rho_ijk[1], rho_ijk[2]);

    double r_dot_rho = rebx_dot_prod(r_ijk,rho_ijk);

    // check if tidal phase lag is larger than pi / 4
    double epsilon = acos(r_dot_rho/r/rho);
    if (fabs(epsilon) > (PI/50)) {
        if (fabs(epsilon) > (PI/4)) {
            fprintf(stderr, "REBOUNDx WARNING: triaxial_torque: Tidal phase lag is greater than 45 degrees; tidal effects are no longer accurate.\n");
        }
        // else {
        //     fprintf(stderr, "REBOUNDx WARNING: triaxial_torque: Tidal phase lag is greater than 1/100 of a rotation; tidal effects may not be accurate.\n");
        // }
    }
   
    prefac = 3*k2*sim->G*primary->m*primary->m*pow(R,5)*r_dot_rho / (pow(rho,2)*pow(r,8));

    rebx_cross_prod(rho_ijk,r_ijk,rho_cross_r_ijk);
    M_ijk[0] += prefac*rho_cross_r_ijk[0];
    M_ijk[1] += prefac*rho_cross_r_ijk[1];
    M_ijk[2] += prefac*rho_cross_r_ijk[2];

    // if (prefac*rho_cross_r[2] != 0.0) {
    //     printf("omega = %.15e\n", omega);
    //     printf("n = %.15e\n", o.n);
    // printf("tide: %.15e, %.15e, %.15e\n", prefac*rho_cross_r[0],prefac*rho_cross_r[1],prefac*rho_cross_r[2]);  // DEBUG
    // }
}

/* updates spin vector, omega, and ijk in lockstep using 4th order Runge Kutta.
If calc_torque_bool = 0, torques (including both triaxial and tidal) NOT calculated, otherwise torques calculated */
static void rebx_update_spin_ijk(struct reb_simulation* const sim, int index, double* const ix, double* const iy, double* const iz, 
    double* const jx, double* const jy, double* const jz, double* const kx, double* const ky, double* const kz, double* const si, double* const sj,
    double* const sk, double* const omega, const double Ii, const double Ij, const double Ik, const double tidal_dt, const double k2, const double R, const double dt){

    // Array for principal moments
    const double I_ijk[3] = {Ii,Ij,Ik};

    // Declare matrices for all R-K calculations
    double rk_M_ijk[4][3] = {}; // ijk components of torque on body, Needs all values initialized to zero because multiple functions add to the values
    double rk_omega_ijk[4][3]; // ijk components of spin vector (total omega vector)
    double rk_ijk_xyz[4][3][3]; // xyz components of each ijk vector
    double rk_ijk_ijk[4][3][3]; // components of each ijk vector in the old ijk basis
    double rk_domega_dts_ijk[4][3]; // matrix for all calculations of domega/dt in old ijk basis
    double rk_dijk_dts_ijk[4][3][3]; // matrix for all calculations of d{ijk}/dt in old ijk basis
    /* second dimension: vector (i_hat,j_hat,k_hat)
    third dimension: component of vector (i,j,k) */
    double rk_dts[4] = {0.0, 0.5*dt, 0.5*dt, dt}; // array of sub-timesteps for each RK calculation

    // initialize omega_0
    rk_omega_ijk[0][0] = *omega**si;
    rk_omega_ijk[0][1] = *omega**sj;
    rk_omega_ijk[0][2] = *omega**sk;

    // initialize xyz components of i,j,k
    rk_ijk_xyz[0][0][0] = *ix;
    rk_ijk_xyz[0][0][1] = *iy;
    rk_ijk_xyz[0][0][2] = *iz;
    rk_ijk_xyz[0][1][0] = *jx;
    rk_ijk_xyz[0][1][1] = *jy;
    rk_ijk_xyz[0][1][2] = *jz;
    rk_ijk_xyz[0][2][0] = *kx;
    rk_ijk_xyz[0][2][1] = *ky;
    rk_ijk_xyz[0][2][2] = *kz;

    // initialize ijk components of i,j,k
    rk_ijk_ijk[0][0][0] = 1.0;
    rk_ijk_ijk[0][0][1] = 0.0;
    rk_ijk_ijk[0][0][2] = 0.0;
    rk_ijk_ijk[0][1][0] = 0.0;
    rk_ijk_ijk[0][1][1] = 1.0;
    rk_ijk_ijk[0][1][2] = 0.0;
    rk_ijk_ijk[0][2][0] = 0.0;
    rk_ijk_ijk[0][2][1] = 0.0;
    rk_ijk_ijk[0][2][2] = 1.0;

    // Runge-Kutta Calculations
    for (int i = 0; i < 4; i++) {
        // Pre-calcs, skip if first iteration
        if (i != 0) {
            for (int j=0; j < 3; j++) {
                rk_omega_ijk[i][j] = rk_domega_dts_ijk[i-1][j]*rk_dts[i] + rk_omega_ijk[0][j];
                for (int k=0; k < 3; k++) {
                    rk_ijk_ijk[i][j][k] = rk_dijk_dts_ijk[i-1][j][k]*rk_dts[i] + rk_ijk_ijk[0][j][k];
                }
                rebx_ijk_to_xyz(rk_ijk_ijk[i][j],rk_ijk_xyz[i][j],rk_ijk_xyz[0]);
            }
        }

        // Calcs
        rebx_calc_triax_torque(sim,index,rk_M_ijk[i],I_ijk,rk_ijk_xyz[i],rk_dts[i],dt); // [DEBUG]
        rebx_calc_tidal_torque(sim,index,rk_M_ijk[i],rk_omega_ijk[i],rk_ijk_xyz[i],tidal_dt,k2,R,rk_dts[i],dt);

        rebx_domega_dt(rk_omega_ijk[i],rk_M_ijk[i],I_ijk,rk_domega_dts_ijk[i]);
        rebx_dijk_dt(rk_ijk_ijk[i],rk_omega_ijk[i],rk_dijk_dts_ijk[i]);
    }
    
    // calculate domega_ijk, d{ijk}
    double domega_ijk[3]; // in old ijk basis
    double dijk_ijk[3][3]; // in old ijk basis
    for (int i = 0; i < 3; i++){
        // domega_ijk[i] = rk_domega_dts_ijk[0][i] * dt; // [DEBUG]
        domega_ijk[i] = (rk_domega_dts_ijk[0][i] + 2*rk_domega_dts_ijk[1][i] + 2*rk_domega_dts_ijk[2][i] + rk_domega_dts_ijk[3][i]) * dt / 6;
        for (int j = 0; j < 3; j++){
            // dijk_ijk[i][j] = rk_dijk_dts_ijk[0][i][j] * dt; // [DEBUG]
            dijk_ijk[i][j] = (rk_dijk_dts_ijk[0][i][j] + 2*rk_dijk_dts_ijk[1][i][j] + 2*rk_dijk_dts_ijk[2][i][j] + rk_dijk_dts_ijk[3][i][j]) * dt / 6;
        }
    }

    // printf("%f\n", domega_ijk[0]); // [DEBUG]
    // printf("%f\n", domega_ijk[1]); // [DEBUG]
    // printf("%f\n", domega_ijk[2]); // [DEBUG]

    // update omega and s vector
    double omega_ijk[3] = {rk_omega_ijk[0][0] + domega_ijk[0],rk_omega_ijk[0][1] + domega_ijk[1],rk_omega_ijk[0][2] + domega_ijk[2]};
    *omega = sqrt(rebx_dot_prod(omega_ijk,omega_ijk));
    *si = omega_ijk[0] / *omega;
    *sj = omega_ijk[1] / *omega;
    *sk = omega_ijk[2] / *omega;

    // calculate new i,j,k vectors, convert to xyz basis
    double ijk_ijk[3][3];
    double ijk_xyz[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ijk_ijk[i][j] = rk_ijk_ijk[0][i][j] + dijk_ijk[i][j];
        }
        rebx_ijk_to_xyz(ijk_ijk[i],ijk_xyz[i],rk_ijk_xyz[0]);
    }

    // re-normalize
    double i_mag = sqrt(rebx_dot_prod(ijk_xyz[0],ijk_xyz[0]));
    double j_mag = sqrt(rebx_dot_prod(ijk_xyz[1],ijk_xyz[1]));
    double k_mag = sqrt(rebx_dot_prod(ijk_xyz[2],ijk_xyz[2]));

    *ix = ijk_xyz[0][0] / i_mag;
    *iy = ijk_xyz[0][1] / i_mag;
    *iz = ijk_xyz[0][2] / i_mag;
    *jx = ijk_xyz[1][0] / j_mag;
    *jy = ijk_xyz[1][1] / j_mag;
    *jz = ijk_xyz[1][2] / j_mag;
    *kx = ijk_xyz[2][0] / k_mag;
    *ky = ijk_xyz[2][1] / k_mag;
    *kz = ijk_xyz[2][2] / k_mag;
}

// runs checks on parameters. Returns 1 if error, 0 otherwise.
static int rebx_validate_params(struct reb_simulation* const sim, const double* const Ii, const double* const Ij,
    const double* const Ik, double* const omega, double* const ix, double* const iy, double* const iz, double* const jx,
    double* const jy, double* const jz, double* const kx, double* const ky, double* const kz, double* const si, double* const sj,
    double* const sk, const double* const tidal_dt, const double* const k2, const double* const R) {

    double tolerance = 1.e-15;

    double i_diff = fabs(*ix**ix + *iy**iy + *iz**iz - 1);
    double j_diff = fabs(*jx**jx + *jy**jy + *jz**jz - 1);
    double k_diff = fabs(*kx**kx + *ky**ky + *kz**kz - 1);
    double s_diff = fabs(*si**si + *sj**sj + *sk**sk - 1);

    double i_dot_j = *ix**jx + *iy**jy + *iz**jz;
    double i_dot_k = *ix**kx + *iy**ky + *iz**kz;
    double j_dot_k = *jx**kx + *jy**ky + *jz**kz;

    double i_cross_j_x = *iy**jz - *iz**jy;
    double i_cross_j_y = *iz**jx - *ix**jz;
    double i_cross_j_z = *ix**jy - *iy**jx;

    if (i_diff > tolerance || j_diff > tolerance || k_diff > tolerance || s_diff > tolerance) {
        reb_error(sim, "REBOUNDx Error: triaxial_torque: Vectors i, j, k, and s must be unit vectors.\n");
        return 1;
    }
    if (i_dot_j != 0 || i_dot_k != 0 || j_dot_k != 0) {
        reb_error(sim, "REBOUNDx Error: triaxial_torque: Vectors i, j, and k must be mutually orthogonal.\n");
        return 1;
    }
    if (i_cross_j_x != *kx || i_cross_j_y != *ky || i_cross_j_z != *kz){
        reb_error(sim, "REBOUNDx Error: triaxial_torque: The cross-product of vectors i and j must equal vector k.\n");
        return 1;
    }
    if (*Ii < 0 || *Ij < 0 || *Ik < 0){
        reb_error(sim, "REBOUNDx Error: triaxial_torque: Principal moments Ii, Ij, Ik must be >= 0.\n");
        return 1;
    } 
    if (*omega < 0){
        reb_error(sim, "REBOUNDx Error: triaxial_torque: Spin rate omega must be >= 0. Negative spin rates can be created by modifying the s unit vector accordingly.\n");
        return 1;
    }
    if (*R < 0){
        reb_error(sim, "REBOUNDx Error: triaxial_torque: The mean radius R must be >= 0.\n");
        return 1;
    }
    if (*k2 < 0 || *k2 > 1.5){
        reb_error(sim, "REBOUNDx Error: triaxial_torque: k2 must be >= 0 and <= 1.5.\n");
        return 1;
    }

    // else
    return 0;
}

void rebx_triaxial_torque(struct reb_simulation* const sim, struct rebx_operator* const triaxial_torque, const double dt){
    const int _N_real = sim->N - sim->N_var;
	for(int i=0; i<_N_real; i++){
		struct reb_particle* const p = &sim->particles[i];

        // check required params
        // check one first in attempt to speed up computation time
        const double* const R = rebx_get_param(sim->extras, p->ap, "tt_R");
        if (R == NULL) {
            continue;
        }
        
        // then check the remaining parameters all together
        const double* const Ii = rebx_get_param(sim->extras, p->ap, "tt_Ii");
        const double* const Ij = rebx_get_param(sim->extras, p->ap, "tt_Ij");
        const double* const Ik = rebx_get_param(sim->extras, p->ap, "tt_Ik");
        double* const omega = rebx_get_param(sim->extras, p->ap, "tt_omega");
        double* const ix = rebx_get_param(sim->extras, p->ap, "tt_ix");
        double* const iy = rebx_get_param(sim->extras, p->ap, "tt_iy");
        double* const iz = rebx_get_param(sim->extras, p->ap, "tt_iz");
        double* const jx = rebx_get_param(sim->extras, p->ap, "tt_jx");
        double* const jy = rebx_get_param(sim->extras, p->ap, "tt_jy");
        double* const jz = rebx_get_param(sim->extras, p->ap, "tt_jz");
        double* const kx = rebx_get_param(sim->extras, p->ap, "tt_kx");
        double* const ky = rebx_get_param(sim->extras, p->ap, "tt_ky");
        double* const kz = rebx_get_param(sim->extras, p->ap, "tt_kz");
        double* const si = rebx_get_param(sim->extras, p->ap, "tt_si");
        double* const sj = rebx_get_param(sim->extras, p->ap, "tt_sj");
        double* const sk = rebx_get_param(sim->extras, p->ap, "tt_sk");
        const double* const tidal_dt = rebx_get_param(sim->extras, p->ap, "tt_tidal_dt");
        const double* const k2 = rebx_get_param(sim->extras, p->ap, "tt_k2");
        if (Ii==NULL || Ij==NULL || Ik==NULL || omega==NULL || ix==NULL || iy==NULL || iz==NULL || jx==NULL || jy==NULL || jz==NULL 
            || kx==NULL || ky==NULL || kz==NULL || si==NULL || sj==NULL || sk==NULL || tidal_dt==NULL || k2==NULL) {
            continue;
        }

        // check validity of parameters if first timestep
        if (sim->t <= dt){
            if (rebx_validate_params(sim,Ii,Ij,Ik,omega,ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk,tidal_dt,k2,R) == 1) {
                return;
            }
        }
        // printf("time: %.5e\n", sim->t);
        rebx_update_spin_ijk(sim,i,ix,iy,iz,jx,jy,jz,kx,ky,kz,si,sj,sk,omega,*Ii,*Ij,*Ik,*tidal_dt,*k2,*R,dt);
    }
}

/*
// computes time-derivative of vectors i,j,k (components in old ijk basis)
static void rebx_dijk_dt_acc(double ijk_ijk[3][3], double omega_ijk[3], double dijk_dts_ijk[3][3], double sim_dt){
    double omega = sqrt(omega_ijk[0]*omega_ijk[0] + omega_ijk[1]*omega_ijk[1] + omega_ijk[2]*omega_ijk[2]);
    double s_ijk[3] = {omega_ijk[0]/omega, omega_ijk[1]/omega, omega_ijk[2]/omega};
    double dtheta = omega*sim_dt;
    double sin_dtheta = sin(dtheta);
    double cos_dtheta = sqrt(1 - sin_dtheta*sin_dtheta);
    double s_cross_i[3];
    double s_cross_j[3];
    double s_cross_k[3];

    rebx_cross_prod(s_ijk,ijk_ijk[0],s_cross_i);
    rebx_cross_prod(s_ijk,ijk_ijk[1],s_cross_j);
    rebx_cross_prod(s_ijk,ijk_ijk[2],s_cross_k);
    
    // di/dt
    dijk_dts_ijk[0][0] = (sin_dtheta*s_cross_i[0] - (1-cos_dtheta)*ijk_ijk[0][0])/sim_dt;
    dijk_dts_ijk[0][1] = (sin_dtheta*s_cross_i[1] - (1-cos_dtheta)*ijk_ijk[0][1])/sim_dt;
    dijk_dts_ijk[0][2] = (sin_dtheta*s_cross_i[2] - (1-cos_dtheta)*ijk_ijk[0][2])/sim_dt;

    // dj/dt
    dijk_dts_ijk[1][0] = (sin_dtheta*s_cross_j[0] - (1-cos_dtheta)*ijk_ijk[1][0])/sim_dt;
    dijk_dts_ijk[1][1] = (sin_dtheta*s_cross_j[1] - (1-cos_dtheta)*ijk_ijk[1][1])/sim_dt;
    dijk_dts_ijk[1][2] = (sin_dtheta*s_cross_j[2] - (1-cos_dtheta)*ijk_ijk[1][2])/sim_dt;

    // dk/dt
    dijk_dts_ijk[2][0] = (sin_dtheta*s_cross_k[0] - (1-cos_dtheta)*ijk_ijk[2][0])/sim_dt;
    dijk_dts_ijk[2][1] = (sin_dtheta*s_cross_k[1] - (1-cos_dtheta)*ijk_ijk[2][1])/sim_dt;
    dijk_dts_ijk[2][2] = (sin_dtheta*s_cross_k[2] - (1-cos_dtheta)*ijk_ijk[2][2])/sim_dt;
}
*/
