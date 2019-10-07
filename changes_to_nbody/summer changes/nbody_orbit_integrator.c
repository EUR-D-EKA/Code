/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

Copyright (c) 2016-2018 Siddhartha Shelton

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "nbody_priv.h"
#include "nbody_orbit_integrator.h"
#include "nbody_potential.h"
#include "nbody_density.c"
#include "nbody_io.h"
#include "nbody_coordinates.h"
#include "nbody_defaults.h"
#include "milkyway_util.h"
#include "nbody_caustic.h"
#include "nbody_bessel.h"

/* Simple orbit integrator in user-defined potential
    Written for BOINC Nbody
    willeb 10 May 2010 */
/* Altered to be consistent with nbody integrator.
 * shelton June 25 2018 */

//added July 20th, 2019

real isotropic_v_dispersion(real a, real v0, real r){
    return (v0 * v0) / 2 * (4 * r * (a + r) *(a+r) * mw_log(1 + r / a) - a * r * (5 * a + 4 * r)) / (a * a * (2 * a + r));
}

mwvector Dynamical_Friction(mwvector pos, mwvector vel, real m, const Potential* pot){
    real pi = 3.1415926;
    real k = 1.428;
    real r = mw_pow((pos.x*pos.x+pos.y*pos.y+pos.z*pos.z), 0.5);
    real lambda = r / 1.6 / k;
    int a = 12;
    int v0 = 73;

    real sigma = isotropic_v_dispersion(a, v0, r);
    real X = mw_pow((vel.x*vel.x+vel.y*vel.y+vel.z*vel.z), 0.5) / pow((2 * sigma), 0.5);

    // density = hernquistSphericalDensity(pot->sph, r) + plummerSpherical(pot->sph, r) + miyamotoNagaiDiskAccel(pot->disk, pos, r)
    //  + hernquistHaloDensity(pot->halo, r) + plummerHaloDensity(pot->halo, r) + NFWMHaloDensity(pot->halo, r) + wilkinsonHalo(pot->halo, r) 
    //  + KVHalo(pot->halo, r);

    real density = plummerSpherical(&(pot->sphere[0]), r) + miyamotoNagaiDiskAccel(&(pot->disk), pos, r) + plummerHaloDensity(&(pot->halo), r);

    real F = -4 * pi * m * m * mw_log(lambda) * density / v0 * 2 * (erf(X) - 2 * X / mw_pow(pi, 0.5) * mw_pow(-1*X, 0.5));
    mwvector result;
    result.x = (F * vel.x / pow((vel.x*vel.x+vel.y*vel.y+vel.z*vel.z), 0.5));
    result.y = (F * vel.y / pow((vel.x*vel.x+vel.y*vel.y+vel.z*vel.z), 0.5));
    result.z = (F * vel.z / pow((vel.x*vel.x+vel.y*vel.y+vel.z*vel.z), 0.5));
    return result;
}
//

void nbReverseOrbit(mwvector* finalPos,
                    mwvector* finalVel,
                    const Potential* pot,
                    mwvector pos,
                    mwvector vel,
                    real tstop,
                    real dt)
{
    mwvector acc, v, x, mw_shift_x, mw_shift_v, temp_acc, temp_v, temp_x;
    real t;
    real dt_half = dt / 2.0;
    const real LMC_mass = mw_pow(10.00, 11.00)/222288.47;
    
    // Set the initial conditions
    x = pos;
    v = vel;
    mw_shift_x.x = 0;
    mw_shift_x.y = 0;
    mw_shift_x.z = 0;

    mw_shift_v.x = 0;
    mw_shift_v.y = 0;
    mw_shift_v.z = 0;
    mw_incnegv(v);

    
    // Get the initial acceleration
    //acc = nbExtAcceleration(pot, x);
    acc.x = X(nbExtAcceleration(pot, x)) - X(Dynamical_Friction(x, v, LMC_mass, pot));
    acc.y = Y(nbExtAcceleration(pot, x)) - Y(Dynamical_Friction(x, v, LMC_mass, pot));
    acc.z = Z(nbExtAcceleration(pot, x)) - Z(Dynamical_Friction(x, v, LMC_mass, pot));

    for (t = 0; t <= tstop; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt);
        
        // Compute the new acceleration
        //acc = nbExtAcceleration(pot, x);
        mwvector pos_sub;
        pos_sub.x = X(x) - X(mw_shift_x);
        pos_sub.y = Y(x) - Y(mw_shift_x);
        pos_sub.x = Z(x) - Z(mw_shift_x);

        mwvector vel_sub;
        vel_sub.x = X(v) - X(mw_shift_v);
        vel_sub.y = Y(v) - Y(mw_shift_v);
        vel_sub.x = Z(v) - Z(mw_shift_v);

        acc.x = X(nbExtAcceleration(pot, pos_sub)) - X(Dynamical_Friction(pos_sub, vel_sub, LMC_mass, pot));
        acc.y = Y(nbExtAcceleration(pot, pos_sub)) - Y(Dynamical_Friction(pos_sub, vel_sub, LMC_mass, pot));
        acc.z = Z(nbExtAcceleration(pot, pos_sub)) - Z(Dynamical_Friction(pos_sub, vel_sub, LMC_mass, pot));

        // Update mw_shift
        temp_acc.x = X(nbExtAcceleration(pot, x)) - X(Dynamical_Friction(x, v, LMC_mass, pot));
        temp_acc.y = Y(nbExtAcceleration(pot, x)) - Y(Dynamical_Friction(x, v, LMC_mass, pot));
        temp_acc.z = Z(nbExtAcceleration(pot, x)) - Z(Dynamical_Friction(x, v, LMC_mass, pot));
        temp_v = v;
        temp_x = x;

        mw_incaddv_s(temp_v, temp_acc, dt);
        mw_incaddv_s(temp_x, temp_v, dt);

        mw_shift_x.x = (X(temp_x) + X(x))/2.0;
        mw_shift_x.y = (Y(temp_x) + Y(x))/2.0;
        mw_shift_x.z = (Z(temp_x) + Z(x))/2.0;


        mw_shift_v.x = (X(temp_v) + X(v))/2.0;
        mw_shift_v.y = (Y(temp_v) + Y(v))/2.0;
        mw_shift_v.z = (Z(temp_v) + Z(v))/2.0;

        mw_incaddv_s(v, acc, dt_half);
    }
    
    
    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);
    
    *finalPos = x;
    *finalVel = v;
}

void nbPrintReverseOrbit(mwvector* finalPos,
                         mwvector* finalVel,
                         const Potential* pot,
                         mwvector pos,
                         mwvector vel,
                         real tstop,
                         real tstopforward,
                         real dt)
{
    mwvector acc, v, x;
    mwvector v_for, x_for;
    mwvector lbr;
    real t;
    real dt_half = dt / 2.0;
    const real LMC_mass = mw_pow(10.00, 11.00)/222288.47;

    // Set the initial conditions
    x = pos;
    v = vel;
    x_for = pos;
    v_for = vel;
    
    mw_incnegv(v);

    // Get the initial acceleration
    //acc = nbExtAcceleration(pot, x);
    acc.x = X(nbExtAcceleration(pot, x)) - X(Dynamical_Friction(x, v, LMC_mass, pot));
    acc.y = X(nbExtAcceleration(pot, x)) - Y(Dynamical_Friction(x, v, LMC_mass, pot));
    acc.z = X(nbExtAcceleration(pot, x)) - Z(Dynamical_Friction(x, v, LMC_mass, pot));

    FILE * fp;
    fp = fopen("reverse_orbit.out", "w");
    // Loop through time
    for (t = 0; t <= tstop; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v, acc, dt_half);
        mw_incaddv_s(x, v, dt);
        
        
        // Compute the new acceleration
        //acc = nbExtAcceleration(pot, x);
        acc.x = X(nbExtAcceleration(pot, x)) - X(Dynamical_Friction(x, v, LMC_mass, pot));
        acc.y = Y(nbExtAcceleration(pot, x)) - Y(Dynamical_Friction(x, v, LMC_mass, pot));
        acc.z = Z(nbExtAcceleration(pot, x)) - Z(Dynamical_Friction(x, v, LMC_mass, pot));
        mw_incaddv_s(v, acc, dt_half);
        
        lbr = cartesianToLbr(x, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", X(x), Y(x), Z(x), X(lbr), Y(lbr), Z(lbr), X(v), Y(v), Z(v));
    }

    fclose(fp);
    fp = fopen("forward_orbit.out", "w");
    acc = nbExtAcceleration(pot, x_for);
    for (t = 0; t <= tstopforward; t += dt)
    {
        // Update the velocities and positions
        mw_incaddv_s(v_for, acc, dt_half);
        mw_incaddv_s(x_for, v_for, dt);
        
        
        // Compute the new acceleration
        acc = nbExtAcceleration(pot, x_for);
        mw_incaddv_s(v_for, acc, dt_half);
        
        lbr = cartesianToLbr(x_for, DEFAULT_SUN_GC_DISTANCE);
        fprintf(fp, "%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n", X(x_for), Y(x_for), Z(x_for), X(lbr), Y(lbr), Z(lbr), X(v_for), Y(v_for), Z(v_for));
    }
    fclose(fp);
    
    /* Report the final values (don't forget to reverse the velocities) */
    mw_incnegv(v);

    *finalPos = x;
    *finalVel = v;
}
