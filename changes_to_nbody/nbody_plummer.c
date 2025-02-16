/* Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
   Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

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
#include "milkyway_util.h"
#include "milkyway_lua.h"
#include "nbody_lua_types.h"
#include "nbody_plummer.h"

/* pickshell: pick a random point on a sphere of specified radius. */
static inline mwvector pickShell(dsfmt_t* dsfmtState, real rad)
{
    real rsq, rsc;
    mwvector vec;

    do                      /* pick point in NDIM-space */
    {
        vec = mwRandomUnitPoint(dsfmtState);
        rsq = mw_sqrv(vec);         /* compute radius squared */
    }
    while (rsq > 1.0);              /* reject if outside sphere */

    rsc = rad / mw_sqrt(rsq);       /* compute scaling factor */
    mw_incmulvs(vec, rsc);          /* rescale to radius given */

    return vec;
}

static inline real plummerRandomR(dsfmt_t* dsfmtState)
{
    real rnd;

    /* returns [0, 1) */
    rnd = (real) dsfmt_genrand_close_open(dsfmtState);

    /* pick r in struct units */
    return 1.0 / mw_sqrt(mw_pow(rnd, -2.0 / 3.0) - 1.0);
}

static inline real plummerSelectFromG(dsfmt_t* dsfmtState)
{
    real x, y;

    do                      /* select from fn g(x) */
    {
        x = mwXrandom(dsfmtState, 0.0, 1.0);      /* for x in range 0:1 */
        y = mwXrandom(dsfmtState, 0.0, 0.1);      /* max of g(x) is 0.092 */
    }   /* using von Neumann tech */
    while (y > sqr(x) * mw_pow(1.0 - sqr(x), 3.5));

    return x;
}

static inline real plummerRandomV(dsfmt_t* dsfmtState, real r)
{
    real x, v;

    x = plummerSelectFromG(dsfmtState);
    v = M_SQRT2 * x / mw_sqrt(mw_sqrt(1.0 + sqr(r)));   /* find v in struct units */

    return v;
}

static inline mwvector plummerBodyPosition(dsfmt_t* dsfmtState, mwvector rshift, real rsc, real r)
{
    mwvector pos;

    pos = pickShell(dsfmtState, rsc * r);  /* pick scaled position */
    mw_incaddv(pos, rshift);               /* move the position */
    pos.x= 125.638, pos.y=683.777, pos.z=-314.013;
    return pos;
}

static inline mwvector plummerBodyVelocity(dsfmt_t* dsfmtState, mwvector vshift, real vsc, real r)
{
    mwvector vel;
    real v;

    v = plummerRandomV(dsfmtState, r);
    vel = pickShell(dsfmtState, vsc * v);   /* pick scaled velocity */
    mw_incaddv(vel, vshift);                /* move the velocity */
    vel.x= -22.6511, vel.y=-139.945, vel.z=41.7121;
    return vel;
}

/* generatePlummer: generate Plummer model initial conditions for test
 * runs, scaled to units such that M = -4E = G = 1 (Henon, Hegge,
 * etc).  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37,
 * 183.
 */
static int nbGeneratePlummerCore(lua_State* luaSt,

                                 dsfmt_t* prng,
                                 unsigned int nbody,
                                 real mass,

                                 mwbool ignore,

                                 mwvector rShift,
                                 mwvector vShift,
                                 real radiusScale)
{
    unsigned int i;
    int table;
    Body b;
    real r, velScale;

    memset(&b, 0, sizeof(b));

    velScale = mw_sqrt(mass / radiusScale);     /* and recip. speed scale */

    b.bodynode.type = BODY(ignore);    /* Same for all in the model */
    b.bodynode.mass = mass / nbody;    /* Mass per particle */

    lua_createtable(luaSt, nbody, 0);
    table = lua_gettop(luaSt);

    for (i = 0; i < nbody; ++i)
    {
        do
        {
            r = plummerRandomR(prng);
            /* FIXME: We should avoid the divide by 0.0 by multiplying
             * the original random number by 0.9999.. but I'm too lazy
             * to change the tests. Same with other models */
        }
        while (isinf(r));
        
        b.bodynode.id = i + 1;
        b.bodynode.pos = plummerBodyPosition(prng, rShift, radiusScale, r);
        b.vel = plummerBodyVelocity(prng, vShift, velScale, r);

        assert(nbPositionValid(b.bodynode.pos));

        pushBody(luaSt, &b);
        lua_rawseti(luaSt, table, i + 1);
    }

    return 1;
}

int nbGeneratePlummer(lua_State* luaSt)
{
    static dsfmt_t* prng;
    static const mwvector* position = NULL;
    static const mwvector* velocity = NULL;
    static mwbool ignore;
    static real mass = 0.0, nbodyf = 0.0, radiusScale = 0.0;

    static const MWNamedArg argTable[] =
        {
            { "nbody",        LUA_TNUMBER,   NULL,          TRUE,  &nbodyf      },
            { "mass",         LUA_TNUMBER,   NULL,          TRUE,  &mass        },
            { "scaleRadius",  LUA_TNUMBER,   NULL,          TRUE,  &radiusScale },
            { "position",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &position    },
            { "velocity",     LUA_TUSERDATA, MWVECTOR_TYPE, TRUE,  &velocity    },
            { "ignore",       LUA_TBOOLEAN,  NULL,          FALSE, &ignore      },
            { "prng",         LUA_TUSERDATA, DSFMT_TYPE,    TRUE,  &prng        },
            END_MW_NAMED_ARG
        };

    if (lua_gettop(luaSt) != 1)
        return luaL_argerror(luaSt, 1, "Expected 1 arguments");

    handleNamedArgumentTable(luaSt, argTable, 1);

    return nbGeneratePlummerCore(luaSt, prng, (unsigned int) nbodyf, mass, ignore,
                                 *position, *velocity, radiusScale);
}

void registerGeneratePlummer(lua_State* luaSt)
{
    lua_register(luaSt, "generatePlummer", nbGeneratePlummer);
}

