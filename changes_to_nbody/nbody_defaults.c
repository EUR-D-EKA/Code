/*
Copyright (C) 2011  Matthew Arsenault
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

#include "nbody_types.h"
#include "nbody_defaults.h"

const NBodyCtx defaultNBodyCtx =
{
    /* Grr lack of C99 named struct initializers in MSVC */
    /* .eps2            */  0.0,
    /* .theta           */  -1.0,  /* Invalid */

    /* .timestep        */  0.0,
    /* .timeEvolve      */  0.0,
    /* .treeRSize       */  DEFAULT_TREE_ROOT_SIZE,
    /* .sunGCDist       */  DEFAULT_SUN_GC_DISTANCE,

    /* .criterion       */  DEFAULT_CRITERION,
    /* .potentialType   */  EXTERNAL_POTENTIAL_DEFAULT,

    /* .MultiOutput     */  FALSE,
    /* .OutputFreq      */  1,
    
    /* .BestLikeStart   */  DEFAULT_BEST_LIKELIHOOD_START,
    /* .useBestLike     */  DEFAULT_USE_BEST_LIKELIHOOD,
    
    /* .BetaSigma       */  DEFAULT_SIGMA_CUTOFF,
    /* .VelSigma        */  DEFAULT_SIGMA_CUTOFF,
    /* .IterMax         */  DEFAULT_SIGMA_ITER,
    /* .BetaCorrect     */  DEFAULT_DISP_CORRECTION,
    /* .VelCorrect      */  DEFAULT_DISP_CORRECTION,
       
    /* .Nstep_control   */  FALSE,
    /* .Ntsteps         */  0,
    /* .useVelDisp      */  DEFAULT_USE_VEL_DISP,
    /* .useBetaDisp     */  DEFAULT_USE_BETA_DISP,
    /* .useQuad         */  DEFAULT_USE_QUADRUPOLE_MOMENTS,
    /* .allowIncest     */  DEFAULT_ALLOW_INCEST,
    /* .quietErrors     */  DEFAULT_QUIET_ERRORS,


    /* .checkpointT     */  NOBOINC_DEFAULT_CHECKPOINT_PERIOD,
    /* .nStep           */  0,

    /* .pot             */  EMPTY_POTENTIAL,
    /* .LMC             */  FALSE
};

const HistogramParams defaultHistogramParams =
{
    /* .phi         */  histogramPhi,
    /* .theta       */  histogramTheta,
    /* .psi         */  histogramPsi,
    /* .lambdaStart */  histogramlambdaStart,
    /* .lambdaEnd   */  histogramlambdaEnd,
    /* .lambdaBins  */  histogramlambdaBins,
    /* .lambdaStart */  histogrambetaStart,
    /* .lambdaEnd   */  histogrambetaEnd,
    /* .lambdaBins  */  histogrambetaBins
};



