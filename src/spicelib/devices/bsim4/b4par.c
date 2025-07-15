/* ******************************************************************************
   *  BSIM4 4.8.3 released on 05/19/2025                                        *
   *  BSIM4 Model Equations                                                     *
   ******************************************************************************

   ******************************************************************************
   *  Copyright (c) 2025 University of California                               *
   *                                                                            *
   *  Project Directors: Prof. Sayeef Salahuddin and Prof. Chenming Hu          *
   *  Developers list: https://www.bsim.berkeley.edu/models/bsim4/auth_bsim4/   *
   ******************************************************************************/

/*
Licensed under Educational Community License, Version 2.0 (the "License"); you may
not use this file except in compliance with the License. You may obtain a copy of the license at
http://opensource.org/licenses/ECL-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations
under the License.

BSIM4 model is supported by the members of Silicon Integration Initiative's Compact Model Coalition. A link to the most recent version of this
standard can be found at: http://www.si2.org/cmc
*/

#include "ngspice/ngspice.h"
#include "ngspice/ifsim.h"
#include "bsim4def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/fteext.h"

int
BSIM4param(
int param,
IFvalue *value,
GENinstance *inst,
IFvalue *select)
{
    double scale;

    BSIM4instance *here = (BSIM4instance*)inst;

    NG_IGNORE(select);

    if (!cp_getvar("scale", CP_REAL, &scale, 0))
        scale = 1;

    switch(param)
    {   case BSIM4_W:
	here->BSIM4w = value->rValue * scale;
            here->BSIM4wGiven = TRUE;
            break;
        case BSIM4_L:
            here->BSIM4l = value->rValue * scale;
            here->BSIM4lGiven = TRUE;
            break;
        case BSIM4_M:
            here->BSIM4m = value->rValue;
            here->BSIM4mGiven = TRUE;
            break;
        case BSIM4_MULT_I:
            here->BSIM4mult_i = value->rValue;
            here->BSIM4mult_iGiven = TRUE;
            break;
        case BSIM4_MULT_Q:
            here->BSIM4mult_q = value->rValue;
            here->BSIM4mult_qGiven = TRUE;
            break;
        case BSIM4_MULT_FN:
            here->BSIM4mult_fn = value->rValue;
            here->BSIM4mult_fnGiven = TRUE;
            break;
        case BSIM4_NF:
            here->BSIM4nf = value->rValue;
            here->BSIM4nfGiven = TRUE;
            break;
        case BSIM4_MIN:
            here->BSIM4min = value->iValue;
            here->BSIM4minGiven = TRUE;
            break;
        case BSIM4_AS:
            here->BSIM4sourceArea = value->rValue * scale * scale;
            here->BSIM4sourceAreaGiven = TRUE;
            break;
        case BSIM4_AD:
            here->BSIM4drainArea = value->rValue * scale * scale;
            here->BSIM4drainAreaGiven = TRUE;
            break;
        case BSIM4_PS:
            here->BSIM4sourcePerimeter = value->rValue * scale;
            here->BSIM4sourcePerimeterGiven = TRUE;
            break;
        case BSIM4_PD:
            here->BSIM4drainPerimeter = value->rValue * scale;
            here->BSIM4drainPerimeterGiven = TRUE;
            break;
        case BSIM4_NRS:
            here->BSIM4sourceSquares = value->rValue;
            here->BSIM4sourceSquaresGiven = TRUE;
            break;
        case BSIM4_NRD:
            here->BSIM4drainSquares = value->rValue;
            here->BSIM4drainSquaresGiven = TRUE;
            break;
        case BSIM4_OFF:
            here->BSIM4off = value->iValue;
            break;
        case BSIM4_SA:
            here->BSIM4sa = value->rValue;
            here->BSIM4saGiven = TRUE;
            break;
        case BSIM4_SB:
            here->BSIM4sb = value->rValue;
            here->BSIM4sbGiven = TRUE;
            break;
        case BSIM4_SD:
            here->BSIM4sd = value->rValue;
            here->BSIM4sdGiven = TRUE;
            break;
        case BSIM4_SCA:
            here->BSIM4sca = value->rValue;
            here->BSIM4scaGiven = TRUE;
            break;
        case BSIM4_SCB:
            here->BSIM4scb = value->rValue;
            here->BSIM4scbGiven = TRUE;
            break;
        case BSIM4_SCC:
            here->BSIM4scc = value->rValue;
            here->BSIM4sccGiven = TRUE;
            break;
        case BSIM4_SC:
            here->BSIM4sc = value->rValue;
            here->BSIM4scGiven = TRUE;
            break;
        case BSIM4_RBSB:
            here->BSIM4rbsb = value->rValue;
            here->BSIM4rbsbGiven = TRUE;
            break;
        case BSIM4_RBDB:
            here->BSIM4rbdb = value->rValue;
            here->BSIM4rbdbGiven = TRUE;
            break;
        case BSIM4_RBPB:
            here->BSIM4rbpb = value->rValue;
            here->BSIM4rbpbGiven = TRUE;
            break;
        case BSIM4_RBPS:
            here->BSIM4rbps = value->rValue;
            here->BSIM4rbpsGiven = TRUE;
            break;
        case BSIM4_RBPD:
            here->BSIM4rbpd = value->rValue;
            here->BSIM4rbpdGiven = TRUE;
            break;
        case BSIM4_DELVTO:
            here->BSIM4delvto = value->rValue;
            here->BSIM4delvtoGiven = TRUE;
            break;
        case BSIM4_XGW:
            here->BSIM4xgw = value->rValue;
            here->BSIM4xgwGiven = TRUE;
            break;
        case BSIM4_NGCON:
            here->BSIM4ngcon = value->rValue;
            here->BSIM4ngconGiven = TRUE;
            break;
        case BSIM4_TRNQSMOD:
            here->BSIM4trnqsMod = value->iValue;
            here->BSIM4trnqsModGiven = TRUE;
            break;
        case BSIM4_ACNQSMOD:
            here->BSIM4acnqsMod = value->iValue;
            here->BSIM4acnqsModGiven = TRUE;
            break;
        case BSIM4_RBODYMOD:
            here->BSIM4rbodyMod = value->iValue;
            here->BSIM4rbodyModGiven = TRUE;
            break;
        case BSIM4_RGATEMOD:
            here->BSIM4rgateMod = value->iValue;
            here->BSIM4rgateModGiven = TRUE;
            break;
        case BSIM4_GEOMOD:
            here->BSIM4geoMod = value->iValue;
            here->BSIM4geoModGiven = TRUE;
            break;
        case BSIM4_RGEOMOD:
            here->BSIM4rgeoMod = value->iValue;
            here->BSIM4rgeoModGiven = TRUE;
            break;
        case BSIM4_IC_VDS:
            here->BSIM4icVDS = value->rValue;
            here->BSIM4icVDSGiven = TRUE;
            break;
        case BSIM4_IC_VGS:
            here->BSIM4icVGS = value->rValue;
            here->BSIM4icVGSGiven = TRUE;
            break;
        case BSIM4_IC_VBS:
            here->BSIM4icVBS = value->rValue;
            here->BSIM4icVBSGiven = TRUE;
            break;
        case BSIM4_IC:
            /* FALLTHROUGH added to suppress GCC warning due to
             * -Wimplicit-fallthrough flag */
            switch(value->v.numValue)
            {   case 3:
                    here->BSIM4icVBS = *(value->v.vec.rVec+2);
                    here->BSIM4icVBSGiven = TRUE;
                    /* FALLTHROUGH */
                case 2:
                    here->BSIM4icVGS = *(value->v.vec.rVec+1);
                    here->BSIM4icVGSGiven = TRUE;
                    /* FALLTHROUGH */
                case 1:
                    here->BSIM4icVDS = *(value->v.vec.rVec);
                    here->BSIM4icVDSGiven = TRUE;
                    break;
                default:
                    return(E_BADPARM);
            }
            break;
        default:
            return(E_BADPARM);
    }
    return(OK);
}
