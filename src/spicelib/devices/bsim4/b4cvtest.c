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
#include "ngspice/cktdefs.h"
#include "bsim4def.h"
#include "ngspice/trandefs.h"
#include "ngspice/const.h"
#include "ngspice/devdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"

int
BSIM4convTest(
GENmodel *inModel,
CKTcircuit *ckt)
{
BSIM4model *model = (BSIM4model*)inModel;
BSIM4instance *here;
double delvbd, delvbs, delvds, delvgd, delvgs;
double delvdbd, delvsbs;
double delvbd_jct, delvbs_jct;
double vds, vgs, vgd, vgdo, vbs, vbd;
double vdbd, vdbs, vsbs;
double cbhat, cdhat, Idtot, Ibtot;
double vses, vdes, vdedo, delvses, delvded, delvdes;
double Isestot, cseshat, Idedtot, cdedhat;
double Igstot, cgshat, Igdtot, cgdhat, Igbtot, cgbhat;
double tol0, tol1, tol2, tol3, tol4, tol5, tol6;

    for (; model != NULL; model = BSIM4nextModel(model))
    {    for (here = BSIM4instances(model); here != NULL ;
              here=BSIM4nextInstance(here)) 
         {    vds = model->BSIM4type
                  * (*(ckt->CKTrhsOld + here->BSIM4dNodePrime)
                  - *(ckt->CKTrhsOld + here->BSIM4sNodePrime));
              vgs = model->BSIM4type
                  * (*(ckt->CKTrhsOld + here->BSIM4gNodePrime)
                  - *(ckt->CKTrhsOld + here->BSIM4sNodePrime));
              vbs = model->BSIM4type
                  * (*(ckt->CKTrhsOld + here->BSIM4bNodePrime)
                  - *(ckt->CKTrhsOld + here->BSIM4sNodePrime));
              vdbs = model->BSIM4type
                   * (*(ckt->CKTrhsOld + here->BSIM4dbNode)
                   - *(ckt->CKTrhsOld + here->BSIM4sNodePrime));
              vsbs = model->BSIM4type
                   * (*(ckt->CKTrhsOld + here->BSIM4sbNode)
                   - *(ckt->CKTrhsOld + here->BSIM4sNodePrime));
              vses = model->BSIM4type
                   * (*(ckt->CKTrhsOld + here->BSIM4sNode)
                   - *(ckt->CKTrhsOld + here->BSIM4sNodePrime));
              vdes = model->BSIM4type
                   * (*(ckt->CKTrhsOld + here->BSIM4dNode)
                   - *(ckt->CKTrhsOld + here->BSIM4sNodePrime));

              vgdo = *(ckt->CKTstate0 + here->BSIM4vgs)
                    - *(ckt->CKTstate0 + here->BSIM4vds);
              vbd = vbs - vds;
              vdbd = vdbs - vds;
              vgd = vgs - vds;

              delvbd = vbd - *(ckt->CKTstate0 + here->BSIM4vbd);
              delvdbd = vdbd - *(ckt->CKTstate0 + here->BSIM4vdbd);
              delvgd = vgd - vgdo;

              delvds = vds - *(ckt->CKTstate0 + here->BSIM4vds);
              delvgs = vgs - *(ckt->CKTstate0 + here->BSIM4vgs);
              delvbs = vbs - *(ckt->CKTstate0 + here->BSIM4vbs);
              delvsbs = vsbs - *(ckt->CKTstate0 + here->BSIM4vsbs);

              delvses = vses - (*(ckt->CKTstate0 + here->BSIM4vses));
              vdedo = *(ckt->CKTstate0 + here->BSIM4vdes)
                    - *(ckt->CKTstate0 + here->BSIM4vds);
              delvdes = vdes - *(ckt->CKTstate0 + here->BSIM4vdes);
              delvded = vdes - vds - vdedo;

              delvbd_jct = (!here->BSIM4rbodyMod) ? delvbd : delvdbd;
              delvbs_jct = (!here->BSIM4rbodyMod) ? delvbs : delvsbs;

              if (here->BSIM4mode >= 0)
              {   Idtot = here->BSIM4cd + here->BSIM4csub - here->BSIM4cbd
                        + here->BSIM4Igidl;
                  cdhat = Idtot - here->BSIM4gbd * delvbd_jct
                        + (here->BSIM4gmbs + here->BSIM4gbbs + here->BSIM4ggidlb) * delvbs
                        + (here->BSIM4gm + here->BSIM4gbgs + here->BSIM4ggidlg) * delvgs
                        + (here->BSIM4gds + here->BSIM4gbds + here->BSIM4ggidld) * delvds;

                  Igstot = here->BSIM4Igs + here->BSIM4Igcs;
                  cgshat = Igstot + (here->BSIM4gIgsg + here->BSIM4gIgcsg) * delvgs
                         + here->BSIM4gIgcsd * delvds + here->BSIM4gIgcsb * delvbs;

                  Igdtot = here->BSIM4Igd + here->BSIM4Igcd;
                  cgdhat = Igdtot + here->BSIM4gIgdg * delvgd + here->BSIM4gIgcdg * delvgs
                         + here->BSIM4gIgcdd * delvds + here->BSIM4gIgcdb * delvbs;

                  Igbtot = here->BSIM4Igb;
                  cgbhat = here->BSIM4Igb + here->BSIM4gIgbg * delvgs + here->BSIM4gIgbd
                         * delvds + here->BSIM4gIgbb * delvbs;
              }
              else
              {   Idtot = here->BSIM4cd + here->BSIM4cbd - here->BSIM4Igidl; /* bugfix */
                  cdhat = Idtot + here->BSIM4gbd * delvbd_jct + here->BSIM4gmbs
                        * delvbd + here->BSIM4gm * delvgd
                        - (here->BSIM4gds + here->BSIM4ggidls) * delvds
                        - here->BSIM4ggidlg * delvgs - here->BSIM4ggidlb * delvbs;

                  Igstot = here->BSIM4Igs + here->BSIM4Igcd;
                  cgshat = Igstot + here->BSIM4gIgsg * delvgs + here->BSIM4gIgcdg * delvgd
                         - here->BSIM4gIgcdd * delvds + here->BSIM4gIgcdb * delvbd;

                  Igdtot = here->BSIM4Igd + here->BSIM4Igcs;
                  cgdhat = Igdtot + (here->BSIM4gIgdg + here->BSIM4gIgcsg) * delvgd
                         - here->BSIM4gIgcsd * delvds + here->BSIM4gIgcsb * delvbd;

                  Igbtot = here->BSIM4Igb;
                  cgbhat = here->BSIM4Igb + here->BSIM4gIgbg * delvgd - here->BSIM4gIgbd
                         * delvds + here->BSIM4gIgbb * delvbd;
              }

              Isestot = here->BSIM4gstot * (*(ckt->CKTstate0 + here->BSIM4vses));
              cseshat = Isestot + here->BSIM4gstot * delvses
                      + here->BSIM4gstotd * delvds + here->BSIM4gstotg * delvgs
                      + here->BSIM4gstotb * delvbs;

              Idedtot = here->BSIM4gdtot * vdedo;
              cdedhat = Idedtot + here->BSIM4gdtot * delvded
                      + here->BSIM4gdtotd * delvds + here->BSIM4gdtotg * delvgs
                      + here->BSIM4gdtotb * delvbs;

              /*
               *  Check convergence
               */

              if ((here->BSIM4off == 0)  || (!(ckt->CKTmode & MODEINITFIX)))
              {   tol0 = ckt->CKTreltol * MAX(fabs(cdhat), fabs(Idtot))
                       + ckt->CKTabstol;
                  tol1 = ckt->CKTreltol * MAX(fabs(cseshat), fabs(Isestot))
                       + ckt->CKTabstol;
                  tol2 = ckt->CKTreltol * MAX(fabs(cdedhat), fabs(Idedtot))
                       + ckt->CKTabstol;
                  tol3 = ckt->CKTreltol * MAX(fabs(cgshat), fabs(Igstot))
                       + ckt->CKTabstol;
                  tol4 = ckt->CKTreltol * MAX(fabs(cgdhat), fabs(Igdtot))
                       + ckt->CKTabstol;
                  tol5 = ckt->CKTreltol * MAX(fabs(cgbhat), fabs(Igbtot))
                       + ckt->CKTabstol;

                  if ((fabs(cdhat - Idtot) >= tol0) || (fabs(cseshat - Isestot) >= tol1)
                      || (fabs(cdedhat - Idedtot) >= tol2))
                  {   ckt->CKTnoncon++;
                      return(OK);
                  }

                  if ((fabs(cgshat - Igstot) >= tol3) || (fabs(cgdhat - Igdtot) >= tol4)
                      || (fabs(cgbhat - Igbtot) >= tol5))
                  {   ckt->CKTnoncon++;
                      return(OK);
                  }

                  Ibtot = here->BSIM4cbs + here->BSIM4cbd
                        - here->BSIM4Igidl - here->BSIM4Igisl - here->BSIM4csub;
                  if (here->BSIM4mode >= 0)
                  {   cbhat = Ibtot + here->BSIM4gbd * delvbd_jct
                            + here->BSIM4gbs * delvbs_jct - (here->BSIM4gbbs + here->BSIM4ggidlb)
                            * delvbs - (here->BSIM4gbgs + here->BSIM4ggidlg) * delvgs
                            - (here->BSIM4gbds + here->BSIM4ggidld) * delvds
                            - here->BSIM4ggislg * delvgd - here->BSIM4ggislb* delvbd + here->BSIM4ggisls * delvds ;
                  }
                  else
                  {   cbhat = Ibtot + here->BSIM4gbs * delvbs_jct + here->BSIM4gbd
                                 * delvbd_jct - (here->BSIM4gbbs + here->BSIM4ggislb) * delvbd
                                 - (here->BSIM4gbgs + here->BSIM4ggislg) * delvgd
                                 + (here->BSIM4gbds + here->BSIM4ggisld - here->BSIM4ggidls) * delvds
                                 - here->BSIM4ggidlg * delvgs - here->BSIM4ggidlb * delvbs;
                  }
                  tol6 = ckt->CKTreltol * MAX(fabs(cbhat),
                  fabs(Ibtot)) + ckt->CKTabstol;
                  if (fabs(cbhat - Ibtot) > tol6)
                  {   ckt->CKTnoncon++;
                      return(OK);
                  }
              }
         }
    }
    return(OK);
}
