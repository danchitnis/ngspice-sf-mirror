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
#include "ngspice/sperror.h"
#include "ngspice/devdefs.h"
#include "ngspice/suffix.h"
#include "ngspice/wordlist.h"
#include "ngspice/cpextern.h"

int
BSIM4checkModel(
BSIM4model *model,
BSIM4instance *here,
CKTcircuit *ckt)
{
struct bsim4SizeDependParam *pParam;
int Fatal_Flag = 0;
FILE *fplog;
wordlist* wl, *wlstart;

    if (cp_getvar("ng_nomodcheck", CP_BOOL, NULL, 0))
        return(0);

    static char modname[BSIZE_SP];
    size_t mlen = strlen(model->BSIM4modName);

    if (mlen < BSIZE_SP) {
        /* Check the model named model->BSIM4modName only once,
           because BSIM4checkModel() is called for each instance. */
        if (!strncmp(modname, model->BSIM4modName, mlen))
            return(0);
        strcpy(modname, model->BSIM4modName);
    }

    pParam = here->pParam;

    wl = wlstart = TMALLOC(wordlist, 1);
    wl->wl_prev = NULL;
    wl->wl_next = NULL;
    wl->wl_word = tprintf("\nChecking parameters for BSIM 4.8 model %s\n", model->BSIM4modName);

    if ((strcmp(model->BSIM4version, "4.8.0")) && (strncmp(model->BSIM4version, "4.80", 4)) && (strncmp(model->BSIM4version, "4.8", 3)) &&
        (strcmp(model->BSIM4version, "4.8.1")) && (strncmp(model->BSIM4version, "4.81", 4)) &&
        (strcmp(model->BSIM4version, "4.8.3")) && (strncmp(model->BSIM4version, "4.83", 4)))
    {
        printf("Warning: This model supports BSIM4 version 4.8\n");
        printf("You specified a wrong version number. Working now with BSIM4.8.3\n");
        wl_append_word(&wl, &wl, tprintf("Warning: This model supports BSIM4 version 4.8\n"));
        wl_append_word(&wl, &wl, tprintf("You specified a wrong version number. Working now with BSIM4.8.3.\n"));
    }

    if ((here->BSIM4rgateMod == 2) || (here->BSIM4rgateMod == 3))
    {   if ((here->BSIM4trnqsMod == 1) || (here->BSIM4acnqsMod == 1))
        {
            wl_append_word(&wl, &wl, tprintf("Warning: You've selected both Rg and charge deficit NQS; select one only.\n"));
        }
    }

    if (model->BSIM4toxe <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Toxe = %g is not positive.\n", model->BSIM4toxe));
        Fatal_Flag = 1;
    }
    if (here->BSIM4toxp <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Toxp = %g is not positive.\n", here->BSIM4toxp));
        Fatal_Flag = 1;
    }
    if (model->BSIM4eot <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: EOT = %g is not positive.\n", model->BSIM4eot));
        Fatal_Flag = 1;
    }
    if(model->BSIM4tempeot <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: TEMPEOT = %g is not positive.\n", model->BSIM4tempeot));
        Fatal_Flag = 1;
    }
    if (model->BSIM4epsrgate < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Epsrgate = %g is not positive.\n", model->BSIM4epsrgate));
        Fatal_Flag = 1;
    }
    if (model->BSIM4epsrsub < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Epsrsub = %g is not positive.\n", model->BSIM4epsrsub));
        Fatal_Flag = 1;
    }
    if (model->BSIM4easub < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Easub = %g is not positive.\n", model->BSIM4easub));
        Fatal_Flag = 1;
    }
    if (model->BSIM4ni0sub <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Easub = %g is not positive.\n", model->BSIM4ni0sub));
        Fatal_Flag = 1;
    }

    if (model->BSIM4toxm <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Toxm = %g is not positive.\n", model->BSIM4toxm));
        Fatal_Flag = 1;
    }

    if (model->BSIM4toxref <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Toxref = %g is not positive.\n", model->BSIM4toxref));
        Fatal_Flag = 1;
    }

    if (pParam->BSIM4lpe0 < -pParam->BSIM4leff)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Lpe0 = %g is less than -Leff.\n",
                    pParam->BSIM4lpe0));
        Fatal_Flag = 1;
    }
    if (model->BSIM4lintnoi > pParam->BSIM4leff/2)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Lintnoi = %g is too large - Leff for noise is negative.\n",
                model->BSIM4lintnoi));
        Fatal_Flag = 1;
    }
    if (pParam->BSIM4lpeb < -pParam->BSIM4leff)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Lpeb = %g is less than -Leff.\n",
                    pParam->BSIM4lpeb));
        Fatal_Flag = 1;
    }
    if (pParam->BSIM4ndep <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Ndep = %g is not positive.\n",
           pParam->BSIM4ndep));
        Fatal_Flag = 1;
    }
    if (pParam->BSIM4phi <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Phi = %g is not positive. Please check Phin and Ndep\n",
                pParam->BSIM4phi));
        wl_append_word(&wl, &wl, tprintf("       Phin = %g  Ndep = %g \n",
                pParam->BSIM4phin, pParam->BSIM4ndep));
        Fatal_Flag = 1;
    }

    if (pParam->BSIM4nsub <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Nsub = %g is not positive.\n",
           pParam->BSIM4nsub));
        Fatal_Flag = 1;
    }
    if (pParam->BSIM4ngate < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Ngate = %g Ngate is not positive.\n",
           pParam->BSIM4ngate));
        Fatal_Flag = 1;
    }
    if (pParam->BSIM4ngate > 1.e25)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Ngate = %g Ngate is too high\n",
           pParam->BSIM4ngate));
        Fatal_Flag = 1;
    }
    if (pParam->BSIM4xj <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Xj = %g is not positive.\n", pParam->BSIM4xj));
        Fatal_Flag = 1;
    }

    if (pParam->BSIM4dvt1 < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Dvt1 = %g is negative.\n", pParam->BSIM4dvt1));
        Fatal_Flag = 1;
    }

    if (pParam->BSIM4dvt1w < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Dvt1w = %g is negative.\n", pParam->BSIM4dvt1w));
        Fatal_Flag = 1;
    }

    if (pParam->BSIM4w0 == -pParam->BSIM4weff)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: (W0 + Weff) = 0 causing divided-by-zero.\n"));
        Fatal_Flag = 1;
        }

    if (pParam->BSIM4dsub < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Dsub = %g is negative.\n", pParam->BSIM4dsub));
        Fatal_Flag = 1;
    }
    if (pParam->BSIM4b1 == -pParam->BSIM4weff)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: (B1 + Weff) = 0 causing divided-by-zero.\n"));
        Fatal_Flag = 1;
    }
    if (here->BSIM4u0temp <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: u0 at current temperature = %g is not positive.\n",
           here->BSIM4u0temp));
        Fatal_Flag = 1;
    }

    if (pParam->BSIM4delta < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Delta = %g is less than zero.\n", pParam->BSIM4delta));
        Fatal_Flag = 1;
    }

    if (here->BSIM4vsattemp <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Vsat at current temperature = %g is not positive.\n",
           here->BSIM4vsattemp));
        Fatal_Flag = 1;
    }

    if (pParam->BSIM4pclm <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Pclm = %g is not positive.\n", pParam->BSIM4pclm));
        Fatal_Flag = 1;
    }

    if (pParam->BSIM4drout < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Drout = %g is negative.\n", pParam->BSIM4drout));
        Fatal_Flag = 1;
    }

    if (here->BSIM4m <= 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: multiplier = %g is not positive.\n", here->BSIM4m));
        Fatal_Flag = 1;
    }

    if (here->BSIM4nf < 1.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Number of finger = %g is smaller than one.\n", here->BSIM4nf));
        Fatal_Flag = 1;
    }

    if((here->BSIM4sa > 0.0) && (here->BSIM4sb > 0.0) &&
    ((here->BSIM4nf == 1.0) || ((here->BSIM4nf > 1.0) && (here->BSIM4sd > 0.0))) )
    {   if (model->BSIM4saref <= 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Fatal: SAref = %g is not positive.\n",model->BSIM4saref));
            Fatal_Flag = 1;
        }
        if (model->BSIM4sbref <= 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Fatal: SBref = %g is not positive.\n",model->BSIM4sbref));
            Fatal_Flag = 1;
        }
    }

    if ((here->BSIM4l + model->BSIM4xl) <= model->BSIM4xgl)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: The parameter xgl must be smaller than Ldrawn+XL.\n"));
        Fatal_Flag = 1;
    }
    if (here->BSIM4ngcon < 1.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: The parameter ngcon cannot be smaller than one.\n"));
        Fatal_Flag = 1;
    }
    if ((here->BSIM4ngcon != 1.0) && (here->BSIM4ngcon != 2.0))
    {   here->BSIM4ngcon = 1.0;
        wl_append_word(&wl, &wl, tprintf("Warning: Ngcon must be equal to one or two; reset to 1.0.\n"));
    }

    if (model->BSIM4gbmin < 1.0e-20)
    {
        wl_append_word(&wl, &wl, tprintf("Warning: Gbmin = %g is too small.\n", model->BSIM4gbmin));
    }

    /* Check saturation parameters */
    if (pParam->BSIM4fprout < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: fprout = %g is negative.\n", pParam->BSIM4fprout));
        Fatal_Flag = 1;
    }
    if (pParam->BSIM4pdits < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: pdits = %g is negative.\n", pParam->BSIM4pdits));
        Fatal_Flag = 1;
    }
    if (model->BSIM4pditsl < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: pditsl = %g is negative.\n", model->BSIM4pditsl));
        Fatal_Flag = 1;
    }

    /* Check gate current parameters */
    if (model->BSIM4igbMod) {
      if (pParam->BSIM4nigbinv <= 0.0)
      {
          wl_append_word(&wl, &wl, tprintf("Fatal: nigbinv = %g is non-positive.\n", pParam->BSIM4nigbinv));
          Fatal_Flag = 1;
      }
      if (pParam->BSIM4nigbacc <= 0.0)
      {
          wl_append_word(&wl, &wl, tprintf("Fatal: nigbacc = %g is non-positive.\n", pParam->BSIM4nigbacc));
          Fatal_Flag = 1;
      }
    }
    if (model->BSIM4igcMod) {
      if (pParam->BSIM4nigc <= 0.0)
      {
                       wl_append_word(&wl, &wl, tprintf("Fatal: nigc = %g is non-positive.\n", pParam->BSIM4nigc));
          Fatal_Flag = 1;
      }
      if (pParam->BSIM4poxedge <= 0.0)
      {
                       wl_append_word(&wl, &wl, tprintf("Fatal: poxedge = %g is non-positive.\n", pParam->BSIM4poxedge));
          Fatal_Flag = 1;
      }
      if (pParam->BSIM4pigcd <= 0.0)
      {
                       wl_append_word(&wl, &wl, tprintf("Fatal: pigcd = %g is non-positive.\n", pParam->BSIM4pigcd));
          Fatal_Flag = 1;
      }
    }

    /* Check capacitance parameters */
    if (pParam->BSIM4clc < 0.0)
    {
       wl_append_word(&wl, &wl, tprintf("Fatal: Clc = %g is negative.\n", pParam->BSIM4clc));
       Fatal_Flag = 1;
    }

    /* Check overlap capacitance parameters */
    if (pParam->BSIM4ckappas < 0.02)
    {
        wl_append_word(&wl, &wl, tprintf("Warning: ckappas = %g is too small.\n", pParam->BSIM4ckappas));
        pParam->BSIM4ckappas = 0.02;
    }
    if (pParam->BSIM4ckappad < 0.02)
    {
        wl_append_word(&wl, &wl, tprintf("Warning: ckappad = %g is too small.\n", pParam->BSIM4ckappad));
        pParam->BSIM4ckappad = 0.02;
    }

    if (model->BSIM4vtss < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Vtss = %g is negative.\n",
            model->BSIM4vtss));
        Fatal_Flag = 1;
    }
    if (model->BSIM4vtsd < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Vtsd = %g is negative.\n",
            model->BSIM4vtsd));
        Fatal_Flag = 1;
    }
    if (model->BSIM4vtssws < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Vtssws = %g is negative.\n",
            model->BSIM4vtssws));
        Fatal_Flag = 1;
    }
    if (model->BSIM4vtsswd < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Vtsswd = %g is negative.\n",
            model->BSIM4vtsswd));
        Fatal_Flag = 1;
    }
    if (model->BSIM4vtsswgs < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Vtsswgs = %g is negative.\n",
            model->BSIM4vtsswgs));
        Fatal_Flag = 1;
    }
    if (model->BSIM4vtsswgd < 0.0)
    {
        wl_append_word(&wl, &wl, tprintf("Fatal: Vtsswgd = %g is negative.\n",
            model->BSIM4vtsswgd));
        Fatal_Flag = 1;
    }


    if (model->BSIM4paramChk ==1)
    {
            /* Check L and W parameters */
        if (pParam->BSIM4leff <= 1.0e-9)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Leff = %g <= 1.0e-9. Recommended Leff >= 1e-8 \n",
                pParam->BSIM4leff));
        }

        if (pParam->BSIM4leffCV <= 1.0e-9)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Leff for CV = %g <= 1.0e-9. Recommended LeffCV >=1e-8 \n",
                pParam->BSIM4leffCV));
        }

            if (pParam->BSIM4weff <= 1.0e-9)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Weff = %g <= 1.0e-9. Recommended Weff >=1e-7 \n",
               pParam->BSIM4weff));
        }

        if (pParam->BSIM4weffCV <= 1.0e-9)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Weff for CV = %g <= 1.0e-9. Recommended WeffCV >= 1e-7 \n",
                pParam->BSIM4weffCV));
        }

        /* Check threshold voltage parameters */
        if (model->BSIM4toxe < 1.0e-10)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Toxe = %g is less than 1A. Recommended Toxe >= 5A\n", model->BSIM4toxe));
        }
        if (here->BSIM4toxp < 1.0e-10)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Toxp = %g is less than 1A. Recommended Toxp >= 5A\n", here->BSIM4toxp));
        }
        if (model->BSIM4toxm < 1.0e-10)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Toxm = %g is less than 1A. Recommended Toxm >= 5A\n", model->BSIM4toxm));
        }

        if (pParam->BSIM4ndep <= 1.0e12)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Ndep = %g may be too small.\n",
                   pParam->BSIM4ndep));
        }
        else if (pParam->BSIM4ndep >= 1.0e21)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Ndep = %g may be too large.\n",
                   pParam->BSIM4ndep));
        }

        if (pParam->BSIM4nsub <= 1.0e14)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Nsub = %g may be too small.\n",
                   pParam->BSIM4nsub));
        }
        else if (pParam->BSIM4nsub >= 1.0e21)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Nsub = %g may be too large.\n",
                   pParam->BSIM4nsub));
        }

        if ((pParam->BSIM4ngate > 0.0) &&
            (pParam->BSIM4ngate <= 1.e18))
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Ngate = %g is less than 1.E18cm^-3.\n",
                   pParam->BSIM4ngate));
        }

            if (pParam->BSIM4dvt0 < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Dvt0 = %g is negative.\n", pParam->BSIM4dvt0));
        }

        if (fabs(1.0e-8 / (pParam->BSIM4w0 + pParam->BSIM4weff)) > 10.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: (W0 + Weff) may be too small.\n"));
            }

            /* Check subthreshold parameters */
        if (pParam->BSIM4nfactor < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Nfactor = %g is negative.\n", pParam->BSIM4nfactor));
        }
        if (pParam->BSIM4cdsc < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Cdsc = %g is negative.\n", pParam->BSIM4cdsc));
        }
        if (pParam->BSIM4cdscd < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Cdscd = %g is negative.\n", pParam->BSIM4cdscd));
        }
            /* Check DIBL parameters */
        if (here->BSIM4eta0 < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Eta0 = %g is negative.\n", here->BSIM4eta0));
        }
        /* Check GIDL parameters */
        if (model->BSIM4gidlclamp >= 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: gidlclamp = %g is zero or positive.\n", model->BSIM4gidlclamp));
        }
            /* Check Abulk parameters */
        if (fabs(1.0e-8 / (pParam->BSIM4b1 + pParam->BSIM4weff)) > 10.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: (B1 + Weff) may be too small.\n"));
        }


            /* Check Saturation parameters */
        if (pParam->BSIM4a2 < 0.01)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: A2 = %g is too small. Set to 0.01.\n",
               pParam->BSIM4a2));
            pParam->BSIM4a2 = 0.01;
        }
        else if (pParam->BSIM4a2 > 1.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: A2 = %g is larger than 1. A2 is set to 1 and A1 is set to 0.\n",
               pParam->BSIM4a2));
            pParam->BSIM4a2 = 1.0;
            pParam->BSIM4a1 = 0.0;
        }

        if (pParam->BSIM4prwg < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Prwg = %g is negative. Set to zero.\n",
                   pParam->BSIM4prwg));
            pParam->BSIM4prwg = 0.0;
        }

        if (pParam->BSIM4rdsw < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Rdsw = %g is negative. Set to zero.\n",
               pParam->BSIM4rdsw));
            pParam->BSIM4rdsw = 0.0;
            pParam->BSIM4rds0 = 0.0;
        }

        if (pParam->BSIM4rds0 < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Rds at current temperature = %g is negative. Set to zero.\n",
               pParam->BSIM4rds0));
            pParam->BSIM4rds0 = 0.0;
        }

        if (pParam->BSIM4rdswmin < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Rdswmin at current temperature = %g is negative. Set to zero.\n",
                   pParam->BSIM4rdswmin));
            pParam->BSIM4rdswmin = 0.0;
        }

        if (pParam->BSIM4pscbe2 <= 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Pscbe2 = %g is not positive.\n", pParam->BSIM4pscbe2));
        }

        if (pParam->BSIM4vsattemp < 1.0e3)
        {
           wl_append_word(&wl, &wl, tprintf("Warning: Vsat at current temperature = %g may be too small.\n", pParam->BSIM4vsattemp));
        }

        if((model->BSIM4lambdaGiven) && (pParam->BSIM4lambda > 0.0) )
        {
            if (pParam->BSIM4lambda > 1.0e-9)
            {
               wl_append_word(&wl, &wl, tprintf("Warning: Lambda = %g may be too large.\n", pParam->BSIM4lambda));
            }
        }

        if((model->BSIM4vtlGiven) && (pParam->BSIM4vtl > 0.0) )
        {
            if (pParam->BSIM4vtl < 6.0e4)
            {
               wl_append_word(&wl, &wl, tprintf("Warning: Thermal velocity vtl = %g may be too small.\n", pParam->BSIM4vtl));
            }

            if (pParam->BSIM4xn < 3.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: back scattering coeff xn = %g is too small. Reset to 3.0 \n", pParam->BSIM4xn));
                pParam->BSIM4xn = 3.0;
            }

            if (model->BSIM4lc < 0.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: back scattering coeff lc = %g is too small. Reset to 0.0\n", model->BSIM4lc));
                pParam->BSIM4lc = 0.0;
            }
        }

        if (pParam->BSIM4pdibl1 < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Pdibl1 = %g is negative.\n", pParam->BSIM4pdibl1));
        }
        if (pParam->BSIM4pdibl2 < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Pdibl2 = %g is negative.\n", pParam->BSIM4pdibl2));
        }

        /* Check stress effect parameters */
        if((here->BSIM4sa > 0.0) && (here->BSIM4sb > 0.0) &&
        ((here->BSIM4nf == 1.0) || ((here->BSIM4nf > 1.0) && (here->BSIM4sd > 0.0))) )
        {   if (model->BSIM4lodk2 <= 0.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: LODK2 = %g is not positive.\n",model->BSIM4lodk2));
            }
            if (model->BSIM4lodeta0 <= 0.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: LODETA0 = %g is not positive.\n",model->BSIM4lodeta0));
            }
        }

        /* Check gate resistance parameters */
        if (here->BSIM4rgateMod == 1)
        {   if (model->BSIM4rshg <= 0.0)
            wl_append_word(&wl, &wl, tprintf("Warning: rshg should be positive for rgateMod = 1.\n"));
        }
        else if (here->BSIM4rgateMod == 2)
        {   if (model->BSIM4rshg <= 0.0)
                wl_append_word(&wl, &wl, tprintf("Warning: rshg <= 0.0 for rgateMod = 2.\n"));
            else if (pParam->BSIM4xrcrg1 <= 0.0)
                     wl_append_word(&wl, &wl, tprintf("Warning: xrcrg1 <= 0.0 for rgateMod = 2.\n"));
        }
        if (here->BSIM4rgateMod == 3)
        {   if (model->BSIM4rshg <= 0.0)
                wl_append_word(&wl, &wl, tprintf("Warning: rshg should be positive for rgateMod = 3.\n"));
            else if (pParam->BSIM4xrcrg1 <= 0.0)
                     wl_append_word(&wl, &wl, tprintf("Warning: xrcrg1 should be positive for rgateMod = 3.\n"));
        }

         /* Check body resistance parameters */

        if (model->BSIM4rbps0 <= 0.0)
        {
                       wl_append_word(&wl, &wl, tprintf("Fatal: RBPS0 = %g is not positive.\n", model->BSIM4rbps0));
            Fatal_Flag = 1;
        }
        if (model->BSIM4rbpd0 <= 0.0)
        {
                       wl_append_word(&wl, &wl, tprintf("Fatal: RBPD0 = %g is not positive.\n", model->BSIM4rbpd0));
            Fatal_Flag = 1;
        }
        if (model->BSIM4rbpbx0 <= 0.0)
        {
                       wl_append_word(&wl, &wl, tprintf("Fatal: RBPBX0 = %g is not positive.\n", model->BSIM4rbpbx0));
            Fatal_Flag = 1;
        }
        if (model->BSIM4rbpby0 <= 0.0)
        {
                       wl_append_word(&wl, &wl, tprintf("Fatal: RBPBY0 = %g is not positive.\n", model->BSIM4rbpby0));
            Fatal_Flag = 1;
        }
        if (model->BSIM4rbdbx0 <= 0.0)
        {
                       wl_append_word(&wl, &wl, tprintf("Fatal: RBDBX0 = %g is not positive.\n", model->BSIM4rbdbx0));
            Fatal_Flag = 1;
        }
        if (model->BSIM4rbdby0 <= 0.0)
        {
                       wl_append_word(&wl, &wl, tprintf("Fatal: RBDBY0 = %g is not positive.\n", model->BSIM4rbdby0));
            Fatal_Flag = 1;
        }
        if (model->BSIM4rbsbx0 <= 0.0)
        {
                       wl_append_word(&wl, &wl, tprintf("Fatal: RBSBX0 = %g is not positive.\n", model->BSIM4rbsbx0));
            Fatal_Flag = 1;
        }
        if (model->BSIM4rbsby0 <= 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Fatal: RBSBY0 = %g is not positive.\n", model->BSIM4rbsby0));
            Fatal_Flag = 1;
        }

        /* Check capacitance parameters */
        if (pParam->BSIM4noff < 0.1)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Noff = %g is too small.\n", pParam->BSIM4noff));
        }

        if (pParam->BSIM4voffcv < -0.5)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Voffcv = %g is too small.\n", pParam->BSIM4voffcv));
        }
        if (pParam->BSIM4moin < 5.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Moin = %g is too small.\n", pParam->BSIM4moin));
        }
        if (pParam->BSIM4moin > 25.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Moin = %g is too large.\n", pParam->BSIM4moin));
        }
        if(model->BSIM4capMod ==2) {
            if (pParam->BSIM4acde < 0.1)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: Acde = %g is too small.\n", pParam->BSIM4acde));
            }
            if (pParam->BSIM4acde > 1.6)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: Acde = %g is too large.\n", pParam->BSIM4acde));
            }
        }

        /* Check overlap capacitance parameters */
        if (model->BSIM4cgdo < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: cgdo = %g is negative. Set to zero.\n", model->BSIM4cgdo));
            model->BSIM4cgdo = 0.0;
            }
        if (model->BSIM4cgso < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: cgso = %g is negative. Set to zero.\n", model->BSIM4cgso));
            model->BSIM4cgso = 0.0;
            }
        if (model->BSIM4cgbo < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: cgbo = %g is negative. Set to zero.\n", model->BSIM4cgbo));
            model->BSIM4cgbo = 0.0;
            }
            if (model->BSIM4tnoiMod == 1){
            wl_append_word(&wl, &wl, tprintf("Warning: TNOIMOD=%d is not supported and may be removed from future version.\n", model->BSIM4tnoiMod));
            }
        if (model->BSIM4idovvdsc <= 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: idovvdsc = %g is zero or negative.\n", model->BSIM4idovvdsc));
        }
        if ((strcmp(model->BSIM4version, "4.8.1")) && (strncmp(model->BSIM4version, "4.81", 4)) && (strncmp(model->BSIM4version, "4.8", 3)) &&
            (strcmp(model->BSIM4version, "4.8.2")) && (strncmp(model->BSIM4version, "4.82", 4)) &&
            (strcmp(model->BSIM4version, "4.8.3")) && (strncmp(model->BSIM4version, "4.83", 4)))
        {
            /* v4.7 */
            if (model->BSIM4tnoiMod == 1 || model->BSIM4tnoiMod == 2) {
                if (model->BSIM4tnoia < 0.0) {
                wl_append_word(&wl, &wl, tprintf("Warning: tnoia = %g is negative. Set to zero.\n", model->BSIM4tnoia));
                model->BSIM4tnoia = 0.0;
                }
                if (model->BSIM4tnoib < 0.0) {
                wl_append_word(&wl, &wl, tprintf("Warning: tnoib = %g is negative. Set to zero.\n", model->BSIM4tnoib));
                model->BSIM4tnoib = 0.0;
                }
                if (model->BSIM4rnoia < 0.0) {
                wl_append_word(&wl, &wl, tprintf("Warning: rnoia = %g is negative. Set to zero.\n", model->BSIM4rnoia));
                model->BSIM4rnoia = 0.0;
                }
                if (model->BSIM4rnoib < 0.0) {
                wl_append_word(&wl, &wl, tprintf("Warning: rnoib = %g is negative. Set to zero.\n", model->BSIM4rnoib));
                model->BSIM4rnoib = 0.0;
                }
            }

            /* v4.7 */
            if (model->BSIM4tnoiMod == 2) {
                    if (model->BSIM4tnoic < 0.0) {

                wl_append_word(&wl, &wl, tprintf("Warning: tnoic = %g is negative. Set to zero.\n", model->BSIM4tnoic));
                model->BSIM4tnoic = 0.0;
                }
                if (model->BSIM4rnoic < 0.0) {

                wl_append_word(&wl, &wl, tprintf("Warning: rnoic = %g is negative. Set to zero.\n", model->BSIM4rnoic));
                model->BSIM4rnoic = 0.0;
                }
            }
        }
        else
        {
            if (model->BSIM4tnoiMod == 1){
                if (model->BSIM4tnoia < 0.0) {

                    wl_append_word(&wl, &wl, tprintf("Warning: tnoia = %g is negative. Set to zero.\n", model->BSIM4tnoia));
                    model->BSIM4tnoia = 0.0;
                }
                if (model->BSIM4tnoib < 0.0) {

                    wl_append_word(&wl, &wl, tprintf("Warning: tnoib = %g is negative. Set to zero.\n", model->BSIM4tnoib));
                    model->BSIM4tnoib = 0.0;
                }
                if (model->BSIM4rnoia < 0.0) {

                    wl_append_word(&wl, &wl, tprintf("Warning: rnoia = %g is negative. Set to zero.\n", model->BSIM4rnoia));
                    model->BSIM4rnoia = 0.0;
                }
                if (model->BSIM4rnoib < 0.0) {

                    wl_append_word(&wl, &wl, tprintf("Warning: rnoib = %g is negative. Set to zero.\n", model->BSIM4rnoib));
                    model->BSIM4rnoib = 0.0;
                }
            }
        }
        /* Limits of Njs and Njd modified in BSIM4.7 */
        if (model->BSIM4SjctEmissionCoeff < 0.1) {
            wl_append_word(&wl, &wl, tprintf("Warning: Njs = %g is less than 0.1. Setting Njs to 0.1.\n", model->BSIM4SjctEmissionCoeff));
            model->BSIM4SjctEmissionCoeff = 0.1;
        }
        else if (model->BSIM4SjctEmissionCoeff < 0.7) {
            wl_append_word(&wl, &wl, tprintf("Warning: Njs = %g is less than 0.7.\n", model->BSIM4SjctEmissionCoeff));
        }
        if (model->BSIM4DjctEmissionCoeff < 0.1) {
            wl_append_word(&wl, &wl, tprintf("Warning: Njd = %g is less than 0.1. Setting Njd to 0.1.\n", model->BSIM4DjctEmissionCoeff));
            model->BSIM4DjctEmissionCoeff = 0.1;
        }
        else if (model->BSIM4DjctEmissionCoeff < 0.7) {
            wl_append_word(&wl, &wl, tprintf("Warning: Njd = %g is less than 0.7.\n", model->BSIM4DjctEmissionCoeff));
        }

        if (model->BSIM4njtsstemp < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Njts = %g is negative at temperature = %g.\n",
                       model->BSIM4njtsstemp, ckt->CKTtemp));
        }
        if (model->BSIM4njtsswstemp < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Njtssw = %g is negative at temperature = %g.\n",
                model->BSIM4njtsswstemp, ckt->CKTtemp));
        }
        if (model->BSIM4njtsswgstemp < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Njtsswg = %g is negative at temperature = %g.\n",
                model->BSIM4njtsswgstemp, ckt->CKTtemp));
        }

        if (model->BSIM4njtsdGiven && model->BSIM4njtsdtemp < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Njtsd = %g is negative at temperature = %g.\n",
                model->BSIM4njtsdtemp, ckt->CKTtemp));
        }
        if (model->BSIM4njtsswdGiven && model->BSIM4njtsswdtemp < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Njtsswd = %g is negative at temperature = %g.\n",
                model->BSIM4njtsswdtemp, ckt->CKTtemp));
        }
        if (model->BSIM4njtsswgdGiven && model->BSIM4njtsswgdtemp < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: Njtsswgd = %g is negative at temperature = %g.\n",
                model->BSIM4njtsswgdtemp, ckt->CKTtemp));
        }

        if (model->BSIM4ntnoi < 0.0)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: ntnoi = %g is negative. Set to zero.\n", model->BSIM4ntnoi));
            model->BSIM4ntnoi = 0.0;
        }

        /* diode model */
        if (model->BSIM4SbulkJctBotGradingCoeff >= 0.99)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: MJS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctBotGradingCoeff));
            model->BSIM4SbulkJctBotGradingCoeff = 0.99;
        }
        if (model->BSIM4SbulkJctSideGradingCoeff >= 0.99)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: MJSWS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctSideGradingCoeff));
            model->BSIM4SbulkJctSideGradingCoeff = 0.99;
        }
        if (model->BSIM4SbulkJctGateSideGradingCoeff >= 0.99)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: MJSWGS = %g is too big. Set to 0.99.\n", model->BSIM4SbulkJctGateSideGradingCoeff));
            model->BSIM4SbulkJctGateSideGradingCoeff = 0.99;
        }

        if (model->BSIM4DbulkJctBotGradingCoeff >= 0.99)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: MJD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctBotGradingCoeff));
            model->BSIM4DbulkJctBotGradingCoeff = 0.99;
        }
        if (model->BSIM4DbulkJctSideGradingCoeff >= 0.99)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: MJSWD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctSideGradingCoeff));
            model->BSIM4DbulkJctSideGradingCoeff = 0.99;
        }
        if (model->BSIM4DbulkJctGateSideGradingCoeff >= 0.99)
        {
            wl_append_word(&wl, &wl, tprintf("Warning: MJSWGD = %g is too big. Set to 0.99.\n", model->BSIM4DbulkJctGateSideGradingCoeff));
            model->BSIM4DbulkJctGateSideGradingCoeff = 0.99;
        }
        if (model->BSIM4wpemod == 1)
        {
            if (model->BSIM4scref <= 0.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: SCREF = %g is not positive. Set to 1e-6.\n", model->BSIM4scref));
                model->BSIM4scref = 1e-6;
            }
            /*Move these checks to temp.c for sceff calculation*/
            /*
            if (here->BSIM4sca < 0.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: SCA = %g is negative. Set to 0.0.\n", here->BSIM4sca);
                here->BSIM4sca = 0.0;
            }
            if (here->BSIM4scb < 0.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: SCB = %g is negative. Set to 0.0.\n", here->BSIM4scb);
                here->BSIM4scb = 0.0;
            }
            if (here->BSIM4scc < 0.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: SCC = %g is negative. Set to 0.0.\n", here->BSIM4scc);
                here->BSIM4scc = 0.0;
            }
            if (here->BSIM4sc < 0.0)
            {
                wl_append_word(&wl, &wl, tprintf("Warning: SC = %g is negative. Set to 0.0.\n", here->BSIM4sc);
                here->BSIM4sc = 0.0;
            }
            */
        }
    }

    if (wlstart->wl_next) {
        if ((fplog = fopen("bsim4.out", "w")) != NULL) {
            while (wlstart) {
                fprintf(fplog, "%s", wlstart->wl_word);
                fprintf(stderr, "%s", wlstart->wl_word);
                wlstart = wlstart->wl_next;
            }
            fclose(fplog);
        }
        else {
            while (wlstart) {
                fprintf(stderr, "%s", wlstart->wl_word);
                wlstart = wlstart->wl_next;
            }
        }
    }

    wl_free(wlstart);
    return(Fatal_Flag);
}

