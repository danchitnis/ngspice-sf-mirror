/***********************************************************************

 HiSIM (Hiroshima University STARC IGFET Model)
 Copyright (C) 2014 Hiroshima University & STARC

 MODEL NAME : HiSIM_HV 
 ( VERSION : 2  SUBVERSION : 2  REVISION : 0 ) 
 Model Parameter 'VERSION' : 2.20
 FILE : hsmhvnoi.c

 DATE : 2014.6.11

 released by 
                Hiroshima University &
                Semiconductor Technology Academic Research Center (STARC)
***********************************************************************/

/**********************************************************************

The following source code, and all copyrights, trade secrets or other
intellectual property rights in and to the source code in its entirety,
is owned by the Hiroshima University and the STARC organization.

All users need to follow the "HISIM_HV Distribution Statement and
Copyright Notice" attached to HiSIM_HV model.

-----HISIM_HV Distribution Statement and Copyright Notice--------------

Software is distributed as is, completely without warranty or service
support. Hiroshima University or STARC and its employees are not liable
for the condition or performance of the software.

Hiroshima University and STARC own the copyright and grant users a perpetual,
irrevocable, worldwide, non-exclusive, royalty-free license with respect 
to the software as set forth below.   

Hiroshima University and STARC hereby disclaims all implied warranties.

Hiroshima University and STARC grant the users the right to modify, copy,
and redistribute the software and documentation, both within the user's
organization and externally, subject to the following restrictions

1. The users agree not to charge for Hiroshima University and STARC code
itself but may charge for additions, extensions, or support.

2. In any product based on the software, the users agree to acknowledge
Hiroshima University and STARC that developed the software. This
acknowledgment shall appear in the product documentation.

3. The users agree to reproduce any copyright notice which appears on
the software on any copy or modification of such made available
to others."

Toshimasa Asahara, President, Hiroshima University
Mitiko Miura-Mattausch, Professor, Hiroshima University
Katsuhiro Shimohigashi, President&CEO, STARC
June 2008 (revised October 2011) 
*************************************************************************/

#include "ngspice/ngspice.h"
#include "hsmhv2def.h"
#include "ngspice/cktdefs.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/noisedef.h"
#include "ngspice/suffix.h"
#include "ngspice/const.h"  /* jwan */
#include "hsmhv2evalenv.h"
/* #include "hsmhv2macro.h" */

/*
 * HSMHV2noise (mode, operation, firstModel, ckt, data, OnDens)
 *    This routine names and evaluates all of the noise sources
 *    associated with MOSFET's.  It starts with the model *firstModel and
 *    traverses all of its insts.  It then proceeds to any other models
 *    on the linked list.  The total output noise density generated by
 *    all of the MOSFET's is summed with the variable "OnDens".
 */

int HSMHV2noise (
     int mode, int operation,
     GENmodel *inModel,
     CKTcircuit *ckt,
     Ndata *data,
     double *OnDens)
{
  HSMHV2model *model = (HSMHV2model *)inModel;
  HSMHV2instance *here;
  double tempOnoise=0.0 ;
  double tempInoise=0.0 ;
  double noizDens[HSMHV2NSRCS] ;
  double lnNdens[HSMHV2NSRCS] ;
  int i;
  double G =0.0 ;
  double TTEMP = 0.0 ;

  /* define the names of the noise sources */
  static char * HSMHV2nNames[HSMHV2NSRCS] = {
    /* Note that we have to keep the order
       consistent with the index definitions 
       in hsmhvdefs.h */
    ".rd",              /* noise due to rd */
    ".rs",              /* noise due to rs */
    ".id",              /* noise due to id */
    ".1ovf",            /* flicker (1/f) noise */
    ".ign",             /* induced gate noise component at the drain node */
    ""                  /* total transistor noise */
  };
  
  for ( ;model != NULL; model = HSMHV2nextModel(model) ) {
    for ( here = HSMHV2instances(model); here != NULL;
	  here = here->HSMHV2nextInstance ) {
      switch (operation) {
      case N_OPEN:
	/* see if we have to to produce a summary report */
	/* if so, name all the noise generators */
	  
	if (((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0) {
	  switch (mode) {
	  case N_DENS:
	    for ( i = 0; i < HSMHV2NSRCS; i++ ) { 
	      NOISE_ADD_OUTVAR(ckt, data, "onoise.%s%s", here->HSMHV2name, HSMHV2nNames[i]);
	    }
	    break;
	  case INT_NOIZ:
	    for ( i = 0; i < HSMHV2NSRCS; i++ ) {
	      NOISE_ADD_OUTVAR(ckt, data, "onoise_total.%s%s", here->HSMHV2name, HSMHV2nNames[i]);
	      NOISE_ADD_OUTVAR(ckt, data, "inoise_total.%s%s", here->HSMHV2name, HSMHV2nNames[i]);
	    }
	    break;
	  }
	}
	break;
      case N_CALC:
	switch (mode) {
	case N_DENS:

         /* temperature */
         TTEMP = ckt->CKTtemp;
         if ( here->HSMHV2_dtemp_Given ) { TTEMP = TTEMP + here->HSMHV2_dtemp ; }
         TTEMP = TTEMP + *(ckt->CKTstate0 + here->HSMHV2deltemp) ;

         /* rs/rd thermal noise */
	  if ( model->HSMHV2_corsrd == 1 || model->HSMHV2_corsrd == 3 ||  model->HSMHV2_cordrift == 1 ) {
	    NevalSrc(&noizDens[HSMHV2RDNOIZ], NULL,
		     ckt, N_GAIN,
		     here->HSMHV2dNodePrime, here->HSMHV2dNode,
		     0.0);
	    noizDens[HSMHV2RDNOIZ] *= 4 * C_KB * TTEMP * here->HSMHV2drainConductance ;
            lnNdens[HSMHV2RDNOIZ] = log( MAX(noizDens[HSMHV2RDNOIZ],N_MINLOG) );
	    
	    NevalSrc(&noizDens[HSMHV2RSNOIZ], NULL,
		     ckt, N_GAIN,
		     here->HSMHV2sNodePrime, here->HSMHV2sNode,
		     0.0);
	    noizDens[HSMHV2RSNOIZ] *= 4 * C_KB * TTEMP * here->HSMHV2sourceConductance ;
            lnNdens[HSMHV2RSNOIZ] = log( MAX(noizDens[HSMHV2RSNOIZ],N_MINLOG) );
	  } else {
	    noizDens[HSMHV2RDNOIZ] = 0e0 ;
	    lnNdens[HSMHV2RDNOIZ] = N_MINLOG ;
	    noizDens[HSMHV2RSNOIZ] = 0e0 ;
	    lnNdens[HSMHV2RSNOIZ] = N_MINLOG ;
	  }

	  /* channel thermal noise */
	  NevalSrc(&noizDens[HSMHV2IDNOIZ], NULL,
		   ckt, N_GAIN,
		   here->HSMHV2dNodePrime, here->HSMHV2sNodePrime,
		   0.0);
	  switch( model->HSMHV2_noise ) {
	  case 1:
	    /* HiSIMHV model */
	    G = here->HSMHV2_noithrml ;
	    noizDens[HSMHV2IDNOIZ] *= 4 * C_KB * TTEMP * G ;
	    lnNdens[HSMHV2IDNOIZ] = log( MAX(noizDens[HSMHV2IDNOIZ],N_MINLOG) );
	    break;
	  }

	  /* flicker noise */
	  NevalSrc(&noizDens[HSMHV2FLNOIZ], NULL,
		   ckt, N_GAIN,
		   here->HSMHV2dNodePrime, here->HSMHV2sNodePrime, 
		   0.0);
	  switch ( model->HSMHV2_noise ) {
	  case 1:
	    /* HiSIM model */
	    noizDens[HSMHV2FLNOIZ] *= here->HSMHV2_noiflick / pow(data->freq, model->HSMHV2_falph) ; 
	    lnNdens[HSMHV2FLNOIZ] = log(MAX(noizDens[HSMHV2FLNOIZ], N_MINLOG));
	    break;
	  }	  

	  /* induced gate noise */
	  NevalSrc(&noizDens[HSMHV2IGNOIZ], NULL,
		   ckt, N_GAIN, 
		   here->HSMHV2dNodePrime, here->HSMHV2sNodePrime, 
		   0.0);
	  switch ( model->HSMHV2_noise ) {
	  case 1:
	    /* HiSIM model */
	    noizDens[HSMHV2IGNOIZ] *= here->HSMHV2_noiigate * here->HSMHV2_noicross * here->HSMHV2_noicross * data->freq * data->freq;
	    lnNdens[HSMHV2IGNOIZ] = log(MAX(noizDens[HSMHV2IGNOIZ], N_MINLOG));
	    break;
	  }

	  /* total */
	  noizDens[HSMHV2TOTNOIZ] = noizDens[HSMHV2RDNOIZ] + noizDens[HSMHV2RSNOIZ]
	    + noizDens[HSMHV2IDNOIZ] + noizDens[HSMHV2FLNOIZ] + noizDens[HSMHV2IGNOIZ];
	  lnNdens[HSMHV2TOTNOIZ] = log(MAX(noizDens[HSMHV2TOTNOIZ], N_MINLOG));
	  
	  *OnDens += noizDens[HSMHV2TOTNOIZ];
	  
	  if ( data->delFreq == 0.0 ) {
	    /* if we haven't done any previous 
	       integration, we need to initialize our
	       "history" variables.
	    */
	    
	    for ( i = 0; i < HSMHV2NSRCS; i++ ) 
	      here->HSMHV2nVar[LNLSTDENS][i] = lnNdens[i];
	    
	    /* clear out our integration variables
	       if it's the first pass
	    */
	    if (data->freq == ((NOISEAN*) ckt->CKTcurJob)->NstartFreq) {
	      for (i = 0; i < HSMHV2NSRCS; i++) {
		here->HSMHV2nVar[OUTNOIZ][i] = 0.0;
		here->HSMHV2nVar[INNOIZ][i] = 0.0;
	      }
	    }
	  }
	  else {
	    /* data->delFreq != 0.0,
	       we have to integrate.
	    */
	    for ( i = 0; i < HSMHV2NSRCS; i++ ) {
	      if ( i != HSMHV2TOTNOIZ ) {
		tempOnoise = 
		  Nintegrate(noizDens[i], lnNdens[i],
			     here->HSMHV2nVar[LNLSTDENS][i], data);
		tempInoise = 
		  Nintegrate(noizDens[i] * data->GainSqInv, 
			     lnNdens[i] + data->lnGainInv,
			     here->HSMHV2nVar[LNLSTDENS][i] + data->lnGainInv,
			     data);
		here->HSMHV2nVar[LNLSTDENS][i] = lnNdens[i];
		data->outNoiz += tempOnoise;
		data->inNoise += tempInoise;
		if ( ((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0 ) {
		  here->HSMHV2nVar[OUTNOIZ][i] += tempOnoise;
		  here->HSMHV2nVar[OUTNOIZ][HSMHV2TOTNOIZ] += tempOnoise;
		  here->HSMHV2nVar[INNOIZ][i] += tempInoise;
		  here->HSMHV2nVar[INNOIZ][HSMHV2TOTNOIZ] += tempInoise;
		}
	      }
	    }
	  }
	  if ( data->prtSummary ) {
	    for (i = 0; i < HSMHV2NSRCS; i++) {
	      /* print a summary report */
	      data->outpVector[data->outNumber++] = noizDens[i];
	    }
	  }
	  break;
	case INT_NOIZ:
	  /* already calculated, just output */
	  if ( ((NOISEAN*)ckt->CKTcurJob)->NStpsSm != 0 ) {
	    for ( i = 0; i < HSMHV2NSRCS; i++ ) {
	      data->outpVector[data->outNumber++] = here->HSMHV2nVar[OUTNOIZ][i];
	      data->outpVector[data->outNumber++] = here->HSMHV2nVar[INNOIZ][i];
	    }
	  }
	  break;
	}
	break;
      case N_CLOSE:
	/* do nothing, the main calling routine will close */
	return (OK);
	break;   /* the plots */
      }       /* switch (operation) */
    }    /* for here */
  }    /* for model */
  
  return(OK);
}



