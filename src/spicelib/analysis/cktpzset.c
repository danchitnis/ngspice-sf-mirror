/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

/* CKTpzSetup(ckt)
 * iterate through all the various
 * pzSetup functions provided for the circuit elements in the
 * given circuit, setup ...
 */

#include "ngspice/ngspice.h"
#include "ngspice/smpdefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "ngspice/sperror.h"

int
CKTpzSetup(CKTcircuit *ckt, int type)
{
    PZAN *job = (PZAN *) ckt->CKTcurJob;

    SMPmatrix *matrix;
    int error;
    int i, solution_col, balance_col;
    int input_pos, input_neg, output_pos, output_neg;

    NIdestroy(ckt);
    error = NIinit(ckt);
    if (error)
	return(error);
    matrix = ckt->CKTmatrix;

    /* Really awful . . . */
    ckt->CKTnumStates = 0;

    for (i = 0; i < DEVmaxnum; i++) {
        if (DEVices[i] && DEVices[i]->DEVpzSetup != NULL && ckt->CKThead[i] != NULL) {
            error = DEVices[i]->DEVpzSetup (matrix, ckt->CKThead[i],
		ckt, &ckt->CKTnumStates);
            if (error != OK)
	        return(error);
        }
    }

    solution_col = 0;
    balance_col = 0;

    input_pos = job->PZin_pos;
    input_neg = job->PZin_neg;

    if (type == PZ_DO_ZEROS) {
	/* Vo/Ii in Y */
	output_pos = job->PZout_pos;
	output_neg = job->PZout_neg;
    } else if (job->PZinput_type == PZ_IN_VOL) {
	/* Vi/Ii in Y */
	output_pos = job->PZin_pos;
	output_neg = job->PZin_neg;
    } else {
	/* Denominator */
	output_pos = 0;
	output_neg = 0;
	input_pos = 0;
	input_neg = 0;
    }

    if (output_pos) {
	solution_col = output_pos;
	if (output_neg)
	    balance_col = output_neg;
    } else {
	solution_col = output_neg;
	SWAP(int, input_pos, input_neg);
    }

    if (input_pos)
	job->PZdrive_pptr = SMPmakeElt(matrix, input_pos, solution_col);
    else
	job->PZdrive_pptr = NULL;

    if (input_neg)
	job->PZdrive_nptr = SMPmakeElt(matrix, input_neg, solution_col);
    else
	job->PZdrive_nptr = NULL;

    job->PZsolution_col = solution_col;
    job->PZbalance_col = balance_col;

    job->PZnumswaps = 1;

#ifdef KLU
    if (matrix->CKTkluMODE)
    {
        fprintf (stdout, "Using KLU as Direct Linear Solver\n") ;

        /* Convert the COO Storage to CSC for KLU and Fill the Binding Table */
        SMPconvertCOOtoCSC (matrix) ;

        /* KLU Pointers Assignment */
        for (i = 0 ; i < DEVmaxnum ; i++)
            if (DEVices [i] && DEVices [i]->DEVbindCSC && ckt->CKThead [i])
                DEVices [i]->DEVbindCSC (ckt->CKThead [i], ckt) ;

        /* ReOrder */
        error = SMPpreOrder (matrix) ;
        if (error) {
            fprintf (stderr, "Error during ReOrdering\n") ;
        }

        /* Conversion from Real Circuit Matrix to Complex Circuit Matrix */
        for (i = 0 ; i < DEVmaxnum ; i++)
            if (DEVices [i] && DEVices [i]->DEVbindCSCComplex && ckt->CKThead [i])
                DEVices [i]->DEVbindCSCComplex (ckt->CKThead [i], ckt) ;

        matrix->SMPkluMatrix->KLUmatrixIsComplex = KLUMatrixComplex ;

        /* Input Pos */
        if ((input_pos > 0) && (solution_col > 0))
        {
            BindElement j, *matched ;

            j.COO = job->PZdrive_pptr ;
            j.CSC = NULL ;
            j.CSC_Complex = NULL ;
            matched = (BindElement *) bsearch (&j, matrix->SMPkluMatrix->KLUmatrixBindStructCOO,
                      matrix->SMPkluMatrix->KLUmatrixLinkedListNZ, sizeof (BindElement), BindCompare) ;
            if (matched == NULL) {
                printf ("Ptr %p not found in BindStruct Table\n", job->PZdrive_pptr) ;
            }
            job->PZdrive_pptr = matched->CSC_Complex ;
        }

        /* Input Neg */
        if ((input_neg > 0) && (solution_col > 0))
        {
            BindElement j, *matched ;

            j.COO = job->PZdrive_nptr ;
            j.CSC = NULL ;
            j.CSC_Complex = NULL ;
            matched = (BindElement *) bsearch (&j, matrix->SMPkluMatrix->KLUmatrixBindStructCOO,
                      matrix->SMPkluMatrix->KLUmatrixLinkedListNZ, sizeof (BindElement), BindCompare) ;
            if (matched == NULL) {
                printf ("Ptr %p not found in BindStruct Table\n", job->PZdrive_nptr) ;
            }
            job->PZdrive_nptr = matched->CSC_Complex ;
        }
    } else {
        fprintf (stdout, "Using SPARSE 1.3 as Direct Linear Solver\n") ;
    }
#endif

    error = NIreinit(ckt);
    if (error)
	return(error);

    return OK;
}

