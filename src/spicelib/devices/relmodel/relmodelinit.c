/**********
Author: Francesco Lannutti - July 2015
**********/

#include "ngspice/config.h"

#include "ngspice/devdefs.h"

#include "relmodelitf.h"
#include "relmodelext.h"
#include "relmodelinit.h"


SPICEdev RELMODELinfo = {
    {
        "RELMODEL",
        "MOSFET Reliability Analysis Addiction Model",

        &RELMODELnSize,
        &RELMODELnSize,
        RELMODELnames,

        &RELMODELpTSize,
        RELMODELpTable,

        &RELMODELmPTSize,
        RELMODELmPTable,

#ifdef XSPICE
        /*----  Fixed by SDB 5.2.2003 to enable XSPICE/tclspice integration  -----*/
                NULL,  /* This is a SPICE device, it has no MIF info data */

                0,     /* This is a SPICE device, it has no MIF info data */
                NULL,  /* This is a SPICE device, it has no MIF info data */

                0,     /* This is a SPICE device, it has no MIF info data */
                NULL,  /* This is a SPICE device, it has no MIF info data */

                0,     /* This is a SPICE device, it has no MIF info data */
                NULL,  /* This is a SPICE device, it has no MIF info data */
                /*---------------------------  End of SDB fix   -------------------------*/
#endif

        DEV_DEFAULT
    },

    NULL,           /* DEVparam       */
    RELMODELmParam, /* DEVmodParam    */
    NULL,           /* DEVload        */
    RELMODELsetup,  /* DEVsetup       */
    NULL,           /* DEVunsetup     */
    NULL,           /* DEVpzSetup     */
    NULL,           /* DEVtemperature */
    NULL,           /* DEVtrunc       */
    NULL,           /* DEVfindBranch  */
    NULL,           /* DEVacLoad      */
    NULL,           /* DEVaccept      */
    NULL,           /* DEVdestroy     */
    NULL,           /* DEVmodDelete   */
    NULL,           /* DEVdelete      */
    NULL,           /* DEVsetic       */
    NULL,           /* DEVask         */
    RELMODELmAsk,   /* DEVmodAsk      */
    NULL,           /* DEVpzLoad      */
    NULL,           /* DEVconvTest    */
    NULL,           /* DEVsenSetup    */
    NULL,           /* DEVsenLoad     */
    NULL,           /* DEVsenUpdate   */
    NULL,           /* DEVsenAcLoad   */
    NULL,           /* DEVsenPrint    */
    NULL,           /* DEVsenTrunc    */
    NULL,           /* DEVdisto       */
    NULL,           /* DEVnoise       */
    NULL,           /* DEVsoaCheck    */
    &RELMODELiSize, /* DEVinstSize    */
    &RELMODELmSize, /* DEVmodSize     */
    NULL            /* DEVreliability */
} ;


SPICEdev *
get_relmodel_info (void)
{
    return &RELMODELinfo ;
}
