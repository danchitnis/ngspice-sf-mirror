/* Copyright 2013 - 2019 Holger Vogt
 *
 * Modified BSD license
 */

/* For comments and explanations see sharedspice.h */

/*******************/
/*   Defines       */
/*******************/

#ifdef _MSC_VER
#define STDIN_FILENO    0
#define STDOUT_FILENO   1
#define STDERR_FILENO   2
#endif

/* If a calling function has high latency times during printing,
   causing memory access errors, you may undef the following line.
   Printing messages are assembled in a wordlist, and sent to the caller
   via a new thread. Delays may occur. */
#define low_latency

/************* About threads in sharedspice.c *************************
   If the calling (main) thread loads a circuit, the .control section
   commands in the input file are executed immediately by the calling
   thread after the ciruit has been parsed and loaded.
   Command bg_run from the calling thread then immediately starts the
   background thread (id. tid) that issues the 'run' command to
   start the simulation in this thread. The main thread returns to the
   caller. .control commands typically are executed prematurely before
   bg_run has returned.
   If the flag 'set controlswait' is given in the .control section,
   all commands following are assembled in the wordlist 'shcontrols',
   a new thread is started (id: tid2) and suspended immediately. Only
   when the background thread tid (and thus the simulation) is ready,
   the tid2 thread is released and the .control commands are executed.
   Before a repeated 'bg_run' is given, or after a 'reset', the command
   'bg_ctrl' has to be sent by the caller to re-start and suspend the
   thread tid2, using the still existing shcontrols.
*/

/**********************************************************************/
/*              Header files for C functions                          */
/**********************************************************************/

#include <stdio.h>
#include <string.h>
#include <setjmp.h>

/* workaround since fputs, putc are replaced by sh_fputs,
sh_putc, through redefinition in ngspice.h */

int myputs(const char* inp, FILE* f);
int myputc(int inp, FILE* f);
int myfputc(int inp, FILE* f);

int
myputs(const char* inp, FILE* f)
{
    return fputs(inp, f);
}

int
myputc(int inp, FILE* f)
{
    return putc(inp, f);
}

int
myfputc(int inp, FILE* f)
{
    return fputc(inp, f);
}

#if defined(__MINGW32__) || defined(_MSC_VER)
#include <windows.h>
#endif

#include "ngspice/ngspice.h"
#include "misc/misc_time.h"
#include "ngspice/randnumb.h"

/*Use Windows threads if on W32 without pthreads*/
#ifndef HAVE_LIBPTHREAD

#if defined(__MINGW32__) || defined(_MSC_VER)
//#if defined(_MSC_VER)
#ifdef SRW
#define mutex_lock(a) AcquireSRWLockExclusive(a)
#define mutex_unlock(a) ReleaseSRWLockExclusive(a)
typedef SRWLOCK mutexType;
#else
#define mutex_lock(a) EnterCriticalSection(a)
#define mutex_unlock(a) LeaveCriticalSection(a)
typedef CRITICAL_SECTION mutexType;
#endif
#define thread_self() GetCurrentThread()
#define threadid_self() GetCurrentThreadId()
typedef HANDLE threadId_t;
#define WIN_THREADS
#define THREADS

#endif

#else

#include <pthread.h>
#define mutex_lock(a) pthread_mutex_lock(a)
#define mutex_unlock(a) pthread_mutex_unlock(a)
#define thread_self() pthread_self()
#define threadid_self() 0  //FIXME t.b.d.
typedef pthread_mutex_t mutexType;
typedef pthread_t threadId_t;
#define THREADS
static pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
static bool cont_condition;

#endif


/* Copied from main.c in ngspice*/
#if defined(__MINGW32__)
#include <stdarg.h>
/* remove type incompatibility with winnt.h*/
#undef BOOLEAN
#include <windef.h>
#include <winbase.h>  /* Sleep */
#elif defined(_MSC_VER)
#include <stdarg.h>
/* remove type incompatibility with winnt.h*/
#undef BOOLEAN
#include <windows.h> /* Sleep */
#include <process.h> /* _getpid */
#define dup _dup
#define dup2 _dup2
#define open _open
#define close _close
#else
#include <unistd.h> /* usleep */
#endif /* __MINGW32__ */

#include "ngspice/iferrmsg.h"
#include "ngspice/ftedefs.h"
#include "ngspice/devdefs.h"
#include "spicelib/devices/dev.h"
#include "spicelib/analysis/analysis.h"
#include "misc/ivars.h"
#include "frontend/resource.h"
#include "frontend/com_measure2.h"
#include "frontend/outitf.h"
#include "ngspice/memory.h"
#include "frontend/com_measure2.h"
#include "frontend/misccoms.h"
#include "ngspice/stringskip.h"
#include "frontend/variable.h"

/* To interupt a spice run */
#include <signal.h>
typedef void (*sighandler)(int);

#include <setjmp.h>
#include "frontend/signal_handler.h"

/*Included for the module to access data*/
#include "ngspice/dvec.h"
#include "ngspice/plot.h"

#ifdef __CYGWIN__
#undef WIN32
#endif
#include  "ngspice/sim.h"

/*For get_output*/
#include <fcntl.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#define S_IRWXU _S_IWRITE
#endif


#ifdef XSPICE
#include "ngspice/evtshared.h"
extern bool wantevtdata;
#endif


/********** includes copied from main.c ************/
#ifdef CIDER
# include "ngspice/numenum.h"
# include "maths/misc/accuracy.h"
#endif

/********** global variables copied from main.c ************/
FILE* slogp = NULL;          /* soa log file ('--soa-log file' command line option) */

/* Frontend and circuit options */
IFsimulator* ft_sim = NULL;

char* errRtn;     /* name of the routine declaring error */
char* errMsg = NULL;     /* descriptive message about what went wrong */
char* cp_program; /* program name 'ngspice' */

char* Infile_Path = NULL; /* Path to netlist input file */

char* hlp_filelist[] = { "ngspice", NULL };


/* Allocate space for global constants declared in const.h
 * and set their values */
double CONSTroot2 = CONSTsqrt2;
double CONSTvt0 = CONSTboltz * REFTEMP / CHARGE;
double CONSTKoverQ = CONSTboltz / CHARGE;
double CONSTe = CONSTnap;

IFfrontEnd* SPfrontEnd = NULL;
int DEVmaxnum = 0;

const bool ft_nutmeg = FALSE;
extern struct comm spcp_coms[];
struct comm* cp_coms = spcp_coms;

/* Main options */
static bool ft_servermode = FALSE;
bool ft_batchmode = FALSE;
bool ft_pipemode = FALSE;
bool rflag = FALSE; /* has rawfile */

/* Frontend options */
bool ft_intrpt = FALSE;     /* Set by the (void) signal handlers. TRUE = we've been interrupted. */
bool ft_setflag = FALSE;    /* TRUE = Don't abort simulation after an interrupt. */
char* ft_rawfile = "rawspice.raw";

#ifdef XSPICE
bool wantevtdata = FALSE;
#endif

bool orflag = FALSE; /* global for -o option */

/* Globals definitions for Machine Accuracy Limits
 * (needed by CIDER)
 */
double BMin;                /* lower limit for B(x) */
double BMax;                /* upper limit for B(x) */
double ExpLim;              /* limit for exponential */
double Accuracy;            /* accuracy of the machine */
double MuLim, MutLim;

IFfrontEnd nutmeginfo = {
    IFnewUid,
    IFdelUid,
    OUTstopnow,
    seconds,
    OUTerror,
    OUTerrorf,
    OUTpBeginPlot,
    OUTpData,
    OUTwBeginPlot,
    OUTwReference,
    OUTwData,
    OUTwEnd,
    OUTendPlot,
    OUTbeginDomain,
    OUTendDomain,
    OUTattributes
};


#ifdef CIDER
/* Global debug flags from CIDER, soon they will become
 * spice variables :)
 */
int ONEacDebug = FALSE;
int ONEdcDebug = TRUE;
int ONEtranDebug = TRUE;
int ONEjacDebug = FALSE;

int TWOacDebug = FALSE;
int TWOdcDebug = TRUE;
int TWOtranDebug = TRUE;
int TWOjacDebug = FALSE;

/* CIDER Global Variable Declarations */

int BandGapNarrowing;
int TempDepMobility, ConcDepMobility, FieldDepMobility, TransDepMobility;
int SurfaceMobility, MatchingMobility, MobDeriv;
int CCScattering;
int Srh, Auger, ConcDepLifetime, AvalancheGen;
int FreezeOut = FALSE;
int OneCarrier;

int MaxIterations = 100;
int AcAnalysisMethod = DIRECT;

double Temp, RelTemp, Vt;
double RefPsi;/* potential at Infinity */
double EpsNorm, VNorm, NNorm, LNorm, TNorm, JNorm, GNorm, ENorm;

/* end cider globals */
#endif /* CIDER */

struct variable* (*if_getparam)(CKTcircuit* ckt, char** name, char* param, int ind, int do_model);
/***********************************************************/

extern IFsimulator SIMinfo;

extern struct comm spcp_coms[ ];
extern void DevInit(void);
extern wordlist *cp_varwl(struct variable *var);
extern void create_circbyline(char *line, bool reset, bool lastline);

static int SIMinit(IFfrontEnd *frontEnd, IFsimulator **simulator);

void exec_controls(wordlist *shcontrols);
void rem_controls(void);

/*The current run (to get variable names, etc)*/
static runDesc *cur_run;

double getisrcval(double time, char *iname);

int sh_fputsll(const char *input, FILE* outf);

int sh_ExecutePerLoop(void);
double getvsrcval(double, char*);
int sh_vecinit(runDesc *run);

ATTRIBUTE_NORETURN void shared_exit(int status);

void sighandler_sharedspice(int num);

void wl_delete_first(wordlist **wlstart, wordlist **wlend);

int add_bkpt(void);
int sharedsync(double*, double*, double, double, double, int, int*, int);

void sh_delete_myvec(void);

#ifdef XSPICE
void shared_send_event(int, double, double, char *, void *, int, int);
void shared_send_dict(int, int, char*, char*);
#endif

#if !defined(low_latency)
static char* outstorage(char*, bool);
static void printsend(void);
#endif

extern wordlist* sourceinfo;

static int totalreset(void);
extern void rem_controls(void);
extern void destroy_wallace(void);
extern void sh_delete_myvec(void);
extern IFsimulator SIMinfo;
extern void spice_destroy_devices(void); /* FIXME need a better place */

extern void destroy_const_plot(void);
extern void com_destroy(wordlist* wl);
extern void com_unalias(wordlist* wl);
extern void com_undefine(wordlist* wl);
extern void com_remcirc(wordlist* wl);
extern void unset_all(void);

#include "ngspice/sharedspice.h"

static SendChar* pfcn;
static void* userptr;
static SendStat* statfcn;
static ControlledExit* ngexit;
static SendData* datfcn;
static SendInitData* datinitfcn;
static BGThreadRunning* bgtr;
static GetVSRCData* getvdat;
static GetISRCData* getidat;
static GetSyncData* getsync;
static pvector_info myvec = NULL;
#ifdef XSPICE
static struct dvec *infovec = NULL;
#endif
char **allvecs = NULL;
char **allplots = NULL;
static bool noprintfwanted = FALSE;
static bool nostatuswanted = FALSE;
static bool nodatawanted = FALSE;
static bool nodatainitwanted = FALSE;
static bool nobgtrwanted = FALSE;
static bool wantvdat = FALSE;
static bool wantidat = FALSE;
static bool wantsync = FALSE;
static NG_BOOL immediate = FALSE;
static NG_BOOL coquit = FALSE;
static jmp_buf errbufm, errbufc;
static int intermj = 1;
#ifdef XSPICE
static SendInitEvtData* sendinitevt;
static SendEvtData* sendevt;
#endif
static void* euserptr;
static wordlist *shcontrols;

// thread IDs
unsigned int main_id, ng_id, command_id;

#ifdef THREADS
mutexType triggerMutex;
mutexType allocMutex;
mutexType fputsMutex;
mutexType vecreallocMutex;
#endif

/* initialization status */
static bool is_initialized = FALSE;
static char* no_init = "Error: ngspice is not initialized!\n   Run ngSpice_Init first";

/* identifier for this ngspice invocation */
int ng_ident = 0;


static struct plot *
get_plot_byname(char* plotname)
{
    struct plot *pl;
    pl = plot_list;
    while (pl) {
        if(cieq(pl->pl_typename, plotname))
            break;
        pl = pl->pl_next;
    }
    return pl;
}

/* -------------------------------------------------------------------------- */
static int
SIMinit(IFfrontEnd* frontEnd, IFsimulator** simulator)
{
    spice_init_devices();
    SIMinfo.numDevices = DEVmaxnum = num_devices();
    SIMinfo.devices = devices_ptr();
    SIMinfo.numAnalyses = spice_num_analysis();

    /* va: we recast, because we use only the public part */
    SIMinfo.analyses = (IFanalysis**)spice_analysis_ptr();


#ifdef CIDER
    /* Evaluates limits of machine accuracy for CIDER */
    evalAccLimits();
#endif /* CIDER */

    SPfrontEnd = frontEnd;
    *simulator = &SIMinfo;

    return OK;
} /* end of function SIMinit */

/******************************************************************/
/*     Main spice command executions and thread control           */
/*****************************************************************/

#ifdef THREADS

static threadId_t tid, printtid, tid2;

static bool fl_running = FALSE;
static NG_BOOL fl_exited = TRUE;

static bool printstopp = FALSE;
static bool ps_exited = TRUE;

#if defined(__MINGW32__) || defined(_MSC_VER)
#define EXPORT_FLAVOR WINAPI
#else
#define EXPORT_FLAVOR
#endif

/* starts a thread to run the controls, started when bg thread finishes */
static void * EXPORT_FLAVOR
_cthread_run(void *controls)
{
    wordlist *wl;
#ifdef HAVE_LIBPTHREAD
    if (!cont_condition)
        printf("Prepared to start controls after bg_run has finished\n");
    pthread_mutex_lock(&triggerMutex);
    cont_condition = FALSE;
    do {
        pthread_cond_wait(&cond, &triggerMutex);
    } while (!cont_condition);
    pthread_mutex_unlock(&triggerMutex);
#endif
    fl_exited = FALSE;
    for (wl = controls; wl; wl = wl->wl_next)
        cp_evloop(wl->wl_word);
    fl_exited = TRUE;
#ifdef HAVE_LIBPTHREAD
    cont_condition = FALSE;
#endif
    wl_free(controls);
    return NULL;
}

/* starts a background thread, e.g. from command bg_run,
   releases controls thread tid2  */
static void * EXPORT_FLAVOR
_thread_run(void *string)
{
    ng_id = threadid_self();
    fl_exited = FALSE;
    /* notify caller that thread is running */
    if (!nobgtrwanted)
        bgtr(fl_exited, ng_ident, userptr);
//  bgtid = thread_self();
    cp_evloop((char *)string);
    FREE(string);
// #ifdef __MINGW32__
//     bgtid.p = NULL;
//     bgtid.x = 0;
// #else
//     bgtid = (threadId_t)0;
// #endif
    fl_exited = TRUE;
    /* notify caller that thread has exited */
    if (!nobgtrwanted)
        bgtr(fl_exited, ng_ident, userptr);
    /* release thread tid2 */
    if (tid2) {
#ifdef HAVE_LIBPTHREAD
        pthread_mutex_lock(&triggerMutex);
        cont_condition = TRUE;
        pthread_cond_signal(&cond);
        pthread_mutex_unlock(&triggerMutex);
        pthread_join(tid2, NULL);
#elif defined _MSC_VER || defined __MINGW32__
        ResumeThread(tid2);
#else

#endif
        tid2 = 0;
    }
    return NULL;
}


/* Stops a running background thread, hopefully */
static int EXPORT_FLAVOR
_thread_stop(void)
{
    int timeout = 0;

    if (fl_running) {
        while (!fl_exited && timeout < 100) {
            /* ft_intrpt is the flag to stop simulation, if set TRUE !
               E.g. SPfrontEnd->IFpauseTest() in dctran.c points to
               OUTstopnow(void), which returns 1, which leads dctran
               to return with -1 and thus terminates the simulation*/
            ft_intrpt = TRUE;
            timeout++;
#if defined(__MINGW32__) || defined(_MSC_VER)
            Sleep(100); // va: windows native
#else
            usleep(10000);
#endif
        }
        if (!fl_exited) {
            fprintf(stderr, "Error: Couldn't stop ngspice\n");
            return EXIT_BAD;
        }
        else
            fprintf(stdout, "Background thread stopped with timeout = %d\n", timeout);

        fl_running = FALSE;
        ft_intrpt = FALSE;
        return EXIT_NORMAL;
    } else {
        fprintf(stderr, "Spice not running\n");
    }
    return EXIT_NORMAL;
}


void
sighandler_sharedspice(int num)
{
    NG_IGNORE(num);
    if (fl_running)
        _thread_stop();
    return;
}

#endif /*THREADS*/

/* create a suspended thread tid2 that is activated when bg_run has finished.
   It executes the .control commands. If the argument is NULL, the thread is
   started with the existing controls (e.g. during command 'reset'. */
void
exec_controls(wordlist *newcontrols)
{
    if (newcontrols && newcontrols->wl_word && !eq(newcontrols->wl_word,"")) {
        shcontrols = newcontrols;
    }
    else {
        tid2 = 0;
        return;
    }

#ifdef THREADS
#ifdef HAVE_LIBPTHREAD
    cont_condition = FALSE;
    usleep(20000); /* wait a little */
    pthread_create(&tid2, NULL, (void * (*)(void *))_cthread_run, (void *)shcontrols);
#elif defined _MSC_VER || defined __MINGW32__
    tid2 = (HANDLE)_beginthreadex(NULL, 0, (unsigned int(__stdcall *)(void *))_cthread_run,
        (void*)shcontrols, CREATE_SUSPENDED, NULL);
#else
    tid2 = CreateThread(NULL, 0, (PTHREAD_START_ROUTINE)_cthread_run, (void*)shcontrols,
        0, NULL);
#endif
#else
    wordlist *wl;
    for (wl = shcontrols; wl; wl = wl->wl_next)
        cp_evloop(wl->wl_word);
    wl_free(shcontrols);
#endif
}

/* free controls after 'quit' */
void rem_controls(void)
{
    wl_free(shcontrols);
}


/* run a ngspice command */
static int
runc(char* command)
{
    char buf[1024] = "";
#ifdef THREADS
#ifndef low_latency
    int timeout = 0;
#endif
    char *string;
    bool fl_bg = FALSE;
    command_id = threadid_self();
    /* run task in background if command is preceeded by "bg_" */
    if (!cieq("bg_halt", command) && !cieq("bg_pstop", command)
        && !cieq("bg_ctrl", command) && ciprefix("bg_", command)) {
        strncpy(buf, command+3, 1024);
        fl_bg = TRUE;
    }
#ifndef low_latency
    /* stop the printf thread 'printsend()' */
    else if (cieq("bg_pstop", command)) {
        while (!ps_exited && timeout < 100) {
            printstopp = TRUE;

#if defined __MINGW32__ || defined _MSC_VER
            Sleep(100); // va: windows native
#else
            usleep(10000);
#endif
            timeout++;
        }
        if (!ps_exited) {
            fprintf(stderr, "Error: Couldn't stop printsend thread\n");
            return EXIT_BAD;
        }
        else
            fprintf(stdout, "Printsend thread stopped with timeout = %d\n", timeout);

        printstopp = FALSE;
        return 2;
    }
#endif
    else
        strncpy(buf, command, 1024);
#else
    strncpy(buf, command, 1024);
#endif

#ifdef THREADS
    /* run in the background */
    if (fl_bg && fl_exited) {
        if (fl_running)
            _thread_stop();
        fl_running = TRUE;
        string = copy(buf);     /*as buf gets freed fairly quickly*/
#ifdef HAVE_LIBPTHREAD
        pthread_create(&tid, NULL, (void * (*)(void *))_thread_run, (void *)string);
        pthread_detach(tid);
#elif defined _MSC_VER || defined __MINGW32__
        tid = (HANDLE)_beginthreadex(NULL, 0, (unsigned int (__stdcall *)(void *))_thread_run,
            (void*)string, 0, NULL);
#else
        tid = CreateThread(NULL, 0, (PTHREAD_START_ROUTINE)_thread_run, (void*)string,
                         0, NULL);
#endif
    } else
        /* bg_halt (pause) a bg run */
        if (!strcmp(buf, "bg_halt")) {
            return _thread_stop();
        /* bg_ctrl prepare running the controls after bg_run */
        } else if (!strcmp(buf, "bg_ctrl")) {
            if (shcontrols)
                exec_controls(wl_copy(shcontrols));
            else
                fprintf(stderr, "Warning: No .control commands available, bg_ctrl skipped\n");
            return 0;
        } else
            /* cannot do anything if ngspice is running in the bg*/
            if (fl_running) {
                if (fl_exited) {
                    _thread_stop();
                    cp_evloop(buf);
                } else {
                    fprintf(stderr, "Warning: cannot execute \"%s\", type \"bg_halt\" first\n", buf);
                }
            } else {
                /*do the command*/
                cp_evloop(buf);
            }
#else
    cp_evloop(buf);
#endif /*THREADS*/
    return 0;
}

/* -------------------------------------------------------------------------- */
/* Read an initialisation file.
dir    is the directory (use NULL or "" for current directory)
name   is the initialisation file's name
Return true on success
SJB 25th April 2005 */
static bool
read_initialisation_file(const char *dir, const char *name)
{
    const char *path;
    bool result = FALSE;

    /* check name */
    if (!name || *name == '\0')
        return FALSE;   /* Fail; name needed */

                        /* contruct the full path */
    if (!dir || *dir == '\0') {
        path = name;
    }
    else {
        path = tprintf("%s" DIR_PATHSEP "%s", dir, name);
        if (!path)
            return FALSE;    /* memory allocation error */
    }

    /* now access the file */
#ifdef HAVE_UNISTD_H
    if (access(path, R_OK) == 0)
        result = TRUE;
#else
    {
        FILE *fp = fopen(path, "r");
        if (fp) {
            fclose(fp);
            result = TRUE;
        }
    }
#endif

    if (result) {
        inp_source(path);
#ifdef TRACE
        printf("Init file: '%s'\n", path);
#endif
    }

    if (path != name)
        tfree(path);

    return result;
}

/* -------------------------------------------------------------------------- */

/**********************************************************/
/* The functions exported explicitely from shared ngspice */
/**********************************************************/



#ifdef THREADS

/* Checks if ngspice is running in the background */
IMPEXP
NG_BOOL
ngSpice_running (void)
{
    return (fl_running && !fl_exited);
}
#endif

/* Set variable no_spinit, if reading 'spinit' is not wanted. */
IMPEXP
int
ngSpice_nospinit(void)
{
    bool t = TRUE;
    cp_vset("no_spinit", CP_BOOL, &t);
    return 0;
}

/* Set variable no_spiceinit, if reading '.spiceinit' is not wanted. */
IMPEXP
int
ngSpice_nospiceinit(void)
{
    bool t = TRUE;
    cp_vset("no_spicenit", CP_BOOL, &t);
    return 0;
}

/* Initialise external voltage source and synchronization */
IMPEXP
int
ngSpice_Init_Sync(GetVSRCData *vsrcdat, GetISRCData *isrcdat, GetSyncData *syncdat, int *ident, void *userData)
{
    getvdat = vsrcdat;
    getidat = isrcdat;
    getsync = syncdat;
    /* set userdata, but don't overwrite with NULL */
    if (userData)
        userptr = userData;
    /* set ngspice shared lib identification number */
    if (ident)
        ng_ident = *ident;
    /* if caller sends NULL, don't try to retrieve voltage */
    if (getvdat) {
        wantvdat = TRUE;
    }
    /* if caller sends NULL, don't try to retrieve current */
    if (getidat) {
        wantidat = TRUE;
    }
    /* if caller sends NULL, don't synchronize */
    if (getsync) {
        wantsync = TRUE;
    }
    return 0;
}


/* Initialise ngspice and setup native methods */
IMPEXP
int
ngSpice_Init(SendChar* printfcn, SendStat* statusfcn, ControlledExit* ngspiceexit,
             SendData* sdata, SendInitData* sinitdata, BGThreadRunning* bgtrun, void* userData)
{
    sighandler old_sigsegv = NULL;

    struct variable* sourcepathvar;

    pfcn = printfcn;
    /* if caller sends NULL, don't send printf strings */
    if (!pfcn)
        noprintfwanted = TRUE;
    userptr = userData;
    statfcn = statusfcn;
    /* if caller sends NULL, don't send status data */
    if (!statfcn)
        nostatuswanted = TRUE;
    ngexit = ngspiceexit;
    datfcn = sdata;
    /* if caller sends NULL, don't send data */
    if (!datfcn)
        nodatawanted = TRUE;
    /* if caller sends NULL, don't initialize and send data */
    datinitfcn = sinitdata;
    if (!datinitfcn)
        nodatawanted = nodatainitwanted = TRUE;
    bgtr = bgtrun;
    if (!bgtr)
        nobgtrwanted = TRUE;
    immediate = FALSE;

    cp_nocc = TRUE;

#ifdef THREADS
    /* init the mutexes */
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_init(&triggerMutex, NULL);
    pthread_mutex_init(&allocMutex, NULL);
    pthread_mutex_init(&fputsMutex, NULL);
    pthread_mutex_init(&vecreallocMutex, NULL);
    cont_condition = FALSE;
#else
#ifdef SRW
    InitializeSRWLock(&triggerMutex);
    InitializeSRWLock(&allocMutex);
    InitializeSRWLock(&fputsMutex);
    InitializeSRWLock(&vecreallocMutex);
#else
    InitializeCriticalSection(&triggerMutex);
    InitializeCriticalSection(&allocMutex);
    InitializeCriticalSection(&fputsMutex);
    InitializeCriticalSection(&vecreallocMutex);
#endif
#endif
    // Id of primary thread
    main_id =  threadid_self();
#endif

    if (!cp_getvar("nosighandling", CP_BOOL, NULL, 0))
        old_sigsegv = signal(SIGSEGV, (SIGNAL_FUNCTION) sigsegvsh);

    ft_rawfile = NULL;
    ivars(NULL);

    cp_in = stdin;
    cp_out = stdout;
    cp_err = stderr;

    /*timer*/
    init_time();

    /*IFsimulator struct initilised*/
    SIMinit(&nutmeginfo, &ft_sim);

    /* program name*/
    cp_program = ft_sim->simulator;

    /* initialze random number generator with seed = 1 */
    int ii = 1;
    cp_vset("rndseed", CP_NUM, &ii);
    com_sseed(NULL);

    /* set a boolean variable to be used in .control sections */
    bool sm = TRUE;
    cp_vset("sharedmode", CP_BOOL, &sm);

    /* set a boolean variable when XSPICE and/or OSDI is enabled,
       to be used in spinit etc. */
#if defined(XSPICE)
    cp_vset("xspice_enabled", CP_BOOL, &sm);
#endif
#if defined(OSDI)
    cp_vset("osdi_enabled", CP_BOOL, &sm);
#endif

    /*parameter fetcher, used in show, alter, altermod */
    if_getparam = spif_getparam_special;

    /* Get startup system limits */
    init_rlimits();

    /*Command prompt stuff */
    ft_cpinit();

    /* Read the user config files */
#ifdef HAVE_PWD_H
    /* Try to source either .spiceinit or ~/.spiceinit. */
    if (access(".spiceinit", 0) == 0) {
        inp_source(".spiceinit");
    } else {
        char *s;
        struct passwd *pw;
        pw = getpwuid(getuid());

        s = tprintf("%s" DIR_PATHSEP "%s", pw->pw_dir, INITSTR);

        if (access(s, 0) == 0)
            inp_source(s);

        tfree(s);
    }
#else /* ~ HAVE_PWD_H */
    /* load user's initialisation file
       try accessing the initialisation file .spiceinit in a user provided
       path read from environmental variable SPICE_USERINIT_DIR,
       if that fails try the alternate name spice.rc, then look into
       the current directory, then the HOME directory, then into USERPROFILE.
       Don't read .spiceinit, if ngSpice_nospiceinit() has been called. */
    if (!cp_getvar("no_spiceinit", CP_BOOL, NULL, 0)) {
        do {
            {
                const char* const userinit = getenv("SPICE_USERINIT_DIR");
                if (userinit) {
                    if (read_initialisation_file(userinit, INITSTR) != FALSE) {
                        break;
                    }
                    if (read_initialisation_file(userinit, ALT_INITSTR) != FALSE) {
                        break;
                    }
                }
            }

            if (read_initialisation_file("", INITSTR) != FALSE) {
                break;
            }
            if (read_initialisation_file("", ALT_INITSTR) != FALSE) {
                break;
            }

            {
                const char* const home = getenv("HOME");
                if (home) {
                    if (read_initialisation_file(home, INITSTR) != FALSE) {
                        break;
                    }
                    if (read_initialisation_file(home, ALT_INITSTR) != FALSE) {
                        break;
                    }
                }
            }

            {
                const char* const usr = getenv("USERPROFILE");
                if (usr) {
                    if (read_initialisation_file(usr, INITSTR) != FALSE) {
                        break;
                    }
                    if (read_initialisation_file(usr, ALT_INITSTR) != FALSE) {
                        break;
                    }
                }
            }
        } while (0); /* end of case that init file is read */
    }
    else {
        fprintf(stdout, "Note: .spiceinit is ignored, because ngSpice_nospiceinit() has been called.\n");
    }

#endif /* ~ HAVE_PWD_H */

    if (!cp_getvar("nosighandling", CP_BOOL, NULL, 0))
        signal(SIGSEGV, old_sigsegv);

    /* initialize display to 'no display at all'*/
    DevInit();

#ifdef FastRand
// initialization and seed for FastNorm Gaussian random generator
    {
        unsigned int rseed = 66;
        initnorm (0, 0);
        if (!cp_getvar("rndseed", CP_NUM, &rseed, 0)) {
            time_t acttime = time(NULL);
            rseed = (unsigned int) acttime;
        }
        initnorm (rseed, 2);
        fprintf (cp_out, "SoS %f, seed value: %ld\n", renormalize(), rseed);
    }
#elif defined (WaGauss)
        initw();
#endif

    fprintf(cp_out,
            "******\n"
            "** %s-%s shared library\n",
            ft_sim->simulator, ft_sim->version);
    if (*Spice_Build_Date != 0)
        fprintf(cp_out, "** Creation Date: %s\n", Spice_Build_Date);
    fprintf(cp_out, "******\n");

    is_initialized = TRUE;

    if(!myvec)
        myvec = TMALLOC(vector_info, sizeof(vector_info));

    /* Read first entry of sourcepath var, set Infile_path for code models */
    if ( cp_getvar("sourcepath", CP_LIST, &sourcepathvar, 0)) {
        Infile_Path = copy(sourcepathvar->va_string);
    }

#if !defined(low_latency)
    /* If caller has sent valid address for pfcn */
    if (!noprintfwanted)

#ifdef HAVE_LIBPTHREAD
        pthread_create(&printtid, NULL, (void * (*)(void *))printsend, (void *)NULL);
#elif defined _MSC_VER || defined __MINGW32__
        printtid = (HANDLE)_beginthreadex(NULL, 0, (unsigned int (__stdcall *)(void *))printsend,
            (void*) NULL, 0, NULL);
#else
        printtid = CreateThread(NULL, 0, (PTHREAD_START_ROUTINE) printsend, NULL,
                         0, NULL);
#endif

#endif

    return 0;
}


/* to be called upon 'quit' */
void
sh_delete_myvec(void)
{
    tfree(myvec);
#ifdef XSPICE
    if (infovec) {
        dvec_free(infovec->v_scale);
        dvec_free(infovec);
    }
#endif
}

/* retrieve a ngspice command from caller and run it
   immediately.
   If NULL is sent, we clear the command memory */
IMPEXP
int  ngSpice_Command(char* comexec)
{
    if (!is_initialized) {
        return 1;
    }

    /* delete existing command memory */
    if (comexec == NULL) {
        cp_resetcontrol(FALSE);
        return 0;
    }
    /* Check if command is reasonable */
    if (*comexec == '\0') {
        fprintf(stderr, "Warning: Received empty string as command, ignored");
        return 1;
    }

    if (ft_ngdebug)
        fprintf(stdout, "\nngSpiceCommand: received command '%s'\n", comexec);

    if ( ! setjmp(errbufc) ) {

        immediate = FALSE;
        intermj = 1;

        if (!is_initialized) {
           fprintf(stderr, no_init);
           return 1;
       }
       runc(comexec);
       /* main thread prepares immediate detaching of dll */
       immediate = TRUE;
       return 0;
    }
    return 1;
}

/* Set the input path for files loaded by code models
   like d_state, file_source, d_source.
   Useful when netlist is sent by ngSpice_Circ and therefore
   Infile_Path cannot be retrieved automatically.
   If NULL is sent, return the current Infile_Path. */
IMPEXP
char *ngCM_Input_Path(const char* path)
{
    /* override existing path */
    if (path) {
        txfree(Infile_Path);
        Infile_Path = copy(path);
    }
    fprintf(stdout, "Note: Codel model file loading path is %s\n", Infile_Path);
    return Infile_Path;
}

/* Return information about a vector to the caller */
IMPEXP
pvector_info  ngGet_Vec_Info(char* vecname)
{
    struct dvec* newvec;

    if (!is_initialized) {
        fprintf(stderr, no_init);
        return NULL;
    }

#ifdef XSPICE
    /* If vector is derived from event data, free it */
    if (infovec) {
        dvec_free(infovec->v_scale);
        dvec_free(infovec);
        infovec = NULL;
    }
#endif

    newvec = vec_get(vecname);

    if (newvec == NULL) {
        fprintf(stderr, "Error: vector %s not found!\n", vecname);
        return NULL;
    }
    if (newvec->v_numdims > 1) {
        fprintf(stderr, "Error: vector %s is multidimensional!\n  This is not yet handled\n!", vecname);
        return NULL;
    }

    myvec->v_name = newvec->v_name;
    myvec->v_type = newvec->v_type;
    myvec->v_flags = newvec->v_flags;
    myvec->v_realdata = newvec->v_realdata;
    myvec->v_compdata = newvec->v_compdata;
    myvec->v_length = newvec->v_length;

#ifdef XSPICE
    /* If we have a vector derived from event data, store its pointer */
    if (newvec->v_scale && newvec->v_scale->v_name && eq(newvec->v_scale->v_name, "step"))
        infovec = newvec;
#endif

    return myvec;
}

/* Receive a circuit from the caller as a
   pointer to an array of char* .
   Last entry in array has to be NULL
*/
IMPEXP
int ngSpice_Circ(char** circa){
    int entries = 0, i;
    char* newline;
    bool reset = FALSE, lastline = FALSE;

    if ( ! setjmp(errbufm) ) {
        intermj = 0;
        immediate = FALSE;
        /* count the entries */
        while (circa[entries]) {
            char* line = skip_ws(circa[entries++]);
            if (ciprefix(".end", line) && (line[4] == '\0' || isspace_c(line[4])))
                break;
        }

        if (ft_ngdebug)
            fprintf(stdout, "\nngspiceCirc: received netlist array with %d entries\n", entries);
        /* create a local copy (to be freed in inpcom.c) */
        for (i = 0; i < entries; i++) {
            newline = copy(circa[i]);
            if (i == 0)
                reset = TRUE;
            else
                reset = FALSE;
            if (i == entries - 1)
                lastline = TRUE;
            create_circbyline(newline, reset, lastline);
        }
        return 0;
    }
    /* upon error */
    return 1;
}


/* return to the caller a pointer to the name of the current plot */
IMPEXP
char* ngSpice_CurPlot(void)
{
    struct plot *pl = plot_cur;
    return pl->pl_typename;
}

/* return to the caller a pointer to an array of all plots created
by ngspice. Last entry in the array is NULL. */
IMPEXP
char** ngSpice_AllPlots(void)
{
    int len = 0, i = 0;
    struct plot *pl = plot_list;
    if (allplots)
        tfree(allplots);

    while (pl) {
        len++;
        pl = pl->pl_next;
    }
    allplots = TMALLOC(char*, len+1);
    pl = plot_list;
    for (i = 0; i < len; i++) {
        allplots[i] = pl->pl_typename;
        pl = pl->pl_next;
    }
    allplots[len] = NULL;
    return allplots;
}

/* return to the caller a pointer to an array of vector names in the plot
named by plotname. Last entry in the array is NULL. */
IMPEXP
char** ngSpice_AllVecs(char* plotname)
{
    struct dvec *d;
    int len = 0, i = 0;
    struct plot *pl;

    if (allvecs)
        tfree(allvecs);

    /* get the plot plotname */
    pl = get_plot_byname(plotname);

    if (pl)
        for (d = pl->pl_dvecs; d; d = d->v_next)
            len++;

    if (len == 0) {
        fprintf(cp_err, "Error: There are no vectors currently active.\n");
        return NULL;
    }

    allvecs = TMALLOC(char*, len + 1);

    for (d = pl->pl_dvecs, i = 0; d; d = d->v_next, i++)
        allvecs[i] = d->v_name;

    allvecs[len] = NULL;

    return allvecs;
}


static double *bkpttmp = NULL;
static int bkpttmpsize = 0;

/* set a breakpoint in ngspice */
IMPEXP
NG_BOOL ngSpice_SetBkpt(double time)
{
    int error;
    CKTcircuit *ckt = NULL;

    if (!ft_curckt || !ft_curckt->ci_ckt) {
        fprintf(cp_err, "Error: no circuit loaded.\n");
        return(FALSE);
    }

    ckt = ft_curckt->ci_ckt;
    if (ckt->CKTbreakSize == 0) {
    /* breakpoints have not yet been set up, so store here preliminary
    and add with fcn add_bkpt() called from DCTran() */
        if (bkpttmp == NULL) {
            bkpttmp = TMALLOC(double, bkpttmpsize + 1);
            if(bkpttmp == NULL)
                return(FALSE);
            bkpttmpsize++;
        }
        else {
            bkpttmp = TREALLOC(double, bkpttmp, bkpttmpsize + 1);
            bkpttmpsize++;
        }
        bkpttmp[bkpttmpsize-1] = time;
        error = 0;
    }
    else
        error = CKTsetBreak(ckt, time);
    if(error)
        return(FALSE);
    return(TRUE);
}

#ifdef XSPICE
/* return callback initialization addresses to caller */
IMPEXP
int  ngSpice_Init_Evt(SendEvtData* sevtdata, SendInitEvtData* sinitevtdata, void* userData)
{
    if (sevtdata)
        wantevtdata = TRUE;
    else
        wantevtdata = FALSE;
    sendinitevt = sinitevtdata;
    sendevt = sevtdata;
    euserptr = userData;
    return(TRUE);
}

/* Get info about the event node vector.
If node_name is NULL, just delete previous data */
IMPEXP
pevt_shared_data ngGet_Evt_NodeInfo(char* node_name)
{
    return EVTshareddata(node_name);
}

/* get a list of all event nodes */
IMPEXP
char** ngSpice_AllEvtNodes(void)
{
    return EVTallnodes();
}
#endif

/* Lock/unlock realloc of result vectors during plotting */
IMPEXP
int ngSpice_LockRealloc(void)
{
    mutex_lock(&vecreallocMutex);
    return 1;
}

IMPEXP
int ngSpice_UnlockRealloc(void)
{
    mutex_unlock(&vecreallocMutex);
    return 1;
}

/* Reset ngspice as far as possible */
IMPEXP
int ngSpice_Reset(void)
{
    if (!is_initialized)
        return 1;
    fprintf(stdout, "Note: Resetting ngspice\n\n");
    return totalreset();
}

/* add the preliminary breakpoints to the list.
   called from dctran.c */
int
add_bkpt(void)
{
    int i;
    int error = 0;
    CKTcircuit *ckt =  ft_curckt->ci_ckt;

    if(bkpttmp && (bkpttmpsize > 0)) {
        for (i = 0; i < bkpttmpsize; i++)
            error = CKTsetBreak(ckt, bkpttmp[i]);
        FREE(bkpttmp);
        bkpttmpsize = 0;
    }
    if(error)
        return(error);
    return(OK);
}


/* use the original vprintf() in the rest of this file
 *   instead of the redirected variant
 */
#undef vfprintf


/*------------------------------------------------------*/
/* Redefine the vfprintf() functions for callback       */
/*------------------------------------------------------*/

/* handling of escape characters (extra \ added) only, if
   'set addescape' is given in .spiceinit */

int
sh_vfprintf(FILE *f, const char *fmt, va_list args)
{
    char buf[1024];
    char *p/*, *s*/;
    int nchars, /*escapes,*/ result;
    size_t size;


    if ((fileno(f) != STDOUT_FILENO && fileno(f) != STDERR_FILENO &&
         f != stderr && f != stdout)
// #ifdef THREADS
//        || (fl_running && bgtid == thread_self())
// #endif
        )
        return vfprintf(f, fmt, args);

    p = buf;

    // size: how much ist left for chars and terminating '\0'
    size = sizeof(buf);
    // assert(size > 0);

    for (;;) {
        va_list ap;

        va_copy(ap, args);
        nchars = vsnprintf(p, size, fmt, ap);
        va_end(ap);

        if(nchars == -1) {           // compatibility to old implementations
            size *= 2;
        } else if (size < (size_t)nchars + 1) {
            size = (size_t)nchars + 1;
        } else {
            break;
        }

        if(p == buf)
            p = TMALLOC(char, size);
        else
            p = TREALLOC(char, p, size);
    }

    /* add / to escape characters, if 'set addescape' is called in .spiceinit */
    if (cp_getvar("addescape", CP_BOOL, NULL, 0)) {
        size_t escapes;
        const char * const escape_chars = "$[]\"\\";
        char *s = p;
        for (escapes = 0; ; escapes++) {
            s = strpbrk(s, escape_chars);
            if (!s)
                break;
            s++;
        }

        if (escapes) {

            size_t new_size = (size_t)nchars + escapes + 1;
            char *src, *dst;

            if (p != buf) {
                p = TREALLOC(char, p, new_size);
            } else if (new_size > sizeof(buf)) {
                p = TMALLOC(char, new_size);
                strcpy(p, buf);
            }

            src = p + nchars;
            dst = src + escapes;

            while (dst > src) {
                char c = *--src;
                *--dst = c;
                if (strchr(escape_chars, c))
                    *--dst = '\\';
            }
        }
    }

    /* use sharedspice.c implementation of fputs (sh_fputs)
       to assess callback function derived from address printfcn received via
       Spice_Init() from caller of ngspice.dll */


    result = sh_fputs(p, f);

    if (p != buf)
        tfree(p);

    return nchars;
}


/*----------------------------------------------------------------------
   Reimplement fprintf() as a call to callback function pfcn
   via sh_vfprintf, sh_fputs, and sh_fputsll
  ----------------------------------------------------------------------*/

int
sh_fprintf(FILE *f, const char *format, ...)
{
    va_list args;
    int rtn;

    va_start (args, format);
    rtn = sh_vfprintf(f, format, args);
    va_end(args);

    return rtn;
}


/*----------------------------------------------------------------------
   Reimplement printf() as a call to callback function pfcn
   via sh_vfprintf, sh_fputs, and sh_fputsll
  ----------------------------------------------------------------------*/

int
sh_printf(const char *format, ...)
{
    va_list args;
    int rtn;

    va_start (args, format);
    rtn = sh_vfprintf(stdout, format, args);
    va_end(args);

    return rtn;
}

int
sh_putc(int inp, FILE* f)
{
    char inpconv[2];
    if ((fileno(f) != STDOUT_FILENO && fileno(f) != STDERR_FILENO &&
    f != stderr && f != stdout))
        return myfputc(inp, f);

    sprintf(inpconv, "%c", inp);
    fputs(inpconv, f);
    return inp;
}

int
sh_fputc(int inp, FILE* f)
{
    char inpconv[2];
    if ((fileno(f) != STDOUT_FILENO && fileno(f) != STDERR_FILENO &&
    f != stderr && f != stdout))
        return myfputc(inp, f);

    sprintf(inpconv, "%c", inp);
    fputs(inpconv, f);
    return inp;
}

/*----------------------------------------------------------------------*/
/* Reimplement fputs() as a call to callback function pfcn              */
/*----------------------------------------------------------------------*/


/* Collect and cat strings.
   If \n is detected, send string
   to caller via pfcn() */
static char* outstringerr = NULL;
static char* outstringout = NULL;

#if defined low_latency || !defined THREADS

/* The strings issued by printf etc. are sent directly to the caller.
   The callback has to be fast enough (low latency). */
int
sh_fputsll(const char *input, FILE* outf)
{
    int result = 0;
    size_t len;
    char *delstring, *newstring, *prstring;
    size_t inputlen = strlen(input);

    /* If caller has sent NULL address for pfcn */
    if (noprintfwanted)
        return -1;

    if (outf == stderr) {
        if (!outstringerr)
            delstring = outstringerr = copy(input);
        else {
            len = strlen(outstringerr);
            delstring = outstringerr = TREALLOC(char, outstringerr, len + inputlen + 2);
            strcat(outstringerr, input);
        }
        if (strchr(input, '\n')) {
            while (outstringerr) {
                newstring = gettok_char(&outstringerr, '\n', FALSE, FALSE);
                if(!newstring)
                    break;
                prstring = tprintf("stderr %s", newstring);

                result = pfcn(prstring, ng_ident, userptr);
                tfree(newstring);
                tfree(prstring);
            }
            /* copy the rest following \n, but without trailing \n to new address */
            if (outstringerr && *outstringerr != '\0')
                outstringerr = copy(outstringerr);
            else
                outstringerr = NULL;
            tfree(delstring);
            return result;
        }
        else if (strchr(input, '\r')) {
            result = pfcn(outstringerr, ng_ident, userptr);
            tfree(outstringerr);
            return result;
        }
    }
    else if (outf == stdout) {
        if (!outstringout)
            delstring = outstringout = copy(input);
        else {
            len = strlen(outstringout);
            delstring = outstringout = TREALLOC(char, outstringout, len + inputlen + 1);
            strcat(outstringout, input);
        }
        if (strchr(input, '\n')) {
            while (outstringout) {
                newstring = gettok_char(&outstringout, '\n', FALSE, FALSE);
                if(!newstring)
                    break;
                prstring = tprintf("stdout %s", newstring);
                result = pfcn(prstring, ng_ident, userptr);
                tfree(newstring);
                tfree(prstring);
            }
            /* copy the rest following \n, but without trailing \n to new address */
            if (outstringout && *outstringout != '\0')
                outstringout = copy(outstringout);
            else
                outstringout = NULL;
            tfree(delstring);
            return result;
        }
        else if (strchr(input, '\r')) {
            result = pfcn(outstringout, ng_ident, userptr);
            tfree(outstringout);
            return result;
        }
    }
    else
        myputs(input, outf);

    return 0;
}

/* provide a lock around printing function.
   May become critical if latency of callback is too high. */
int
sh_fputs(const char *input, FILE* outf)
{
    mutex_lock(&fputsMutex);
    sh_fputsll(input, outf);
    mutex_unlock(&fputsMutex);
    return 0;
}

#else

/* FIFO storage for strings created by fputs and all other printing commands.
   A string will be appended to the FIFO by fcn
   outstorage() by the main thread or the background (bg_) thread.
   Each string is read from top of the FIFO by independent thread using
   again fcn outstoraghe(), top entry is deleted and string is
   sent to caller in an endless loop by fcn printsend() */
static wordlist *wlstart = NULL, *wlend = NULL;
//static bool printstopp = FALSE;


int
sh_fputs(const char *input, FILE* outf)
{
    int result = 0;
    size_t len;
    char *delstring, *newstring, *prstring;
    size_t inputlen = strlen(input);

    /* If caller has sent NULL address for pfcn */
    if (noprintfwanted)
        return -1;

    if (outf == stderr) {
        if (!outstringerr)
            delstring = outstringerr = copy(input);
        else {
            len = strlen(outstringerr);
            delstring = outstringerr = TREALLOC(char, outstringerr, len + inputlen + 2);
            strcat(outstringerr, input);
        }
        if (strchr(input, '\n')) {
            while (outstringerr) {
                newstring = gettok_char(&outstringerr, '\n', FALSE, FALSE);
                if(!newstring)
                    break;
                prstring = tprintf("stderr %s", newstring);
                mutex_lock(&fputsMutex);
                outstorage(prstring, TRUE);
                mutex_unlock(&fputsMutex);
                tfree(newstring);
                prstring = NULL; /* keep prstring here, address is in use */
            }
            /* copy the rest following \n, but without trailing \n to new address */
            if (outstringerr && *outstringerr != '\0')
                outstringerr = copy(outstringerr);
            else
                outstringerr = NULL;
            tfree(delstring);
            return result;
        }
        else if (strchr(input, '\r')) {
            mutex_lock(&fputsMutex);
            outstorage(outstringerr, TRUE);
            mutex_unlock(&fputsMutex);
            outstringerr = NULL;
            return 0;
        }
    }
    if (outf == stdout) {
        if (!outstringout)
            delstring = outstringout = copy(input);
        else {
            len = strlen(outstringout);
            delstring = outstringout = TREALLOC(char, outstringout, len + inputlen + 1);
            strcat(outstringout, input);
        }
        if (strchr(input, '\n')) {
            while (outstringout) {
                newstring = gettok_char(&outstringout, '\n', FALSE, FALSE);
                if(!newstring)
                    break;
                prstring = tprintf("stdout %s", newstring);
                mutex_lock(&fputsMutex);
                outstorage(prstring, TRUE);
                mutex_unlock(&fputsMutex);
                tfree(newstring);
                prstring = NULL;
            }
            /* copy the rest following \n, but without trailing \n to new address */
            if (outstringout && *outstringout != '\0')
                outstringout = copy(outstringout);
            else
                outstringout = NULL;
            tfree(delstring);
            return result;
        }
        else if (strchr(input, '\r')) {
            mutex_lock(&fputsMutex);
            outstorage(outstringout, TRUE);
            mutex_unlock(&fputsMutex);
            outstringout = NULL;
            return result;
        }
    }
    else
        myputs(input, outf);

    return 0;
}

static char *outsend = NULL;

/* Endless loop in its own thread for reading data from FIFO and sending to caller */
static void
printsend(void)
{
    ps_exited = FALSE;
    printstopp = FALSE;
    for (;;) {
#if defined(__MINGW32__) || defined(_MSC_VER)
        Sleep(50);  // loop delay
#else
        usleep(50000);
#endif
        if (printstopp) { // issued by shared_exit()
            // catch the final error message
            mutex_lock(&fputsMutex);
            outsend = outstorage(NULL, FALSE);
            mutex_unlock(&fputsMutex);

            break;
        }
        mutex_lock(&fputsMutex);
        outsend = outstorage(NULL, FALSE);
        mutex_unlock(&fputsMutex);
        if (outsend) {
            /* requires outsend to be copied by the caller,
            because it is freed immediately */
            pfcn(outsend, userptr);
            tfree(outsend);
        }
    }
    ps_exited = TRUE;
}

/* remove the first entry of a wordlist, but keep wl->wl_word */
void wl_delete_first(wordlist **wlstart, wordlist **wlend)
{
    wordlist *wl_temp;

    if (!(*wlstart))
        return;
    if ((*wlstart) && !((*wlstart)->wl_next)) {
        tfree(*wlstart); /* keep wlstart->wl_word */
        (*wlstart) = NULL;
        (*wlend) = NULL;
        return;
    }
    wl_temp = (*wlstart)->wl_next;
    wl_temp->wl_prev = NULL;
    tfree(*wlstart); /* keep wlstart->wl_word */
    (*wlstart) = wl_temp;
}


/* create a wordlist FIFO using global static variables wlstart and wlend.
   wordin has to be malloced on the heap */
char* outstorage(char* wordin, bool write)
{
    char *wordout = NULL;

    if(write)
        wl_append_word(&wlstart, &wlend, wordin);
    else if (wlstart) {
        wordout = wlstart->wl_word;
        wl_delete_first(&wlstart, &wlend);
    }
    return wordout;
}

#endif


/* New progress report to statfcn().
   An update occurs only every DELTATIME milliseconds.
   We may have two threads: main and bg_run */
#define DELTATIME 150
void SetAnalyse(
   const char * Analyse, /*in: analysis type */
   int DecaPercent /*in: 10 times the progress [%]*/
   /*HWND hwAnalyse, in: global handle to analysis window */
) {
    /* If caller has sent NULL address for statfcn */
    if (nostatuswanted)
        return;

    /* check in which thread we are in */
    static unsigned int ng_id1 = 0, ng_id2 = 0;
    bool thread1;

#if defined (USE_OMP) \
|| defined (HAVE_QUERYPERFORMANCECOUNTER) \
|| defined (HAVE_CLOCK_GETTIME) \
|| defined (HAVE_GETTIMEOFDAY) \
|| defined (HAVE_TIMES) \
|| defined (HAVE_GETRUSAGE) \
|| defined (HAVE_FTIME)
    PerfTime timenow;               /* actual time stamp */
    int diffsec, diffmillisec;      /* differences actual minus prev. time stamp */
    int result;                     /* return value from callback function */
    char* s;                        /* outputs to callback function */
    int OldPercent;                 /* Previous progress value */
    char OldAn[128];                /* Previous analysis type */
    char olds[128];                 /* previous output */
    static PerfTime timebefore;     /* previous time stamp */

    /* thread 1 */
    static int OldPercent1 = -2;     /* Previous progress value */
    static char OldAn1[128];         /* Previous analysis type */
    static char olds1[128];          /* previous output */
    static PerfTime timebefore1;     /* previous time stamp */
    /* thread2 */
    static int OldPercent2 = -2;     /* Previous progress value */
    static char OldAn2[128];         /* Previous analysis type */
    static char olds2[128];          /* previous output */
    static PerfTime timebefore2;     /* previous time stamp */

    /*set the two thread ids */
    unsigned int ng_idl = threadid_self();
    if (ng_id1 == 0) {
        ng_id1 = ng_idl;
        strncpy(OldAn1, Analyse, 127); //strcpy(OldAn1, "?"); /* initial value */
    }
    else if (ng_id2 == 0 && ng_id1 != ng_idl) {
        ng_id2 = ng_idl;
        strncpy(OldAn2, Analyse, 127); // strcpy(OldAn2, "?"); /* initial value */
    }

    if (ng_idl == ng_id1) {
        thread1 = TRUE;
        strcpy(OldAn, OldAn1);
        strcpy(olds, olds1);
        OldPercent = OldPercent1;
        timebefore.milliseconds = timebefore1.milliseconds;
        timebefore.seconds = timebefore1.seconds;
    }
    else if (ng_idl == ng_id2) {
        thread1 = FALSE;
        strcpy(OldAn, OldAn2);
        strcpy(olds, olds2);
        OldPercent = OldPercent2;
        timebefore.milliseconds = timebefore2.milliseconds;
        timebefore.seconds = timebefore2.seconds;
    }
    else
        return;

    CKTcircuit *ckt = NULL;

    if (ft_curckt)
        ckt = ft_curckt->ci_ckt;

    if ((DecaPercent == OldPercent) && !strcmp(OldAn, Analyse))
        return;

    /* get current time */
    perf_timer_get_time(&timenow);
    timediff(&timenow, &timebefore, &diffsec, &diffmillisec);

    s = TMALLOC(char, 128);

    if (!strcmp(Analyse, "tran")) {
        if (ckt && (ckt->CKTtime > ckt->CKTfinalTime - ckt->CKTmaxStep)) {
           sprintf(s, "--ready--");
           result = statfcn(s, ng_ident, userptr);
           tfree(s);
           return;
        }
    }

    if (DecaPercent >= 1000){
        /* Because CKTmaxStep may be smaller than 0.1%, we print only when CKTtime is large enough. */
        if (!strcmp(Analyse, "tran") && ckt && (ckt->CKTtime < ckt->CKTfinalTime - ckt->CKTmaxStep)) {
            tfree(s);
            return;
        }
        sprintf( s, "--ready--");
        result = statfcn(s, ng_ident, userptr);
        tfree(s);
        return;
    }
    /* info every one percent of progress:
       actual time, progress,
       to catch linearity of progress of simulation */
    if (ft_ngdebug && !strcmp(Analyse, "tran"))
       if ((int)((double)DecaPercent/10.) > (int)((double)OldPercent/10.)) {
          printf("%3.1f%% percent progress after %4.2f seconds.\n", (double)DecaPercent/10., seconds());
       }
    if(thread1)
        OldPercent1 = DecaPercent;
    else
        OldPercent2 = DecaPercent;
    /* output only into hwAnalyse window and if time elapsed is larger than
       DELTATIME given value, or if analysis has changed, else return */
    if ((diffsec > 0) || (diffmillisec > DELTATIME) || strcmp(OldAn, Analyse)) {
        if (DecaPercent < 0) {
            sprintf( s, "--ready--");
        }
        else if (DecaPercent == 0) {
            sprintf( s, "%s", Analyse);
        }
        else if (!strcmp(Analyse, "shooting")) {
            sprintf( s, "%s: %d", Analyse, DecaPercent);
        }
        else {
            sprintf( s, "%s: %3.1f%%", Analyse, (double)DecaPercent/10.);
        }
        if (thread1) {
            timebefore1.milliseconds = timenow.milliseconds;
            timebefore1.seconds = timenow.seconds;
        }
        else {
            timebefore2.milliseconds = timenow.milliseconds;
            timebefore2.seconds = timenow.seconds;
        }
        /* info when previous analysis period has finished */
        if (strcmp(OldAn, Analyse)) {
            if ((ft_nginfo || ft_ngdebug) && (strcmp(OldAn, "")))
               printf("%s finished after %5.3f seconds.\n", OldAn, seconds());
            if(thread1)
                strncpy(OldAn1, Analyse, 127);
            else
                strncpy(OldAn2, Analyse, 127);
        }
        /* ouput only after a change */
        if (strcmp(olds, s))
            result = statfcn(s, ng_ident, userptr);
        if(thread1)
            strcpy(olds1, s);
        else
            strcpy(olds2, s);
    }
    tfree(s);
#else
    char* s;
    int result;
    static bool havesent = FALSE;
    if (!havesent) {
        s = copy("No usage info available");
        result = statfcn(s, ng_ident, userptr);
        tfree(s);
        havesent = TRUE;
    }
#endif
}

/* a dll or shared library should never exit, if loaded dynamically,
   but ask for graceful shutdown (e.g. being detached) via a callback function */
ATTRIBUTE_NORETURN void shared_exit(int status)
{
    /* alert caller to detach dll (if we are in the main thread),
    or detach after a short sleep, if immediate is true, and we are
    in a worker thread */
    if (immediate)
#if defined(__MINGW32__) || defined(_MSC_VER)
        Sleep(100); // va: windows native
#else
        usleep(10000);
#endif
    /* status >= 1000 tells us that we react on command 'quit'
      hand this information over to caller */
    if (status >= 1000) {
        coquit = TRUE;
        fprintf(stdout, "\nNote: 'quit' asks for resetting or detaching ngspice.dll.\n");
        status -= 1000;
    }
    else {
        coquit = FALSE;
        fprintf(stderr, "Error: ngspice.dll cannot recover and awaits to be reset or detached\n");
    }
#ifndef low_latency
    // set flag to stop the printsend thread
    printstopp = TRUE;
    // leave this thread for 100ms to stop the printsend thread
#if defined __MINGW32__ || defined _MSC_VER
    Sleep(100);
#else
    usleep(100000);
#endif
    // send the final error message already caught in printsend()
    if (outsend) {
        /* requires outsend to be copied by the caller,
        because it is freed immediately */
        pfcn(outsend, ng_ident, userptr);
        tfree(outsend);
    }
#endif
    // if we are in a worker thread, we exit it here
    // detaching then has to be done explicitely by the caller
    if (fl_running && !fl_exited) {
        fl_exited = TRUE;
        bgtr(fl_exited, ng_ident, userptr);
        // set a flag that ngspice wants to be detached
        if(ngexit)
            ngexit(status, FALSE, coquit, ng_ident, userptr);
        // finish and exit the worker thread
#ifdef HAVE_LIBPTHREAD
        pthread_exit(NULL);
#elif defined _MSC_VER || defined __MINGW32__
        _endthreadex(1);
#endif
    }
    // set a flag in caller to detach ngspice.dll
    if(ngexit)
        ngexit(status, immediate, coquit, ng_ident, userptr);

    // jump back to finish the calling function
    if (!intermj)
        longjmp(errbufm,1); /* jump back to ngSpice_Circ() */
    else
        longjmp(errbufc,1); /* jump back to ngSpice_Command() */
}

static int len = 0;
static pvecvaluesall curvecvalsall;

#ifdef olld
static pvecvalues* curvecvals;
static char type_name[128];
int sh_ExecutePerLoop_old(void)
{
    struct dvec *d;
    int i, veclen;
    struct plot *pl = plot_cur;
    /* return immediately if callback not wanted */
    if (nodatawanted)
        return 2;

    /* reset data structure if there is a change in plot type */
    if (strcmp(type_name, pl->pl_typename)) {

        if (curvecvals) {
            for (i = 0; i < len; i++)
                tfree(curvecvals[i]);
            tfree(curvecvals);
        }
        len = 0;
        memset(type_name, 0, 128);
    }

    /* initialize new for every new plot, e.g. if changed from op1 to ac1
       or from tran1 to tran2 */
    if ((pl) && (len == 0)) {
        strcpy(type_name, pl->pl_typename);
        for (d = pl->pl_dvecs; d; d = d->v_next)
            len++;

        if (len == 0) {
            fprintf(cp_err, "Error: There are no vectors currently active.\n");
            return 1;
        }
        /* allocate memory for the number of vectors */
        curvecvals = TMALLOC(pvecvalues, len);

        /* allocate memory for each entry and add vector names once */
        for (d = pl->pl_dvecs, i = 0; d; d = d->v_next, i++) {
            curvecvals[i] = TMALLOC(vecvalues, 1);
            curvecvals[i]->name = copy(d->v_name);
        }
    }
    /* get the data of the last entry to the plot vector */
    veclen = pl->pl_dvecs->v_length - 1;
    for (d = pl->pl_dvecs, i = 0; d; d = d->v_next, i++) {
        /* test if real */
        if (d->v_flags & VF_REAL) {
            curvecvals[i]->is_complex = FALSE;
            curvecvals[i]->creal = d->v_realdata[veclen];
            curvecvals[i]->cimag = 0.;
        }
        else {
            curvecvals[i]->is_complex = TRUE;
            curvecvals[i]->creal = d->v_compdata[veclen].cx_real;
            curvecvals[i]->cimag = d->v_compdata[veclen].cx_imag;
        }
    }
    /* now call the callback function to return the data to the caller */
    if (!nodatawanted)
  //      datfcn(curvecvals, len, ng_ident, userptr);

    return 0;
}
#endif

/* called each time a new data set is written to the output vectors */
int sh_ExecutePerLoop(void)
{
    struct dvec *d;
    int i, veclen;
//  double testval;
    struct plot *pl = plot_cur;
    /* return immediately if callback not wanted */
    if (nodatawanted)
        return 2;

    /* get the data of the last entry to the plot vector */
    veclen = pl->pl_dvecs->v_length - 1;
    /* safeguard against vectors with 0 length (e.g. @c1[i] during ac simulation) */
    if (veclen < 0)
        return 2;
    curvecvalsall->vecindex = veclen;
    for (d = pl->pl_dvecs, i = 0; d; d = d->v_next, i++) {
        /* test if real */
        if (d->v_flags & VF_REAL) {
            curvecvalsall->vecsa[i]->is_complex = FALSE;
//          testval = d->v_realdata[veclen];
            curvecvalsall->vecsa[i]->creal = d->v_realdata[veclen];
            curvecvalsall->vecsa[i]->cimag = 0.;
        }
        else {
            curvecvalsall->vecsa[i]->is_complex = TRUE;
            curvecvalsall->vecsa[i]->creal = d->v_compdata[veclen].cx_real;
            curvecvalsall->vecsa[i]->cimag = d->v_compdata[veclen].cx_imag;
        }
    }
    /* now call the callback function to return the data to the caller */
    datfcn(curvecvalsall, len, ng_ident, userptr);

    return 0;
}


/* called once for a new plot from beginPlot() in outitf.c,
   after the vectors in ngspice for this plot have been set.
   Transfers vector information to the caller via callback datinitfcn()
   and sets transfer structure for use in sh_ExecutePerLoop() */
int sh_vecinit(runDesc *run)
{
    struct dvec *d, *ds;
    int veccount, i;
    static pvecinfoall pvca = NULL;
    pvecinfo *pvc;

    /* return immediately if callback not wanted */
    if (nodatainitwanted)
        return 2;

    cur_run = run;

    len = veccount = cur_run->numData;

    if (veccount == 0) {
        fprintf(cp_err, "Error: There are no vectors currently active.\n");
        return 1;
    }

    /* delete the structs from the previous plot */
    if (pvca) {
        for (i = 0; i < pvca->veccount; i++)
            tfree(pvca->vecs[i]);
        tfree(pvca->vecs);
        tfree(pvca);
        pvca = NULL;
    }

    pvc = TMALLOC(pvecinfo, veccount);
    ds = cur_run->runPlot->pl_scale;
    for (i = 0, d = cur_run->runPlot->pl_dvecs; i < veccount; i++, d = d->v_next) {
        pvc[i] = TMALLOC(vecinfo, 1);
        pvc[i]->number = i;
        pvc[i]->pdvec = (void*)d;
        pvc[i]->pdvecscale = (void*)ds;
        pvc[i]->vecname = d->v_name;
        pvc[i]->is_real = (d->v_flags & VF_REAL);
    }
    pvca = TMALLOC(vecinfoall, 1);
    // the plot
    pvca->title = cur_run->runPlot->pl_title;
    pvca->date = cur_run->runPlot->pl_date;
    pvca->name = cur_run->runPlot->pl_name;
    pvca->type = cur_run->runPlot->pl_typename;
    pvca->veccount = veccount;
    // the data
    pvca->vecs = pvc;
    /* now call the callback function to return the data to the caller */
    datinitfcn(pvca, ng_ident, userptr);

    /* generate the data tranfer structure,
       data will be sent from sh_ExecutePerLoop() via datfcn() */
    if (!curvecvalsall) {
        curvecvalsall = TMALLOC(vecvaluesall, 1);
    }
    else {
        for (i = 0; i < curvecvalsall->veccount; i++)
            tfree(curvecvalsall->vecsa[i]);
        tfree(curvecvalsall->vecsa);
    }

    curvecvalsall->veccount = veccount;
    curvecvalsall->vecsa = TMALLOC(pvecvalues, veccount);

    for (i = 0, d = cur_run->runPlot->pl_dvecs; i < veccount; i++, d = d->v_next) {
        curvecvalsall->vecsa[i] = TMALLOC(vecvalues,1);
        curvecvalsall->vecsa[i]->name = d->v_name;
        if (cieq(d->v_plot->pl_scale->v_name, d->v_name))
            curvecvalsall->vecsa[i]->is_scale = TRUE;
        else
            curvecvalsall->vecsa[i]->is_scale = FALSE;
    }
    return 0;
}


/* issue callback to request external voltage data for source vname */
double
getvsrcval(double time, char *vname)
{
    double vval;
    if (!wantvdat) {
        fprintf(stderr, "Error: No callback supplied for source %s\n", vname);
        shared_exit(EXIT_BAD);
    }
    else {
        /* callback fcn */
        getvdat(&vval, time, vname, ng_ident, userptr);
        return vval;
    }
}


/* issue callback to request external current data for source iname*/
double
getisrcval(double time, char *iname)
{
    double ival;
    if (!wantidat) {
        fprintf(stderr, "Error: No callback supplied for source %s\n", iname);
        shared_exit(EXIT_BAD);
    }
    else {
        /* callback fcn */
        getidat(&ival, time, iname, ng_ident, userptr);
        return ival;
    }
}


/*
    return value 1: continue with new time step, ckt->CKTtime + ckt->CKTdelta will be
                    done next automatically.
                    For time synchronization we may choose our own ckt->CKTdelta, being
                    smaller than the one suggested by ngspice.
    return value 0: will redo the most recent time step. We may subtract olddelta and
                    continue with new ckt-CKTdelta.
                    This is necessary if non-convergence has been detected (redostep = 1).
                    The newly suggested ckt-CKTdelta has already been divided by 8.
                    This is also enforced if the truncation error is too large.
                    The newly suggested ckt-CKTdelta may be accompanied by an increase
                    of integration order.
                    For time synchronization, if the actual, converged ckt-CKTtime is
                    beyond the optimum common time, we subtract olddelta and then choose
                    our own ckt->CKTdelta, being smaller than olddelta.
    Whereas redostep is set by ngspice, the user may decide via the callback function,
    to redo the most recent step because of other reasons. This is accomplished by
    returning a 1 with the callback function.

*/

/*
    ckttime   pointer to ckt->CKTtime, which already has been used trying to achieve
              convergence, after olddelta had been added in the previous step.
    cktdelta  pointer to newly defined ckt->CKTdelta, e.g. by recognizing truncation errors
    olddelta  old ckt->CKTdelta, has already been added in the previous step.
    finalt    final time CKTfinaltime
    delmin    minimum delta CKTdelmin
    redostep  if 0, converged,
              if 1, either no convergence, need to redo with new ckt->CKTdelta
              or ckt->CKTdelta has been reduced by truncation errors too large.
    rejected  pointer to ckt->CKTstat->STATrejected, counts rejected time points.
    loc       location of function call in dctran.c: 0: after breakpoint handling, 1: at end of for loop
*/

int
sharedsync(double *pckttime, double *pcktdelta, double olddelta, double finalt,
           double delmin, int redostep, int *rejected, int loc)
{
    /* standard procedure, cktdelta has been provided by ngspice */
    if (!wantsync) {
        if (redostep) {
            *pckttime -= olddelta;
            (*rejected)++;
            return 1;
        }
        else
            return 0;
    /* synchronization required, to be done by changing cktdelta */
    } else {
        if (redostep) {
            *pckttime -= olddelta;
            (*rejected)++;
            /* use cktdelta as suggested by ngspice or acquire new cktdelta
            via pointer pcktdelta in user supplied callback */
            getsync(*pckttime, pcktdelta, olddelta, redostep, ng_ident, loc, userptr);
            /* never move beyond final time */
            if (*pckttime + *pcktdelta > finalt)
                *pcktdelta = finalt - *pckttime - 1.1 * delmin;
            return 1;
        }
        else {
            /* Use cktdelta as suggested by ngspice or acquire new cktdelta
               via pointer pcktdelta in user supplied callback. Redo the previous
               step if return value from getsync is 1. */
            int retval = getsync(*pckttime, pcktdelta, olddelta, redostep, ng_ident, loc, userptr);
            /* never move beyond final time */
            if (*pckttime + *pcktdelta > finalt)
                *pcktdelta = finalt - *pckttime - 1.1 * delmin;

            /* user has decided to redo the step, ignoring redostep being set to 0
            by ngspice. */
            if (retval) {
                *pckttime -= olddelta;
                (*rejected)++;
            }
            return retval;
        }
    }
}

#ifdef XSPICE
void shared_send_event(int index, double step, double dvalue, char *svalue, void *pvalue, int plen, int mode)
{
    if(wantevtdata)
        sendevt(index, step, dvalue, svalue, pvalue, plen, mode, ng_ident, euserptr);
    return;
}

void shared_send_dict(int index, int no_of_nodes, char* name, char*type)
{
    if (sendinitevt)
        sendinitevt(index, no_of_nodes, name, type, ng_ident, euserptr);
}
#endif

static int totalreset(void)
{

    is_initialized = FALSE;

    // if we are in a worker thread, we exit it here
    // detaching then has to be done explicitely by the caller
    if (fl_running && !fl_exited) {
        fl_exited = TRUE;
        bgtr(fl_exited, ng_ident, userptr);
    // finish and exit the worker thread
#ifdef HAVE_LIBPTHREAD
        pthread_exit(NULL);
#elif defined _MSC_VER || defined __MINGW32__
        _endthreadex(1);
#endif
    }

    /* start to clean up the mess */

    noprintfwanted = FALSE;
    nostatuswanted = FALSE;
    nodatawanted = FALSE;
    nodatainitwanted = FALSE;
    nobgtrwanted = FALSE;
    wantvdat = FALSE;
    wantidat = FALSE;
    wantsync = FALSE;
    immediate = FALSE;
    coquit = FALSE;

    wordlist all = { "all", NULL, NULL };
    wordlist star = { "*", NULL, NULL };

    tfree(Infile_Path);
    wl_free(sourceinfo);
    sourceinfo = NULL;

    com_destroy(&all);
    com_unalias(&star);
    com_undefine(&star);

    cp_remvar("history");
    cp_remvar("noglob");
    cp_remvar("brief");
    cp_remvar("sourcepath");
    cp_remvar("program");
    cp_remvar("prompt");

    destroy_wallace();

    rem_controls();

    while (ft_curckt) {
        com_remcirc(NULL);
    }

    cp_destroy_keywords();
    destroy_ivars();

    tfree(errMsg);

    destroy_const_plot();
    spice_destroy_devices();
    unset_all();
    cp_resetcontrol(FALSE);
    sh_delete_myvec();

#ifdef THREADS
    /* Destroy the mutexes */
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_destroy(&triggerMutex);
    pthread_mutex_destroy(&allocMutex);
    pthread_mutex_destroy(&fputsMutex);
    pthread_mutex_destroy(&vecreallocMutex);
    cont_condition = FALSE;
#else
#ifdef SRW
    /* Do we need to remove the SWR locks? */
//    InitializeSRWLock(&triggerMutex);
//    InitializeSRWLock(&allocMutex);
//    InitializeSRWLock(&fputsMutex);
//    InitializeSRWLock(&vecreallocMutex);
#else
    DeleteCriticalSection(&triggerMutex);
    DeleteCriticalSection(&allocMutex);
    DeleteCriticalSection(&fputsMutex);
    DeleteCriticalSection(&vecreallocMutex);
#endif
#endif
    // Id of primary thread
    main_id = 0;
#endif

    return 0;
};
