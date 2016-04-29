/*
 * numparam.h
 */

/*** interface to spice frontend  subckt.c ***/

#include "numpaif.h"
#include "ngspice/hash.h"

/***** numparam internals ********/

typedef enum {Nodekey = '#'} _nNodekey;  /* Introduces node symbol */
typedef enum {Intro   = '&'} _nIntro;    /* Introduces preprocessor tokens */
typedef enum {Comment = '*'} _nComment;  /* Spice Comment lines */
typedef enum {Psp     = '{'} _nPsp;      /* Ps expression */


/* -----------------------------------------------------------------
 * I believe the entry_t should be a union of type but I need more info.
 * ----------------------------------------------------------------- */

struct nupa_type;

extern const struct nupa_type S_nupa_real;
extern const struct nupa_type S_nupa_string;
extern const struct nupa_type S_nupa_subckt;
extern const struct nupa_type S_nupa_unknown;
extern const struct nupa_type S_nupa_model;

#define  NUPA_REAL     (&S_nupa_real)
#define  NUPA_STRING   (&S_nupa_string)
#define  NUPA_SUBCKT   (&S_nupa_subckt)
#define  NUPA_UNKNOWN  (&S_nupa_unknown)
#define  NUPA_MODEL    (&S_nupa_model)

typedef const struct nupa_type *nupa_type;


typedef struct entry_s {
    nupa_type tp;      /* type: I)nt R)eal S)tring F)unction M)acro P)ointer */
    char *symbol;
    int  level;                 /* subckt nesting level */
    double vl;                  /* float value if defined */
    int  ivl;                   /* int value or string buffer index */
    char *sbbase;               /* string buffer base address if any */
} entry_t;


typedef struct {                /* the input scanner data structure */
    SPICE_DSTRING srcfile;      /* last piece of source file name */
    SPICE_DSTRING lookup_buf;   /* useful temp buffer for quick symbol lookup */
    int srcline;
    int oldline;
    int errcount;
    int max_stack_depth;        /* alloced maximum depth of the symbol stack */
    int stack_depth;            /* current depth of the symbol stack */
    NGHASHPTR *symbols;         /* stack of scopes for symbol lookup */
                                /*  [0] denotes global scope */
    NGHASHPTR inst_symbols;     /* instance qualified symbols - after a pop */
    char **inst_name;           /* name of subcircuit */
    char **dynrefptr;
    char *dyncategory;
    int hs_compatibility;       /* allow extra keywords */
} dico_t;


void initdico(dico_t *);
int donedico(dico_t *);
void dico_free_entry(entry_t *);
bool defsubckt(dico_t *, struct card *, nupa_type categ);
int findsubckt(dico_t *, char *s);
bool nupa_substitute(dico_t *, char *s, char *r, bool err);
bool nupa_assignment(dico_t *, char *s, char mode);
bool nupa_subcktcall(dico_t *, char *s, char *x, char *inst_name);
void nupa_subcktexit(dico_t *);
dico_t *nupa_fetchinstance(void);
entry_t *entrynb(dico_t *dico, char *s);
entry_t *attrib(dico_t *, NGHASHPTR htable, char *t, char op);
void del_attrib(void *);
