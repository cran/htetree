/*
 * commom variables for the causalTree routine
 *
 * Start with things that depend on R.h
 */
#include <R.h>
#include <Rinternals.h>
#include <math.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("causalTree", String)
#else
#define _(String) (String)
#endif

/*
 * Memory defined with R_alloc is removed automatically
 *  That with "CALLOC" I have to remove myself.  Use the
 *  latter for objects that need to persist between the
 *  s_to_ct1 and s_to_ct2 calls
 */
#define ALLOC(a,b)  R_alloc(a,b)
#define CALLOC(a,b) R_chk_calloc((size_t)(a), b)
#define RPARTNA(a) ISNAN(a)

/* done with the R internals */
#define LEFT  (-1)              /*used for the variable "extra" in nodes */
#define RIGHT  1
#define MISSING 0

#ifdef MAINRP
#define EXTERN
#else
#define EXTERN extern
#endif

#ifndef DEBUG
#define DEBUG 1
#endif

EXTERN struct {
    double complexity;
    double alpha;
    double iscale;              /* used to check improvement==0, with error */
    double **ydata;
    double **xdata;
    double *xtemp;
    double *wt;
    double *treatment;
    double **ytemp;
    double *wtemp;              /* temp vector of weights */
    double *trtemp;             /* temp vector of treatment status */
    double *lwt;
    double *ltr;
    double *rwt;                /*scratch double vectors, of length ncat */
    double *rtr;
    double *vcost;              /* variable costs */
    int *numcat;                /* variable type: 0=cont, 1+  =#categories */
    int **sorts;                /* matrix of sort indices */
    int n;                      /* total number of subjects  */
    int num_y;                  /* number of y variables */
    int nvar;                   /* number of predictors */
    int maxpri;
    int maxsur;                 /* max # of primary or surrogate splits to use */
    int usesurrogate;
    int num_unique_cp;
    int min_node;               /* minimum size for any terminal node */
    int min_split;              /*minimum size before we attempt a split */
    int num_resp;               /*length of the response vector */
    int sur_agree;              /*0 =  my style, 1=CART style */
    int maxnode;                /*controls the maximum depth of the tree */
    int *tempvec;               /*to be allocated by the mainline, of length n */
    int *which;
    int *csplit;
    int *left;
    int *right;
    double max_y;               /* maximum absolute value of y */
    double *xvar;               /* variance of predictors: for distance measure in matching method */
    double propensity;          /* propensity score used in this causal Tree */
    int NumHonest;              /* NumHonest for CT-H cross-validation function*/
    int NumXval;                /* number of cross validation data sets */
} ct;

// EXTERN struct cptable *cptable_tail;
// EXTERN int (*ct_init) ();       /*called to initialize a splitting function */
// EXTERN void (*ct_choose) ();    /*set to the splitting function */
// EXTERN void (*ct_eval) ();      /*set to the evaluation routine */
// EXTERN double (*ct_error) ();   /*set to the prediction error routine */
// EXTERN double (*ct_xeval) ();
// EXTERN int nodesize;

EXTERN struct cptable *cptable_tail;
EXTERN int (*ct_init) (int n, double *y[], int maxcat, char **error,
            int *size, int who, double *wt, double *treatment,
            int bucketnum, int bucketMax, double *train_to_est_ratio);       /*called to initialize a splitting function */
EXTERN void (*ct_choose) (int n,double *y[],double *x,int nclass,int edge,double *improve,double *split,int *csplit,double myrisk,double *wt,double *treatment, ...);
EXTERN void (*ct_eval) (int n,double *y[],double *value,double *tr_mean,double *con_mean,double *risk,double *wt,double *treatment,double max_y, ...);
EXTERN double (*ct_error) (char);   /*set to the prediction error routine */
EXTERN double (*ct_xeval) (double *y, ...);
EXTERN int nodesize;
// EXTERN void (*ct_choose) (int n, double *y[], double *x, int nclass, int edge, double *improve, double *split, int *csplit, double myrisk, double *wt, double *treatment, double propensity, int minsize, double alpha, double train_to_est_ratio, int bucketnum, int bucketMax);    /*set to the splitting function */
// EXTERN void (*ct_eval) (int n, double *y[], double *value, double *tr_mean, double *con_mean, double *risk, double *wt, double *treatment, double max_y, double propensity, double alpha, double train_to_est_ratio);      /*set to the evaluation routine */
// EXTERN double (*ct_xeval) (double *y, double wt, double treatment, double tr_mean, double con_mean, double trs, double cons, double alpha, double xtrain_to_est_ratio, double propensity, double xtrain_to_est_ratio, double gamma);

/*
 * The user inputs his complexity parameter as a percentage. and the
 *   printout is also scaled in this way.  The book and the computations all
 *   have an easier time with absolute cp.  So complex = what the user
 *   typed and alpha = complex * (risk of top node) = what is used
 *   internally.
 * The variable 'complex' in node.h is also on the internal scale.
 *
 * Categorical variables must be coded as 1,2,3, ..., and there may be
 *  missing categories.  The upper limit is determined on the fly.
 */
