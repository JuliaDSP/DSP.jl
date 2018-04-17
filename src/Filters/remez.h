#ifndef _SCIPY_PRIVATE_SIGNAL_SIGTOOLS_H_
#define _SCIPY_PRIVATE_SIGNAL_SIGTOOLS_H_


#define BOUNDARY_MASK 12
#define OUTSIZE_MASK 3
#define FLIP_MASK  16
#define TYPE_MASK  (32+64+128+256+512)
#define TYPE_SHIFT 5

#define FULL  2
#define SAME  1
#define VALID 0

#define CIRCULAR 8
#define REFLECT  4
#define PAD      0 

#define MAXTYPES 21


/* Generally useful structures for passing data into and out of
   subroutines.  Used in the generic routines instead of the
   Python Specific structures so that the routines can be easily
   grabbed and used in another scripting language */

typedef int intp;
typedef unsigned int uintp;

typedef struct {
  char *data;
  int elsize;
} Generic_ptr;

typedef struct {
  char *data;
  intp numels;
  int elsize;
  char *zero;        /* Pointer to Representation of zero */
} Generic_Vector;

typedef struct {
  char *data;
  int  nd;
  intp  *dimensions;
  int  elsize;
  intp  *strides;
  char *zero;         /* Pointer to Representation of zero */
} Generic_Array;

typedef void (MultAddFunction) (char *, intp, char *, intp, char *, intp *, intp *, int, intp, int, intp *, intp *, uintp *);


void
scipy_signal_sigtools_linear_filter_module_init(void);

int pre_remez(double *h2, int numtaps, int numbands, double *bands,
              double *response, double *weight, int type, int maxiter,
              int grid_density);

void simple_test(int xlen, double *x);

/*
static int index_out_of_bounds(int *, int *, int );
static long compute_offsets (unsigned long *, long *, int *, int *, int *, int *, int);
static int increment(int *, int, int *);
static void convolveND(Generic_Array *, Generic_Array *, Generic_Array *, MultAddFunction *, int);
static void RawFilter(Generic_Vector, Generic_Vector, Generic_Array, Generic_Array, Generic_Array *, Generic_Array *, BasicFilterFunction *, int);
*/

#endif
