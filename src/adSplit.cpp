/* ----------------------------------------------------------------------
Program: isis
Title: Identifying Splits with clear separation - a class discovery 
  method for gene expression data
Author: Anja von Heydebreck <heydebre@molgen.mpg.de>, 
        Wolfgang Huber <w.huber@dkfz.de>
Maintainer: Anja von Heydebreck <heydebre@molgen.mpg.de>
Depends: R (>= 1.2.3)
Description: The program implements the class discovery method for 
  microarray gene expression data described in "Identifying Splits with 
  Clear Separation: A New Class Discovery Method for Gene Expression Data", 
  Anja von Heydebreck, Wolfgang Huber, Annemarie Poustka, and Martin 
  Vingron, ISMB 2001, Bioinformatics 17, suppl. 1 (2001) p. 107. 
URL: http://www.r-project.org, http://www.molgen.mpg.de/~heydebre/software.html
----------------------------------------------------------------------- */
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include "math.h"

#define CAREFUL
#undef DEBUG

//--------------------------------------------------
// a convenient inline function
//--------------------------------------------------
inline double sqr(double x) { return x*x;}

//------------------------------------------------------------
// A class to encapsulate the parameters of the score function 
// maximization algorithm
//------------------------------------------------------------
class isis_pars {
 public:
  int p;                   // dimension of discriminant space
  int p_offs;              // number of dimensions to "trim"
  int min_nrobj_split;     // minimum number of elements in a group to even 
                           // consider the split

  isis_pars(int* x) {
    p               = x[0];
    p_offs          = x[1];
    min_nrobj_split = x[2];
  }

  isis_pars() {
    p = p_offs = min_nrobj_split = 0;
  }
};

//---------------------------------------------------
// Some classes useful for being thrown as exceptions
//---------------------------------------------------
class ValueOutOfBounds {
private:
  int line, val, minval, maxval;
public:
  ValueOutOfBounds(int l, int v, int mi, int ma) {
    line=l; val=v; minval=mi; maxval=ma;
  }
  void print() {
    printf("Exception in %s, line %d: value %d out of bounds [%d...%d]\n",
	    __FILE__, line, val, minval, maxval);
  }
};
class Tomato {
private:
  char message[256];
  int line;
public:
  Tomato(int l, char* msg) { 
    line=l; 
    if (strlen(msg) > 255) msg[255]=0;
    strcpy(message, msg);
  }
  void print() {
    printf("Exception in %s, line %d: %s\n", __FILE__, line, message);
  }
};
//------------------------------------------------------------
// A very simple class to represent "splits", i.e. binary
// partitions 
//------------------------------------------------------------
class split {
 private:
  char* elements;
 public:
  static int nrobj;         

  split();                  // default constructor
  split(split& s);          // constructor from another split
  ~split();                 // destructor
  split operator=(split s); // assignment
  int operator[](int i);    // return the class membership of element i 
                            // (result is 0 or 1)
  void set(int i, int cl);  // set the class membership of element i (to 
                            // 0 or 1) (using "set" rather than a version 
                            // of "operator=" allows for a somewhat more 
                            // straightforward checking whether "cl" has 
                            // a legal value.)
  void flip(int i);         // flip the class membership of element i
  int n1();                 // number of elements of class "1"
}; 


// the common instance of the static variable "nrobj"
int split::nrobj; 

// constructors
split::split() { 
  elements = new char[nrobj]; 
}
split::split(split& s) { 
  elements = new char[nrobj]; 
  for (int i=0; i<nrobj; i++) elements[i] = s[i];
}
// destructor
split::~split() { 
  delete[] elements; 
}

int split::operator[](int i) { 
#ifdef CAREFUL
  if ((i<0) || (i>=nrobj)) throw ValueOutOfBounds(__LINE__, i, 0, nrobj-1); 
#endif
  return (int) elements[i];
}

void split::set(int i, int cl) { 
#ifdef CAREFUL
  if ((i<0) || (i>=nrobj)) throw ValueOutOfBounds(__LINE__,  i, 0, nrobj-1); 
  if ((cl<0) || (cl>1))    throw ValueOutOfBounds(__LINE__, cl, 0, 1); 
#endif
  elements[i] = cl;
}

split split::operator=(split s) { 
  for (int i=0; i<nrobj; i++) elements[i] = s[i];
#ifdef CAREFUL
  for (int i=0; i<nrobj; i++) 
    if ((elements[i]!=0) && (elements[i]!=1)) throw ValueOutOfBounds(__LINE__, elements[i], 0, 1); 
#endif
  return *this;
}

void split::flip(int i) { 
#ifdef CAREFUL
  if ((i<0) || (i>=nrobj)) throw ValueOutOfBounds(__LINE__,  i, 0, nrobj-1); 
#endif
  elements[i] = 1-elements[i];
}

int split::n1()
{
#ifdef CAREFUL
  for (int i=0; i<nrobj; i++) 
    if (((*this)[i]!=0) && ((*this)[i]!=1))
      throw ValueOutOfBounds(__LINE__, (*this)[i], 0, 1);
#endif
  int n1 = 0;
  for (int i=0; i<nrobj; i++) 
    if ((*this)[i]==1) n1++;
  return n1;
}

// ------------------------------------------------------------
// a datastructure which is convenient for keeping track of the
// array elements' indices when using qsort()
// ------------------------------------------------------------
struct doublewithindex {
   double value;
   int index;
};

// ------------------------------------------------------------
// compare functions for use with qsort()
// ------------------------------------------------------------
int compare_ascending(const void *x, const void *y) {
  double d =  ((doublewithindex*)x) -> value - 
              ((doublewithindex*)y) -> value ;
  return (d>0) ? 1 : ((d<0) ? -1 : 0);
}
int compare_descending(const void *x, const void *y) {
  double d =  ((doublewithindex*)x) -> value - 
              ((doublewithindex*)y) -> value ;
  return (d<0) ? 1 : ((d>0) ? -1 : 0);
}

//--------------------------------------------------------------------
// ttesttwo
// INPUT
// data:      an array of doubles, storing the nrow x ncol data matrix
//            column-wise (i.e. row is the fast index)
// s:         the split
// par:       parameters
//
// OUTPUT
// t:       an array of nrow doubles, to be filled with the t-values
//          memory for the array is assumed to be already allocated
//--------------------------------------------------------------------
void ttesttwo(double* data, int nrow, int ncol, split& s, isis_pars& par, double* t) {
  // number of elements in the two groups
  int n1 = s.n1();
  int n0 = split::nrobj - n1;
  const double factor = (1.0/n0+1.0/n1) / (split::nrobj - 2.0);

  if ((n0 < par.min_nrobj_split) || (n1 < par.min_nrobj_split))
    throw ValueOutOfBounds(__LINE__, n1, par.min_nrobj_split, split::nrobj-par.min_nrobj_split);
 
  for (int r=0; r < nrow; r++) {
    double  x[2] = {0.0, 0.0};   // first moments
    double xx[2] = {0.0, 0.0};   // second moments
    for (int i=0; i<ncol; i++) {
       x[s[i]] +=     data[ i + ncol*r ];
      xx[s[i]] += sqr(data[ i + ncol*r ]);
    }    
    double v0 = (xx[0] - sqr(x[0])/n0);
    double v1 = (xx[1] - sqr(x[1])/n1);
    t[r] =  (x[0]/n0 - x[1]/n1) / sqrt( (v0+v1) * factor ); 
  } // for r
  return;
}

//--------------------------------------------------------------------
// tscore
// INPUT
// data:    an array of doubles, storing the nrow x ncol data matrix
//          column-wise (i.e. row is the fast index)
// s:       the split
// par:     parameters
//
// OUTPUT
// score:   the t-score
//--------------------------------------------------------------------
double tscore(double* data, int nrow, int ncol, split& s, isis_pars& par) {

  // number of elements in the two groups
  int n1 = s.n1();
  int n0 = split::nrobj - n1;
  doublewithindex* twi = new doublewithindex[nrow];

  // Construct an array "list" with the indices of group elements of the smaller group
  int nlist = (n0<=n1) ? n0 : n1;
  int clist = (n0<=n1) ?  0 : 1;
  int* list = new int[nlist];
  int k    = 0;
  for (int i=0; i<ncol; i++) 
    if (s[i]==clist) {
      list[k] = i;
      k++;
    }
  if (k != nlist) 
    throw Tomato(__LINE__, "tscore: internal error, shame on the programmer");

  // For standardized data, selecting by the t-statistic is equivalent to calculating
  // the mean of values in one group, and selecting the rows with the
  // highest absolute value of that. The latter is faster.
  for (int r=0; r < nrow; r++) {
    double sum = 0;
    for (k=0; k < nlist; k++)
      sum += data[ncol*r + list[k]];
    twi[r].value = fabs(sum);
    twi[r].index = r;
  }
  delete[] list;

  // we need the highest scores generated
  qsort(twi, nrow, sizeof(doublewithindex), compare_descending);

  // a: the discriminant axis
  double* a = new double[par.p];  
  // loop over j: discrimant space dimensions
  for (int j=par.p_offs; j < par.p; j++) {
    int r = twi[j].index;
    double  x[2] = {0.0, 0.0};   // first moments
    double xx[2] = {0.0, 0.0};   // second moments
    for (int i=0; i<ncol; i++) {
       x[s[i]] +=     data[ i + ncol*r ];
      xx[s[i]] += sqr(data[ i + ncol*r ]);
    }    
    double v0 = (xx[0] - sqr(x[0])/n0);
    double v1 = (xx[1] - sqr(x[1])/n1);
    a[j] = (x[0]/n0 - x[1]/n1) / (v0+v1); 
  } // for j 

  // cs: the coordinates of the samples on the discriminant axis
  double* cs = new double[ncol];  
  for (int i=0; i < ncol; i++) cs[i] = 0;

  // loop over j: discrimant space dimensions
  for (int j=par.p_offs; j<par.p; j++) {  
    int r = twi[j].index;
    // loop over i: samples
    for (int i=0; i < ncol; i++) {
      cs[i] += a[j] * data[ i + ncol*r ];
    }  // for i
  } // for j

  double score=0;
  ttesttwo(cs, 1, ncol, s, par, &score);

  delete[] cs;
  delete[] a;
  delete[] twi;

  return fabs(score);
}

//---------------------------------------------------------------------------
// A version of tscore for arrays of splits, with an interface very 
// similar to gotomax
//---------------------------------------------------------------------------
void tscore(double* data, int nrow, int ncol, 
            split* s, int nsplit, 
            isis_pars& par, double* score) 
{
  for (int j=0; j<nsplit; j++)
    score[j] = tscore(data, nrow, ncol, s[j], par);
  return;
}

//---------------------------------------------------------------------------
// gotomax
//
// Find local maxima of the score function by steepest ascent hillclimbing,
// starting from initial values given in "s". The hillclimbing path is
// not allowed to go to splits in which one of the parts has less than
// par.min_nrobj_split elements.
//
// INPUT
// data:    an array of doubles, storing the nrow x ncol data matrix
//          column-wise (i.e. row is the fast index)
// s:       an array of splits, of length nsplit
// par:     parameters
//
// OUTPUT
// s:       each input split is replaced by the local maximum that 
//          the hill-climbing converges to
// score:   array of DLD scores at the maxima.
//          this function assumes that memory (nsplit*sizeof(double)) 
//          is allocated!
//----------------------------------------------------------------------
void gotomax(double* data, int nrow, int ncol, 
             split* s, int nsplit, 
             isis_pars& par, double* score) 
{
  split* ms = new split[nsplit];
  double mscore, score_cand, score_best;
  int i_best;

  // For the purpose of calling tscore during the hillclimbing, allow also
  // bipartitions of smaller size. However, we would not accept those as a
  // local maximum.
  isis_pars tscorepar;
  tscorepar.p                =  par.p;
  tscorepar.p_offs           =  par.p_offs; 
  tscorepar.min_nrobj_split  =  par.min_nrobj_split-1;

  for (int j=0; j<nsplit; j++) {
    ms[j]  = s[j];
    mscore = tscore(data, nrow, ncol, ms[j], tscorepar);
    do {
      score_best    = mscore;
      i_best        = -1;
      for (int i=0; i<split::nrobj; i++) {
	ms[j].flip(i);
        int n1 = ms[j].n1();
	if ((n1 >= tscorepar.min_nrobj_split) & (n1 <= split::nrobj - tscorepar.min_nrobj_split)) { 
	  score_cand = tscore(data, nrow, ncol, ms[j], tscorepar);
	  if (score_cand > score_best) {
	    i_best     = i;
	    score_best = score_cand;
	  }
	}
	ms[j].flip(i);
      }
      if (i_best >= 0) {
	ms[j].flip(i_best);
	mscore = score_best;
      }
    } while (i_best >= 0);
    
    int n1 = ms[j].n1();
    if ((n1 >= par.min_nrobj_split) & (n1 <= split::nrobj-par.min_nrobj_split)) { 
      // if local maximum that we found is within allowed range:
      score[j] = mscore;
    } else {
      // if local maximum is outside allowed range, set everything to 0.
      for (int i=0; i<split::nrobj; i++) ms[j].set(i, 0);
      score[j] = 0;
    }
  } // for j

  // sort
  doublewithindex* swi = new doublewithindex[nsplit];
  for (int j=0; j < nsplit; j++) {
    swi[j].value = score[j];
    swi[j].index = j;
  }
  qsort(swi, nsplit, sizeof(doublewithindex), compare_descending);
  for (int j=0; j < nsplit; j++) {
    score[j] = swi[j].value;
    s[j]     = ms[swi[j].index];
  }
  delete[] ms;
  return;
}


//-------------------------------------------------------------
// external interface e.g. to R
//-------------------------------------------------------------
extern "C" void isis(char** func, 
                     double* data, int* nrow, int* ncol, 
                     int* sp, int* nsplit, 
                     int* parameters, 
                     double* t,                   
                     int* returncode) {

  isis_pars par(parameters);
  *returncode = 0;

#ifdef DEBUG
  printf("nrow=%d, ncol=%d, nsplit=%d\n", *nrow, *ncol, *nsplit);
  int i, r;
  for (int j=0; j< *nsplit; j++) {
    for (i=0; i < *ncol; i++) printf("%d ", sp[(*nsplit)*i + j]); printf("\n"); 
  }
  printf("p=%d, p_offs=%d, min_nrobj_split=%d\n",  par.p, par.p_offs, par.min_nrobj_split);
  r=0;         i=0;         printf("data[%4d, %4d] = %f\n", r, i, data[ i + (*ncol)*r ]); 
  r=0;         i=(*ncol)-1; printf("data[%4d, %4d] = %f\n", r, i, data[ i + (*ncol)*r ]); 
  r=(*nrow)-1; i=0;         printf("data[%4d, %4d] = %f\n", r, i, data[ i + (*ncol)*r ]); 
  r=(*nrow)-1; i=(*ncol)-1; printf("data[%4d, %4d] = %f\n", r, i, data[ i + (*ncol)*r ]); 
#endif
  
  split::nrobj = *ncol;
  split* s = new split[*nsplit];

  try {
    // Copy char array "sp" into array of class split "s", and check validity
    // Note: R stores array column-wise (i.e. j is the fast index)
    for (int j=0; j < *nsplit; j++)
      for (int i=0; i < *ncol; i++) 
	s[j].set(i, sp[(*nsplit)*i + j]);  

    /*---------- ttesttwo ----------*/
    if (strcmp(*func, "ttesttwo")==0) {
      if (*nsplit != 1)
        throw ValueOutOfBounds(__LINE__, *nsplit, 1, 1);
      ttesttwo(data, *nrow, *ncol, s[0], par, t);

    /*---------- tscore ----------*/
    } else if (strcmp(*func, "tscore")==0) {
	tscore(data, *nrow, *ncol, s, *nsplit, par, t);

    /*---------- gotomax ----------*/
    } else if (strcmp(*func, "gotomax")==0) {
      gotomax(data, *nrow, *ncol, s, *nsplit, par, t);

      // Copy split array "s" back into char array "sp"
      for (int j=0; j < *nsplit; j++)
        for (int i=0; i < *ncol; i++) 
	  sp[(*nsplit)*i + j] = s[j][i];

    } else
      throw Tomato(__LINE__, "isis() called with unknown function");

  } // try

  catch (ValueOutOfBounds x) {
    x.print();
    *returncode = -1;
  }
  catch(Tomato x) {
    x.print();
    *returncode = -2;
  }

  delete[] s;
  return;
}












