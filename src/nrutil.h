#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#include <Rcpp.h>

Rcpp::DataFrame trstat(int *ncol, int table[], int uwt[], int *sampsz, int *rowm,
                       int *m_trt_switch, int *obsstat, int *obscorr_0, int *obscorr_1,
                       double *pval);

typedef struct arc_tag ARC;
typedef struct snd_tag SUBND;
typedef struct rec_tag REC;
typedef struct node_tag NODE;

int imax(int x, int y);
int imin(int x, int y);
void nrerror(char error_text[]);
int *ivector(long nl, long nh);
double *dvector(long nl, long nh);

void free_ivector(int *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);

void faclog(int sampsz, double *fact);
void calnds(int ncol, int *table, int rowm, int *nnodes, int *colm, int *cumcol,
            int *stpos, int *minvl);
NODE *NDvector(long nl, long nh);
void forind(int ncol, int rowm, int nnodes, int sampsz, int *colm, int *cumcol,
            double *fact, int *stpos, int *minvl, NODE *nodes);
void backind(int ncol, int rowm, int nnodes, int sampsz, int m_trt_switch,
             int obscorr_0, int obscorr_1, int *colm, int *cumcol, int *wt,
             double *fact, int *stpos, int *minvl, NODE *nodes);
void finalpass(int nnodes, int ncol, int *minvl, int *stpos, int obsstat, int lowval,
               int obscorr_0, int obscorr_1, NODE *nodes, double *rtail);
void free_NDvector(NODE *v, long nl, long nh);
void free_SNvector(SUBND *v, long nl, long nh);
void free_arc(ARC *arc);
void freerec(REC *rec);
REC *crerec(int r, double pr); // create record
double addlog(double num1, double num2);

SUBND *SNvector(long nl, long nh);
void init(int lower, int upper, SUBND *subnodes);
SUBND *crenode(int remcorr_0, int remcorr_1, int lp, int sp, int ipl,
               double prarc); // create node
double prbar0(int coltot, int numsucc, double *fact);
void dropnd(SUBND *cursnode);
void corrlpsp(int sampsz, int parm, int stage, int ncol, int rowm, int nnodes,
              int *colm, int *cumcol, int *stpos, int *minvl, NODE *nodes, int *spl,
              int *lpl);
double testmax(double x, double y);

#endif /* _NR_UTILS_H_ */
