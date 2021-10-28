#include <Rcpp.h>
#include <RProgress.h>
#include <stdio.h>
#include <stdlib.h>

#include "nrutil.h"

using namespace Rcpp;

#define NR_END 1
#define FREE_ARG char *

struct rec_tag {
    int pastr;    // past record
    double pastp; // past probability
    REC *nextrec;
};

struct arc_tag {
    int arc;
    double pr; // probability length
    ARC *nextarc;
    SUBND *child;
};

struct snd_tag {
    int parcorr_0; // partial sum for the quadratic term when x = 0
    int parcorr_1; // partial sum for the quadratic term when x = 1
    int spl;       // shortest path length (for the test statistic)
    int lpl;       // longest path length (for the test statistic)
    double tp;     // total probability
    int numpred;   // number of predecessor nodes created
    ARC *arc;
    REC *rec;
};

struct node_tag {
    int splcorr;
    // shortest path length when working backwards for quadratic term
    int lplcorr;
    // longest path length when working backwards for quadratic term
    int upper;
    int lower;
    int numsucc;
    SUBND *subnodes;
};

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
    Rcout << "Numerical Recipes run-time error...\n";
    Rcout << "...now exiting to system...\n";
    stop(error_text);
}

int *ivector(long nl, long nh) // integer vector
/* allocate an int vector with subscript range v[nl..nh] */
{
    int *v;

    v = (int *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
    char message[] = "allocation failure in ivector()";
    if (!v) nrerror(message);

    return v - nl + NR_END;
}

double *dvector(long nl, long nh) // double vector
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;

    v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
    char message[] = "allocation failure in dvector()";
    if (!v) nrerror(message);

    return v - nl + NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    free((FREE_ARG)(v + nl - NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
    free((FREE_ARG)(v + nl - NR_END));
}

/* (C) Copr. 1986-92 Numerical Recipes Software #(k3#12#1i.. */

int imax(int x, int y) {
    if (x > y)
        return (x);
    else
        return (y);
}

int imin(int x, int y) {
    if (x < y)
        return (x);
    else
        return (y);
}

DataFrame trstat(int *ncol, int table[], int uwt[], int *sampsz, int *rowm,
                 int *m_trt_switch, int *obsstat, int *obscorr_0, int *obscorr_1,
                 double *pval) {
    int i, nnodes, high, lowval;
    // lowval = test statistic value (for p-value calculation, add
    // up everything with test stats larger than this value)
    int *colm, *cumcol, *wt, *stpos, *minvl, *table2;
    // colm = num responses in ith cluster?, wt = score for that
    // dose (x-value), stpos = starting position at each stage
    // table2 = sample size for each cluster, stored backwards
    // (ex: for 3 clusters, if colm = [0,2,3,3], then sample size
    // for cluster 1 = 3, cluster 2 = 3, cluster 3 = 2) minvl =
    // min number of responses at that stage (for terminal node
    // (last stage), it should be equal to the first sufficient
    // statistic (3) - in IntroductionExample.txt and it should be
    // zero for all other stages)
    double *fact, normcon; // , mean, var;  // normcon = normalizeng constant -
    // total probability for initializng node, mean and
    // variance only for the asymptotic distribution
    REC *currec;
    SUBND *cursnode;
    NODE *nodes;

    if (*rowm == 0 || *rowm == *sampsz) {
        *pval = 1.0;
        DataFrame df;
        return df;
    }
    table2 = ivector(0, *ncol - 1);
    // vector[0, 1, 2] // numbers using "IntroductionExample.txt"
    // table2 = sample size for each cluster, stored backwards
    // (ex: for 3 clusters, if colm = [0,2,3,3], then sample size
    // for cluster 1 = 3, cluster 2 = 3, cluster 3 = 2)
    colm = ivector(1, *ncol + 1);
    cumcol = ivector(1, *ncol + 1);
    wt = ivector(1, *ncol + 1);
    stpos = ivector(1, *ncol + 1);
    minvl = ivector(1, *ncol + 1);
    fact = dvector(1, *sampsz + 1);

    faclog(*sampsz, fact);

    wt[1] = 0;
    for (i = 2; i <= *ncol + 1; i++) {
        table2[i - 2] = table[*ncol - i + 1];
        wt[i] = uwt[i - 2];
    }

    calnds(*ncol, table2, *rowm, &nnodes, colm, cumcol, stpos, minvl);
    nodes = NDvector(1, nnodes);
    forind(*ncol, *rowm, nnodes, *sampsz, colm, cumcol, fact, stpos, minvl, nodes);
    for (i = 1; i <= nnodes; i++) {
        nodes[i].upper = nodes[nnodes - i + 1].lplcorr;
        nodes[i].lower = nodes[nnodes - i + 1].splcorr;
    }

    calnds(*ncol, table, *rowm, &nnodes, colm, cumcol, stpos, minvl);
    // calculate nodes
    forind(*ncol, *rowm, nnodes, *sampsz, colm, cumcol, fact, stpos, minvl, nodes);

    backind(*ncol, *rowm, nnodes, *sampsz, *m_trt_switch, *obscorr_0, *obscorr_1, colm,
            cumcol, wt, fact, stpos, minvl, nodes);

    lowval = nodes[stpos[1]].subnodes[*obscorr_0 + *obscorr_1].spl;

    normcon = nodes[stpos[1]].subnodes[*obscorr_0 + *obscorr_1].tp;

    *pval = (double)0.0;

    // this is where the final p-value is calculated:
    finalpass(nnodes, *ncol, minvl, stpos, *obsstat, lowval, *obscorr_0, *obscorr_1,
              nodes, pval);
    // cursnode = current subnode
    cursnode = &(nodes[stpos[*ncol + 1]].subnodes[0]);
    // currec = current record - last record that is processed
    currec = cursnode->rec;

    // create 3 vectors todo: name vectors to something meaningful
    IntegerVector v1;
    NumericVector v2;
    NumericVector v3;

    while (currec != NULL) {
        v1.push_back(currec->pastr);
        v2.push_back(currec->pastp);
        v3.push_back(exp(currec->pastp - normcon));
        if (currec->pastr >= *obsstat) *pval = *pval + exp(currec->pastp - normcon);
        currec = currec->nextrec;
    }

    cursnode = &(nodes[stpos[*ncol + 1]].subnodes[0]);
    currec = cursnode->rec;
    freerec(currec);
    high = imin(nodes[stpos[*ncol + 1]].upper, *obscorr_0 + *obscorr_1);
    free_SNvector(nodes[stpos[*ncol + 1]].subnodes, nodes[stpos[*ncol + 1]].lower,
                  high);
    free_NDvector(nodes, 1, nnodes);
    free_ivector(table2, 0, *ncol - 1);
    free_ivector(colm, 1, *ncol + 1);
    free_ivector(cumcol, 1, *ncol + 1);
    free_ivector(wt, 1, *ncol + 1);
    free_ivector(stpos, 1, *ncol + 1);
    free_ivector(minvl, 1, *ncol + 1);
    free_dvector(fact, 1, *sampsz + 1);

    DataFrame dist = DataFrame::create(Named("test.stat") = v1, _["unnorm.probs"] = v2,
                                       _["norm.probs"] = v3);

    return dist;
}

void finalpass(int nnodes, int ncol, int *minvl, int *stpos, int obsstat, int lowval,
               int obscorr_0, int obscorr_1, NODE *nodes, double *rtail) {
    int pos, hlim, llim, r, stage, high, k, j;
    double npr;
    SUBND *cursnode, *succ;
    ARC *carc;
    REC *newrec, *crec, *curr, *nxt;

    newrec = crerec(0, (double)0.0);
    nodes[stpos[1]].subnodes[obscorr_0 + obscorr_1].rec = newrec;
    // new progress bar here?
    double total = (ncol);
    RProgress::RProgress pb_final(
        "Final calculations [:bar] :current/:total (:percent %)", total);
    pb_final.tick(0);

    for (stage = 1; stage <= ncol; stage++) {
        llim = minvl[stage];
        if (stage == 1)
            hlim = llim;
        else
            hlim = llim + stpos[stage - 1] - stpos[stage] - 1;
        for (pos = stpos[stage]; pos <= stpos[stage] + hlim - llim; pos++) {
            if (nodes[pos].subnodes == NULL) continue;
            high = imin(nodes[pos].upper, obscorr_0 + obscorr_1);
            for (k = nodes[pos].lower; k <= high; k++) {
                cursnode = &(nodes[pos].subnodes[k]);
                if (cursnode->parcorr_0 < 0 || cursnode->parcorr_1 < 0) continue;
                carc = cursnode->arc; // carc = current subnode arc
                while (carc != NULL) {
                    succ = carc->child;
                    crec = cursnode->rec;
                    while (crec != NULL) {
                        r = crec->pastr + carc->arc; // from Formula: t_k + r_(1,k+1)
                        // (part 2, part a)
                        if (r + succ->lpl >= lowval) {
                            // from Formula: t_k + r_(1,k+1) +
                            // LP_1(k+1, s_(1,k+1), s_(2,k+1)) <
                            // t_obs (part 2, part a - in reverse (if
                            // NOT part a, then move on))
                            if (crec->pastp < 0.0000001)
                                npr = carc->pr;
                            else
                                npr = crec->pastp + carc->pr;
                            newrec = crerec(r, npr);
                            if (succ->rec == NULL)
                                succ->rec = newrec;
                            else if (newrec->pastr < succ->rec->pastr) {
                                newrec->nextrec = succ->rec;
                                succ->rec = newrec;
                            } else {
                                curr = succ->rec;
                                nxt = curr->nextrec;
                                while (nxt != NULL) {
                                    if (r == curr->pastr ||
                                        (nxt != NULL && (r < nxt->pastr)))
                                        break;
                                    curr = nxt;
                                    nxt = curr->nextrec;
                                }
                                if (r == curr->pastr) {
                                    curr->pastp = addlog(curr->pastp, newrec->pastp);
                                    free(newrec);
                                } else {
                                    newrec->nextrec = nxt;
                                    curr->nextrec = newrec;
                                }
                            }
                        }
                        crec = crec->nextrec;
                        // move on to the next record and
                        // don't transfer this record (if
                        // part 2, part a is true)
                    }
                    carc = carc->nextarc;
                }
                if (cursnode->rec != NULL) freerec(cursnode->rec);
            }
            if (nodes[pos].subnodes != NULL) {
                for (j = nodes[pos].lower; j <= high; j++) {
                    if (nodes[pos].subnodes[j].arc != NULL)
                        free_arc(nodes[pos].subnodes[j].arc);
                }
            }
            free_SNvector(nodes[pos].subnodes, nodes[pos].lower, high);
        }
        pb_final.tick();
    }
}

REC *crerec(int r, double pr) { // create record

    REC *record;

    record = (REC *)malloc(sizeof(REC));
    record->pastr = r;  // t_k = partial rank length
    record->pastp = pr; // h_k = unnormalized sum of probabilites of all
    // partial paths with rank length r_1k
    record->nextrec = NULL;
    return (record);
}

void freerec(REC *rec) {
    void freerec(REC * rec);

    if (rec->nextrec != NULL) freerec(rec->nextrec);
    free(rec);
    return;
}

void free_arc(ARC *arc) {
    void free_arc(ARC * arc);

    if (arc->nextarc != NULL) free_arc(arc->nextarc);
    free(arc);
    return;
}

void backind(int ncol, int rowm, int nnodes, int sampsz, int m_trt_switch,
             int obscorr_0, int obscorr_1, int *colm, int *cumcol, int *wt,
             double *fact, int *stpos, int *minvl, NODE *nodes) {
    int llim, hlim, k, j, npos, parm /*partial sum yi*/, llm, hlm, ill, ihh, count;
    int num, ipl, remcorr_0, remcorr_1, ipl1, nodind, high, up;
    double prarc, newprob;
    SUBND *cursnode;
    ARC *newarc, *carc, *narc;

    npos = stpos[ncol + 1];
    nodes[npos].subnodes = SNvector(0, 0);
    nodes[npos].subnodes[0].parcorr_0 = obscorr_0;
    nodes[npos].subnodes[0].parcorr_1 = obscorr_1;
    nodes[npos].subnodes[0].lpl = 0;
    nodes[npos].subnodes[0].spl = 0;
    nodes[npos].subnodes[0].numpred = 0; // number of predecessor nodes created
    nodes[npos].subnodes[0].tp = 0.0;
    nodes[npos].subnodes[0].arc = NULL;
    nodes[npos].subnodes[0].rec = NULL;

    double total = (ncol + 1) - (m_trt_switch + 2);
    RProgress::RProgress pb_for1(
        "Network calculations 1/2 [:bar] :current/:total (:percent)", total);
    pb_for1.tick(0);

    for (j = ncol + 1; j >= m_trt_switch + 2; j--) {
        llim = minvl[j];
        hlim = llim + stpos[j - 1] - stpos[j] - 1;
        for (parm = llim; parm <= hlim; parm++, npos++) {
            if (nodes[npos].subnodes == NULL) continue;
            int ref_1;
            llm = imax(0, parm - colm[j]);
            // from formula ( (k-1,u) in P(k, s_1k) ): u (the
            // minimum value u can take on for the k-1 stage)
            hlm = imin(parm, cumcol[j - 1]);
            // from formula ( (k-1,u) in P(k,
            // s_1k) ): u (the maximum value u
            // can take on for the k-1 stage)
            ill = stpos[j - 1] + llm - minvl[j - 1];
            ihh = ill + hlm - llm;
            count = parm - llm; // from Formula: s_1k - u
            up = imin(nodes[npos].upper, obscorr_0 + obscorr_1);
            for (nodind = nodes[npos].lower; nodind <= up; nodind++) {
                cursnode = &(nodes[npos].subnodes[nodind]);
                if ((cursnode->parcorr_1) < 0)
                    continue; // for terminal node, this will always be
                // equivalent to s2^1 > 0 or s2^2 > 0
                num = count; // s_1k - u
                // Create the Predecessor Nodes (ex: 5,3,2,2 and 5,4,2,3):
                for (k = ill; k <= ihh; k++, num--) { // k represents node position
                    ipl = num * (colm[j] - num); // last term of v in the Formula (part
                    // (a)): (s_1k - u)(nk - s_1k + u)

                    remcorr_0 = obscorr_0; // v in the Formula (part (a))
                    remcorr_1 = cursnode->parcorr_1 - ipl;
                    if (j == m_trt_switch + 2 && remcorr_1 != 0) continue;
                    if (j != m_trt_switch + 2 &&
                        (remcorr_1 > obscorr_1 || remcorr_1 < 0))
                        continue;
                    (cursnode->numpred)++; // number of predecessor nodes created
                    if (nodes[k].subnodes == NULL) { // initialize subnode
                        high = imin(nodes[k].upper, obscorr_0 + obscorr_1);
                        nodes[k].subnodes = SNvector(nodes[k].lower, high);
                        init(nodes[k].lower, high, nodes[k].subnodes);
                    }
                    newarc = (ARC *)malloc(sizeof(ARC));
                    ipl1 = wt[j] * num; // from Formula: x_k(s_1k - u) = RANK/ARC LENGTH
                    prarc = prbar0(colm[j], num, fact);
                    // colm[j] = column total (n_k), num = (s_1k -
                    // u), from Formula: choose(n_k,(s_1k - u))
                    // PROBABILITY LENGTH OF CONNECTING ARC
                    newarc->child = cursnode;
                    newarc->arc = ipl1; // RANK/ARC LENGTH
                    newarc->pr = prarc; // PROBABILITY LENGTH OF CONNECTING ARC
                    newarc->nextarc = NULL;
                    ref_1 = (obscorr_0 + obscorr_1) - (remcorr_0 + remcorr_1);
                    if ((nodes[k].subnodes[ref_1].lpl) < 0) {
                        nodes[k].subnodes[ref_1].parcorr_0 = remcorr_0;
                        nodes[k].subnodes[ref_1].parcorr_1 = remcorr_1;
                        nodes[k].subnodes[ref_1].lpl = (cursnode->lpl) + ipl1;
                        nodes[k].subnodes[ref_1].spl = (cursnode->spl) + ipl1;
                        nodes[k].subnodes[ref_1].tp = (cursnode->tp) + prarc;
                        // from Formula:  TP(k-1,u,v)=c_ok *
                        // TP(k,s_1k,s_2k) (we add instead of
                        // multiply because we are working with the
                        // LOG probabilities)
                        nodes[k].subnodes[ref_1].arc = newarc;
                    } else { // Else, (k-1,u,v) already exists (i.e., was found
                        // previously to be the predecessor of another
                        // node at the kth stage) and then:
                        if ((cursnode->lpl) + ipl1 > nodes[k].subnodes[ref_1].lpl)
                            // from Formula: max{LP_1(k-1,u,v), r_1k
                            // + LP_1(k,s_1k,s_2k)} (order switched)
                            nodes[k].subnodes[ref_1].lpl = cursnode->lpl + ipl1;
                        // from Formula:
                        // LP_1(k-1,u,v)=max{LP_1(k-1,u,v), r_1k
                        // + LP_1(k,s_1k,s_2k)}
                        if ((cursnode->spl) + ipl1 < nodes[k].subnodes[ref_1].spl)
                            // from Formula: min{SP_1(k-1,u,v), r_1k
                            // + SP_1(k,s_1k,s_2k)} (order switched)
                            nodes[k].subnodes[ref_1].spl = cursnode->spl + ipl1;
                        // from Formula:
                        // SP_1(k-1,u,v)=min{SP_1(k-1,u,v), r_1k
                        // + SP_1(k,s_1k,s_2k)}
                        newprob = prarc + cursnode->tp; // from Formula: c_ok *
                        // TP(k,s_1k,s_2k) -- used below
                        nodes[k].subnodes[ref_1].tp =
                            addlog(nodes[k].subnodes[ref_1].tp, newprob);
                        // from Formula : TP(k-1,u,v)=TP(k-1,u,v)
                        // + c_ok * TP(k,s_1k,s_2k)
                        if (nodes[k].subnodes[ref_1].arc == NULL)
                            nodes[k].subnodes[ref_1].arc = newarc;
                        else {
                            carc = nodes[k].subnodes[ref_1].arc;
                            narc = carc->nextarc;
                            while (narc != NULL) {
                                carc = carc->nextarc;
                                narc = carc->nextarc;
                            }
                            carc->nextarc = newarc;
                        }
                    }
                }
                // if node has no predeccesors then drop
                if ((cursnode->numpred) == 0) dropnd(cursnode);
            }
        }
        pb_for1.tick();
    }

    total = m_trt_switch - 1;
    RProgress::RProgress pb_for2(
        "Network calculations 2/2 [:bar] :current/:total (:percent)", total);
    pb_for2.tick(0);
    for (j = m_trt_switch + 1; j >= 2; j--) {
        llim = minvl[j];
        hlim = llim + stpos[j - 1] - stpos[j] - 1;
        for (parm = llim; parm <= hlim; parm++, npos++) {
            // go through the number of nodes (in the gamma1
            // network) at Stage j - 1
            // at the beginning, we create the terminal node, so
            // nodes[npos].subnodes.parcorr = s2, nodes[npos].subnodes.parcorr_0
            // = s2^1, nodes[npos].subnodes.parcorr_1 = s2^2
            if (nodes[npos].subnodes == NULL)
                continue; // The continue statement forces the next iteration
            // of the loop to take place, skipping any code in
            // between.
            // hit "continue" if there are no quadruples to be created from the
            // Gamma(s_1) network for the node (j-1, parm)
            // int up_1, up_2, high_1, high_2, ref_1, ref_2;
            int ref_1;
            // Grab the Successor Node (ex: 6,4,3,2):
            llm = imax(0, parm - colm[j]);
            // from formula ( (k-1,u) in P(k, s_1k) ): u (the
            // minimum value u can take on for the k-1 stage)
            hlm = imin(parm, cumcol[j - 1]); // from formula ( (k-1,u) in P(k,
            // s_1k) ): u (the maximum value u
            // can take on for the k-1 stage)
            ill = stpos[j - 1] + llm - minvl[j - 1];
            ihh = ill + hlm - llm;
            count = parm - llm; // from Formula: s_1k - u
            up = imin(nodes[npos].upper, obscorr_0 + obscorr_1);
            for (nodind = nodes[npos].lower; nodind <= up; nodind++) {
                cursnode = &(nodes[npos].subnodes[nodind]);
                if ((cursnode->parcorr_0) < 0)
                    continue; // for terminal node, this will always be
                // equivalent to s2^1 > 0 or s2^2 > 0
                num = count; // s_1k - u
                // Create the Predecessor Nodes (ex: 5,3,2,2 and
                // 5,4,2,3):
                for (k = ill; k <= ihh; k++, num--) { // k represents node position
                    ipl = num * (colm[j] - num); // last term of v in the Formula (part
                    // (a)): (s_1k - u)(nk - s_1k + u)

                    remcorr_0 = cursnode->parcorr_0 - ipl;
                    // v in the Formula (part (a))
                    remcorr_1 = 0;
                    if (j == 2 && remcorr_0 != 0) continue;
                    if (j != 2 && (remcorr_0 > obscorr_0 || remcorr_0 < 0)) continue;

                    // From formula: if ... then DROP NODE - DON'T INCLUDE IN
                    // NETWORK
                    (cursnode->numpred)++; // number of predecessor nodes created
                    if (nodes[k].subnodes == NULL) { // initialize subnode
                        high = imin(nodes[k].upper, obscorr_0 + obscorr_1);
                        nodes[k].subnodes = SNvector(nodes[k].lower, high);
                        init(nodes[k].lower, high, nodes[k].subnodes);
                    }
                    newarc = (ARC *)malloc(sizeof(ARC));
                    ipl1 = wt[j] * num; // from Formula: x_k(s_1k - u) = RANK/ARC LENGTH
                    prarc = prbar0(colm[j], num, fact);
                    // colm[j] = column total (n_k), num = (s_1k -
                    // u), from Formula: choose(n_k,(s_1k - u))
                    // PROBABILITY LENGTH OF CONNECTING ARC
                    newarc->child = cursnode;
                    newarc->arc = ipl1; // RANK/ARC LENGTH
                    newarc->pr = prarc; // PROBABILITY LENGTH OF CONNECTING ARC
                    newarc->nextarc = NULL;
                    ref_1 = (obscorr_0 + obscorr_1) - (remcorr_0 + remcorr_1);
                    if ((nodes[k].subnodes[ref_1].lpl) < 0) {
                        nodes[k].subnodes[ref_1].parcorr_0 = remcorr_0;
                        nodes[k].subnodes[ref_1].parcorr_1 = remcorr_1;
                        nodes[k].subnodes[ref_1].lpl = (cursnode->lpl) + ipl1;
                        // from Formula: LP_1(k-1,u,v)=r_1k
                        // + LP_1(k,s_1k,s_2k) (switch order
                        // of addition)
                        nodes[k].subnodes[ref_1].spl = (cursnode->spl) + ipl1;
                        // from Formala: SP_1(k-1,u,v)=r_1k
                        // + SP_1(k,s_1k,s_2k) (switch order
                        // of addition)
                        nodes[k].subnodes[ref_1].tp = (cursnode->tp) + prarc;
                        // from Formula:  TP(k-1,u,v)=c_ok *
                        // TP(k,s_1k,s_2k) (we add instead of
                        // multiply because we are working with the
                        // LOG probabilities)
                        nodes[k].subnodes[ref_1].arc = newarc;
                    } else { // Else, (k-1,u,v) already exists (i.e., was found
                        // previously to be the predecessor of another
                        // node at the kth stage) and then:
                        if ((cursnode->lpl) + ipl1 > nodes[k].subnodes[ref_1].lpl)
                            // from Formula: max{LP_1(k-1,u,v), r_1k
                            // + LP_1(k,s_1k,s_2k)} (order switched)
                            nodes[k].subnodes[ref_1].lpl = cursnode->lpl + ipl1;
                        // from Formula:
                        // LP_1(k-1,u,v)=max{LP_1(k-1,u,v), r_1k
                        // + LP_1(k,s_1k,s_2k)}
                        if ((cursnode->spl) + ipl1 < nodes[k].subnodes[ref_1].spl)
                            // from Formula: min{SP_1(k-1,u,v), r_1k
                            // + SP_1(k,s_1k,s_2k)} (order switched)
                            nodes[k].subnodes[ref_1].spl = cursnode->spl + ipl1;
                        // from Formula:
                        // SP_1(k-1,u,v)=min{SP_1(k-1,u,v), r_1k
                        // + SP_1(k,s_1k,s_2k)}
                        newprob = prarc + cursnode->tp; // from Formula: c_ok *
                        // TP(k,s_1k,s_2k) -- used below
                        nodes[k].subnodes[ref_1].tp =
                            addlog(nodes[k].subnodes[ref_1].tp, newprob);
                        // from Formula : TP(k-1,u,v)=TP(k-1,u,v)
                        // + c_ok * TP(k,s_1k,s_2k)
                        if (nodes[k].subnodes[ref_1].arc == NULL)
                            nodes[k].subnodes[ref_1].arc = newarc;
                        else {
                            carc = nodes[k].subnodes[ref_1].arc;
                            narc = carc->nextarc;
                            while (narc != NULL) {
                                carc = carc->nextarc;
                                narc = carc->nextarc;
                            }
                            carc->nextarc = newarc;
                        }
                    }
                }
                if ((cursnode->numpred) == 0)
                    // if node has no predeccesors (ex: if 5,3,2,2 was not
                    // connected to any nodes at Stage 4, then drop the
                    // 5,3,2,2 node) then drop
                    dropnd(cursnode);
            }
        }
        pb_for2.tick();
    }
}

void dropnd(SUBND *cursnode) {
    SUBND *cnode;
    ARC *carc, *narc;

    cursnode->lpl = -1;
    cursnode->spl = -1;
    cursnode->tp = -1;
    cursnode->parcorr_0 = -1;
    cursnode->parcorr_1 = -1;
    carc = cursnode->arc;
    while (carc != NULL) {
        cnode = carc->child;
        (cnode->numpred)--;
        if ((cnode->numpred) == 0) dropnd(cnode);
        narc = carc->nextarc;
        free(carc);
        carc = narc;
    }
    cursnode->arc = NULL;
    return;
}

void init(int lower, int upper, SUBND *subnodes) {
    int k;

    for (k = lower; k <= upper; k++) {
        subnodes[k].parcorr_0 = -1;
        subnodes[k].parcorr_1 = -1;
        subnodes[k].lpl = -1;
        subnodes[k].spl = -1;
        subnodes[k].tp = -1;
        subnodes[k].numpred = 0;
        subnodes[k].arc = NULL;
        subnodes[k].rec = NULL;
    }
    return;
}

SUBND *crenode(int remcorr_0, int remcorr_1, int lp, int sp, int ipl, double prarc) {
    SUBND *newnode;

    newnode = (SUBND *)malloc(sizeof(SUBND));
    newnode->parcorr_0 = remcorr_0;
    newnode->parcorr_1 = remcorr_1;
    newnode->lpl = lp + ipl;
    newnode->spl = sp + ipl;
    newnode->tp = prarc;
    newnode->numpred = 0;
    return (newnode);
}

void forind(int ncol, int rowm, int nnodes, int sampsz, int *colm, int *cumcol,
            double *fact, int *stpos, int *minvl, NODE *nodes) {
    int llim, hlim, k, j, npos, spl, lpl;
    // shortest path and longest path for the correlation sufficient statistic

    npos = stpos[1];
    nodes[npos].splcorr = 0;
    nodes[npos].lplcorr = 0;
    nodes[npos].subnodes = NULL;
    for (j = 2; j <= ncol + 1; j++) {
        npos = stpos[j];
        llim = minvl[j];                           // min value of s1_k at stage j
        hlim = llim + stpos[j - 1] - stpos[j] - 1; // max value of s1_k at stage j
        for (k = llim; k <= hlim; k++) {
            corrlpsp(sampsz, k, j, ncol, rowm, nnodes, colm, cumcol, stpos, minvl,
                     nodes, &spl, &lpl);
            nodes[npos].splcorr = spl;
            nodes[npos].lplcorr = lpl;
            nodes[npos].numsucc = 0;
            nodes[npos].subnodes = NULL;
            npos++;
        }
    }
}

void corrlpsp(int sampsz, int parm, int stage, int ncol, int rowm, int nnodes,
              int *colm, int *cumcol, int *stpos, int *minvl, NODE *nodes, int *spl,
              int *lpl) {
    int llm, hlm, ill, ihh, isp /*ith shortest path*/, ilp /*ith longest path*/, ipl,
        isp2, ilp2, count, pos;

    llm = imax(0, parm - colm[stage]);
    hlm = imin(parm, cumcol[stage - 1]);
    ill = stpos[stage - 1] + llm - minvl[stage - 1];
    // use arc lengths from nodes at this
    // position to nodes at ihh position
    ihh = ill + hlm - llm;

    count = parm - llm; // parm = value of k (s1_k at stage j);  this part of
    // r_2k formula: s_(1,k+1) - s_1k
    ipl = count * (colm[stage] - count);
    // r_2k = (s_(1,k+1)-s_1k )(n_(k+1)-s_(1,k+1)+s_1k )
    isp = nodes[ill].splcorr + ipl;
    // ith shortest path - starting with arc from ill
    // then will go to arc ihh in the below for loop
    ilp = nodes[ill].lplcorr + ipl;
    if (ill < ihh) {
        for (pos = ill + 1; pos <= ihh; pos++) {
            count--;
            ipl = count * (colm[stage] - count);
            isp2 = nodes[pos].splcorr + ipl;
            ilp2 = nodes[pos].lplcorr + ipl;
            if (isp2 < isp) isp = isp2;
            if (ilp2 > ilp) ilp = ilp2;
        }
    }
    *spl = isp;
    *lpl = ilp;
    return;
}

// calculates nodes
// builds network - before the first pass
void calnds(int ncol, int *table, int rowm, int *nnodes, int *colm, int *cumcol,
            int *stpos, int *minvl) {
    // Formula: R(k-1, s_(1,k-1)) = {(k,u):max(s_(1,k-1), s_1 - sum^N_(l=k)nl)
    // <= u <= min(s_1, s_(1,k-1) + n_k)} where rowm = s_1, cumcol[ncol + 1] =
    // sum^N_(l=k)nl table2 = sample size for each cluster, stored backwards
    // (ex: for 3 clusters, if colm = [0,2,3,3], then sample size for cluster 1
    // = 3, cluster 2 = 3, cluster 3 = 2) colm = sample size for each cluster,
    // stored backwards (ex: for 3 clusters, if colm = [0,2,3,3], then sample
    // size for cluster 1 = 3, cluster 2 = 3, cluster 3 = 2)

    int i, stage, iconst, llim, hlim, npos;
    colm[1] = 0;
    cumcol[1] = 0;

    // build colm and cumcol (cummulative sum of colm) vectors
    for (i = 2; i <= ncol + 1; i++) {
        colm[i] = table[i - 2]; // column total - sample size in each treatment/study
        cumcol[i] = cumcol[i - 1] + colm[i];
        // cumulative column totals - 7 is the total sample size
    }
    // this part of Formula: s_1 - sum^N_(l=k)nl FOR THE INITIAL NODE - at the
    // initial node, this number will be the most negative it can be, then it
    // decreases. That's why we add a positive number to this negative number
    // down in the for loop:
    iconst = rowm - cumcol[ncol + 1];
    npos = 2;
    minvl[ncol + 1] = rowm; // this part of Formula: s_1
    stpos[ncol + 1] = 1;
    // goes through stages 0 to ncol - 1 and finds the minimum partial sum
    // (llim) and maximum parital sum (hlim) for that stage if llim = 0 and hlim
    // = 2, then at that stage (say, stage k), there are 3 successor nodes
    // ([k,0], [k,1], and [k,2]) npos is a cumulative sum of the number of nodes
    // at each stage
    for (stage = ncol - 1; stage >= 0; stage--) {
        llim = imax(0, iconst + cumcol[stage + 1]);
        // lower limit for u -> this part of Formula:
        // max(s_(1,k-1), s_1 - sum^N_(l=k)nl) - see
        // comment above iconst declaration
        // essentially taking iconst + cummulative
        // sample size at each stage
        hlim = imin(rowm, cumcol[stage + 1]);
        // upper limit for u -> this part of Formula:
        // min(s_1, s_(1,k-1) + n_k) ... sort of(?)
        stpos[stage + 1] = npos;
        minvl[stage + 1] = llim;
        npos = npos + (hlim - llim + 1);
    }
    *nnodes = npos - 1;
    return;
}

double prbar0(int coltot, int numsucc, double *fact) {
    double logprob;

    if (coltot == numsucc || numsucc == 0)
        logprob = 0.0;
    else
        logprob = fact[coltot + 1] - fact[numsucc + 1] - fact[coltot - numsucc + 1];
    return (logprob);
}

double testmax(double x, double y) {
    if (x > y)
        return (x);
    else
        return (y);
}

double addlog(double num1, double num2) {
    double total2, t12, t22, tmax2;

    tmax2 = testmax(num1, num2);
    t12 = testmax(num1 - tmax2, -80.0);
    t22 = testmax(num2 - tmax2, -80.0);
    total2 = tmax2 + log(exp(t12) + exp(t22));
    return (total2);
}

void faclog(int sampsz, double *fact) {
    int i;

    fact[1] = 0.0;
    for (i = 1; i <= sampsz; i++) {
        fact[i + 1] = fact[i] + log((double)i);
    }
    return;
}

#define NR_END 1
#define FREE_ARG char *

NODE *NDvector(long nl, long nh)
/* allocate a NODE vector with subscript range v[nl..nh] */
{
    NODE *v;

    v = (NODE *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(NODE)));
    char message[] = "allocation failure in NDvector()";
    if (!v) nrerror(message);
    return v - nl + NR_END;
}

SUBND *SNvector(long nl, long nh)
/* allocate a NODE vector with subscript range v[nl..nh] */
{
    SUBND *v;

    v = (SUBND *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(SUBND)));
    char message[] = "allocation failure in SNvector()";
    if (!v) nrerror(message);
    return v - nl + NR_END;
}

void free_NDvector(NODE *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
    free((FREE_ARG)(v + nl - NR_END));
}

void free_SNvector(SUBND *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
    free((FREE_ARG)(v + nl - NR_END));
}

// [[Rcpp::export]]
DataFrame trstatWrapper(int ncol, IntegerVector table, IntegerVector uwt, int sampsz,
                        int rowm, int m_trt_switch, int obsstat, int obscorr_0,
                        int obscorr_1, double pval) {
    int *table_array = table.begin();
    int *uwt_array = uwt.begin();

    DataFrame distribution =
        trstat(&ncol, table_array, uwt_array, &sampsz, &rowm, &m_trt_switch, &obsstat,
               &obscorr_0, &obscorr_1, &pval);
    return distribution;
}
