/* ============================================================================================
 *  Re-An method for measuring lnZ of the 1D dimerzied Heisenberg model
 *  Reference: https://arxiv.org/abs/2403.08642
 *  @Yiming_Ding, Westlake University
 *  Last updated: Mar 19, 2024
 * ============================================================================================*/
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>

class heisenbergSSE {
public:
    heisenbergSSE(int l, double a_inv, double bt, double jw, int threadId);
    ~heisenbergSSE();
    int L;                  // length of the chain which is assumed to be even in this program
    int nQ;                 // number of qubits (1/2 spins)
    int nBonds;             // number of bonds
    double beta;             // inverse of T
    int n;                  // number of non-identity operators
    int nw;                 // number of operators on the weak bond
    int M;                  // cut-off of SSE
    double Jw;              // ratio between J2 and J1
    unsigned int seed;      // seed for pseudo-random numbers

    int **bSites;           // for saving bonds and the sites they connect
    double *Js;                 // couplings on different bonds
    bool *isWeakBond;
    int *spins;             // for saving spin configuration
    int *opString;

    int *vertexList;
    int *vFirst;            // for each site, record the first site it links
    int *vLast;             // for each site, record the last site it links

    void diagUpdate();
    void makeVertexList();
    void loopUpdate();
    void adjustM();
    void updateConfig();

    inline void iniMeasure();      // initialize relevant quantities for a new measurement
    inline void measure();         // do measurements
    inline void statisticize();    // do calculations

    double zzRatio;
    double alpha_inv;

protected:
    int nMeasure;                               // times of measurements
    inline double getRandProb();
    inline int getRandBond();

    double prob_add_factor;                      // == 0.5 * beta * nBonds; will be divided by (M - n)
    double prob_remove_factor;                   // == 2.0 / (beta * nBonds); will multiply by (M - n + 1)
};

heisenbergSSE::heisenbergSSE(int l, double a_inv, double bt, double jw, int threadId) {
    seed = static_cast<unsigned int>(time(nullptr) + threadId);
    // =============================
    //  Assign the params
    // =============================
    this->L = l;
    this->nQ = l;
    this->nBonds = l;
    this->beta = bt;
    this->Jw = jw;
    this->alpha_inv = a_inv;

    this->n = 0;
    this->nw = 0;
    this->M = 20;

    // =============================
    //  Assign the two factors
    // =============================
    this->prob_add_factor = 0.5 * beta * double(nBonds);
    this->prob_remove_factor = 2.0 / ( beta * double(nBonds));

    // =============================
    //  Make bSites 1D
    // =============================
    Js = new double [nBonds];
    isWeakBond = new bool [nBonds];

    std::vector<std::vector<int>> bSites_vec;
    std::vector<int> bij(2);

    for (int i = 0; i < L; ++i) {
        bij.at(0) = i;
        bij.at(1) = (i + 1) % L;
        bSites_vec.push_back(bij);

        if (i % 2 == 0) {
            Js[i] = Jw;
            isWeakBond[i] = true;
        }
        else {
            Js[i] = 1.0;
            isWeakBond[i] = false;
        }
    }

    bSites = new int *[nBonds];
    for ( int b = 0; b < nBonds; ++b ) {
        bSites[b] = new int [2];
    }
    for ( int b = 0; b < nBonds; ++b ) {
        bSites[b][0] = bSites_vec.at(b).at(0);
        bSites[b][1] = bSites_vec.at(b).at(1);
    }

    // =============================
    //  Make initial spins
    // =============================
    spins = new int[nQ];
    for ( int i = 0; i < nQ; ++i ) {
        getRandProb() > 0.5? spins[i] = 1: spins[i] = -1;
    }

    // ===============================================
    //  Fill 'opString' with -1 elements (identity)
    // ===============================================
    opString = new int [ M ];
    for ( int i = 0; i < M; ++i )
        opString[i] = -1;

    // =============================
    //  Allocate vFirst and vLast
    // =============================
    vFirst = new int  [nQ];
    vLast = new int [nQ];
    vertexList = new int[4 * M];
}

heisenbergSSE::~heisenbergSSE() {
    delete[] bSites;
    delete[] spins;
    delete[] opString;
    delete[] vFirst;
    delete[] vLast;
    delete[] Js;
    delete[] vertexList;
    delete[] isWeakBond;
}

void timer(clock_t t0, clock_t t1) {
    int interval = (int) (t1 - t0);
    auto hour = interval / 3600;
    auto min = (interval - 3600 * hour) / 60;
    auto sec = interval - 3600 * hour - 60 * min;
    std::cout << "◆ Time used: " << hour << " hour " << min << " min " << sec << " sec " << std::endl;
}
