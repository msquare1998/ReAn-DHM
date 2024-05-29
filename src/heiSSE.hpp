/* ============================================================================================
 *  Re-An method for measuring lnZ of the 1D dimerzied Heisenberg model
 *  Reference: https://arxiv.org/abs/2403.08642
 *  @Yiming_Ding, Westlake University
 *  Last updated: Mar 19, 2024
 * ============================================================================================*/
#include "heiSSE.h"

void heisenbergSSE::diagUpdate() {
    /* ---------------------------------
      opString[p] := 2 * b[p] + a[p]
    ----------------------------------*/
    int newBond, theBond, opstring_p;
    double prob_accept;

    for (int p = 0; p < M; ++p) {
        opstring_p = opString[p];

        if (opstring_p < 0) {
            newBond = getRandBond();
            if (spins[bSites[newBond][0]] != spins[bSites[newBond][1]]) {
                if ((prob_add_factor * Js[newBond] >= double(M - n)) or (prob_add_factor * Js[newBond] >= getRandProb() * double(M - n))) {
                    opString[p] = 2 * newBond;
                    n += 1;
                    if (isWeakBond[newBond]) {
                        nw += 1;
                    }
                }
            }
        }

        else if (opstring_p % 2 == 0) {
            theBond = opstring_p / 2;
            prob_accept = prob_remove_factor * double(M - n + 1) / Js[theBond];
            if ((prob_accept >= 1) or (getRandProb() <= prob_accept)) {
                opString[p] = -1;
                n -= 1;
                if (isWeakBond[theBond]) {
                    nw -= 1;
                }
            }
        }

        else {
            theBond = opstring_p / 2;
            spins[bSites[theBond][0]] *= -1;
            spins[bSites[theBond][1]] *= -1;
        }
    }
}

void heisenbergSSE::makeVertexList() {
    int opstring_p, b_p;
    int v_leg0;
    int s0, s1;
    int s0_vLast, s1_vLast;
    
    for (int v = 0; v < 4 * M; ++v)
        vertexList[v] = -1;
    for (int s = 0; s < nQ; ++s) {
        vFirst[s] = -1;
        vLast[s] = -1;
    }

    for (int p = 0; p < M; ++p) {
        opstring_p = opString[p];

        if (opstring_p != -1) {
            b_p = opstring_p / 2;
            s0 = bSites[b_p][0];
            s1 = bSites[b_p][1];
            v_leg0 = 4 * p;
            s0_vLast = vLast[s0];
            s1_vLast = vLast[s1];

            if (s0_vLast > -1) {
                vertexList[s0_vLast] = v_leg0;
                vertexList[v_leg0] = s0_vLast;
            } else
                vFirst[s0] = v_leg0;
            vLast[s0] = v_leg0 + 2;

            if (s1_vLast > -1) {
                vertexList[s1_vLast] = v_leg0 + 1;
                vertexList[v_leg0 + 1] = s1_vLast;
            } else
                vFirst[s1] = v_leg0 + 1;
            vLast[s1] = v_leg0 + 3;
        }
    }

    // PBC correction
    int s_vFirst;
    int s_vLast;
    for (int s = 0; s < nQ; ++s) {
        s_vFirst = vFirst[s];
        s_vLast = vLast[s];

        if (s_vFirst != -1) {
            vertexList[s_vFirst] = s_vLast;
            vertexList[s_vLast] = s_vFirst;
        }
    }
}

void heisenbergSSE::loopUpdate() {
    int v_head, v_tail;

    for (int v = 0; v < 4 * M; v += 2) {
        if (vertexList[v] < 0)
            continue;

        v_head = v;

        if (getRandProb() < 0.5f) {
            do {
                opString[ v_head / 4 ] ^= 1;
                vertexList[ v_head ] = -2;
                v_tail = v_head ^ 1;
                v_head = vertexList[v_tail];
                vertexList[v_tail] = -2;
            } while (v_head != v);
        }

        else {
            do {
                vertexList[v_head] = -1;
                v_tail = v_head ^ 1;
                v_head = vertexList[v_tail];
                vertexList[v_tail] = -1;
            } while (v_head != v);
        }
    }

    for (int i = 0; i < nQ; ++i) {
        if (vFirst[i] != -1) {
            if (vertexList[vFirst[i]] == -2)
                spins[i] *= -1;
        }

        else {
            if (getRandProb() < 0.5)
                spins[i] *= -1;
        }
    }
}

void heisenbergSSE::adjustM() {
    int newM = n + n / 3;
    if (M < newM){
        int *opString_copy = new int [M];
        for (int i = 0; i < M; ++i)
            opString_copy[i] = opString[i];
        delete []opString;
        opString = new int [newM];
        for (int i = 0; i < M; ++i)
            opString[i] = opString_copy[i];
        for (int i = M; i < newM; ++i)
            opString[i] = -1;
        M = newM;
        delete []vertexList;
        vertexList = new int[4 * M];
    }
}

void heisenbergSSE::updateConfig() {
    diagUpdate();
    makeVertexList();
    loopUpdate();
}

double heisenbergSSE::getRandProb() {
    return double((rand_r(&seed) % RAND_MAX)) / double(RAND_MAX);
}

int heisenbergSSE::getRandBond() {
    return (rand_r(&seed) % nBonds);
}

void heisenbergSSE::iniMeasure() {
    zzRatio = 0.0;
    nMeasure = 0;
}

void heisenbergSSE::measure() {
    zzRatio += pow(alpha_inv, nw);
    nMeasure += 1;
}

void heisenbergSSE::statisticize() {
    zzRatio /= double(nMeasure);
}