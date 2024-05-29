/* ============================================================================================
 *  Re-An method for measuring lnZ of the 1D dimerzied Heisenberg model
 *  Reference: https://arxiv.org/abs/2403.08642
 *  @Yiming_Ding, Westlake University
 *  Last updated: Mar 19, 2024
 * ============================================================================================*/
#include "heiSSE.hpp"
#include <pthread.h>

static int L, stepThm, stepStat, nDivision;
static double BETA, Jw, Jw0, Lambda, epsilon;
static double *AlphaList, *JwList;
static double *Ratios;

void doDivision() {
    std::vector<double> AlphaList_vec, JwList_vec;
    double Jw_copy = Jw;
    double alpha;
    
    while (Jw_copy > Jw0) {
        JwList_vec.push_back(Jw_copy);

        alpha = pow(epsilon, 1.0 / (Lambda * BETA * Jw_copy * L));      // For 1D system, estimate n ~ ΛβJL
        if (alpha < epsilon)
            alpha = epsilon;
        else if (alpha > 0.99999)
            alpha = 0.99999;
        AlphaList_vec.push_back(alpha);
        Jw_copy *= alpha;
    }
    nDivision = static_cast<int>(JwList_vec.size());
    AlphaList = new double [nDivision];
    JwList = new double [nDivision];
    for (int i = 0; i < nDivision; ++i) {
        AlphaList[i] = AlphaList_vec.at(i);
        JwList[i] = JwList_vec.at(i);
    }
}

void* threadFunc(void *arg) {
    auto *idx = static_cast<int*>(arg);

    auto *model = new heisenbergSSE(L, AlphaList[*idx], BETA, JwList[*idx], *idx);
    for (int i = 0; i < stepThm; ++i) {
        model -> updateConfig();
        model -> adjustM();
    }
    model -> iniMeasure();
    for (int i = 0; i < stepStat; ++i) {
        model -> updateConfig();
        model -> measure();
    }
    model -> statisticize();
    Ratios[*idx] = model -> zzRatio;
    delete model;

    return nullptr;
}

int main( int argc, char **argv ) {
    clock_t t0 = time( nullptr );
    L = std::stoi(argv[1]);
    epsilon = std::stod(argv[2]);
    BETA = std::stod(argv[3]);
    Lambda = std::stod(argv[4]);
    stepThm = std::stoi(argv[5]);
    stepStat = std::stoi(argv[6]);
    Jw0 = std::stod(argv[7]);
    Jw = std::stod(argv[8]);
    int nThread = std::stoi(argv[9]);
    int nBins = std::stoi(argv[10]);

    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //  Initialization A
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    doDivision();
    int cycle = nDivision / nThread;
    int remain = nDivision - nThread * cycle;

    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //  Report the environment
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    std::cout << "@ L = " << L << " Dimerzied Heisenberg model, beta = " << BETA << ", nBins = " << nBins << std::endl;
    std::cout << "  # number of division = " << nDivision << std::endl;
    std::cout << "  # number of threads = " << nThread << ", cycle = " << cycle << ", remain = " << remain << std::endl;

    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //  Initialization B
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    int nBonds_Js = L / 2;
    Ratios = new double [nDivision];
    auto *ps = new pthread_t [nDivision];
    auto *idx = new int [nDivision];
    for (int i = 0; i < nDivision; ++i)
        idx[i] = i;

    std::string file = "../../data/L" + std::to_string(L) + "_epslion" + argv[2] +
                           "_beta" + argv[3] + "_JwRef" + argv[7] + "_Jw" + argv[8] + "/JwList.dat";
    std::fstream f0; f0.precision(16);

    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    //  QMC simulation
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    std::string file_lnZ = "../../data/L" + std::to_string(L) + "_epslion" + argv[2] +
                           "_beta" + argv[3] + "_JwRef" + argv[7] + "_Jw" + argv[8] + "/lnZ.dat";
    f0.open(file_lnZ, std::ios::app);

    for (int bin = 0; bin < nBins; ++bin) {
        std::cout << "  # Start running bin " << bin << std::endl;
        // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //  Parallelization
        // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        for (int c = 0; c < cycle; ++c) {
	    std::cout << "\t--> Cycle " << c << ": [" << c * nThread << ", " << (c + 1) * nThread << ")\n";
            for (int p = 0; p < nThread; ++p) {
                pthread_create(&ps[c * nThread + p], nullptr, threadFunc, &idx[c * nThread + p]);
            }
            for (int p = 0; p < nThread; ++p) {
                pthread_join(ps[c * nThread +p], nullptr);
            }
        }

        if (remain > 0) {
            std::cout << "\t~~> Extra cycle: [" << cycle * nThread << ", " << nDivision << ")\n";
            for (int p = cycle * nThread; p < nDivision; ++p) {
		        pthread_create(&ps[p], nullptr, threadFunc, &idx[p]);
            }
            for (int p = cycle * nThread; p < nDivision; ++p) {
                pthread_join(ps[p], nullptr);
            }
        }

        // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //  Calculate lnPT
        // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        double lnzLnzSubtrc = 0.0;
        double zzBatchRatio = 1.0;
        for (int i = 0; i < nDivision; ++i) {
            zzBatchRatio *= Ratios[i];
            if (zzBatchRatio < 1e-10) {
                lnzLnzSubtrc += log(zzBatchRatio);
                zzBatchRatio = 1.0;
            }
        }
        lnzLnzSubtrc += log(zzBatchRatio);
        f0 << -BETA * double(nBonds_Js) * (-1.0) - lnzLnzSubtrc << std::endl;
    }

    f0.close();

    // =====================================
    //  Deallocate memory
    // =====================================
    delete[] AlphaList;
    delete[] JwList;
    delete[] Ratios;
    
    clock_t t1 = time( nullptr );
    for ( int i = 0; i < 77; ++i ) std::cout << "="; std::cout << std::endl;
    timer(t0, t1);
    return 0;
}