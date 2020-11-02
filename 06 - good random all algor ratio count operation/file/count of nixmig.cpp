#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <set>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include "MersenneTwister.h"
using namespace std;

const bool DEBUG = false;

vector<int> L_nm; // tank in algorithm is L
vector<int> tmp; // tmpL in algorithm

int cost_nm;        // use count cost of all operation in algorithm
int moveCnt_nm;     // use count only move operation
int addCnt_nm;      // use count only add operation
int deleteCnt_nm;   // use count only delete operation

const int moveCst_nm   = 1;      // move cost
const int addCst_nm    = 1;      // add cost
const int deleteCst_nm = 1;      // delete cost

const double alpha = 5.464;
const double beta = 4.1;

void setCnt_nm(bool all = 1, bool move = 1, bool add = 1, bool deleteC = 1) {
    if (all == true) cost_nm = 0;
    if (move == true) moveCnt_nm = 0;
    if (add == true) addCnt_nm = 0;
    if (deleteC == true) deleteCnt_nm = 0;
}

const int ttl   = 5;
const int thres = 60;
const int overThres = 400;
const double a_const = 0.5;

int exNodes = 0;
int rNodes = 0;
int lc = 0;
int m; // min index

int getIdxMax(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_nm[i] > L_nm[idx])
            idx = i;
    }
    return idx;
}

int getIdxMin(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_nm[i] < L_nm[idx])
            idx = i;
    }
    return idx;
}

int getMax(int n) {
    return L_nm[getIdxMax(n)];
}

int getMin(int n) {
    return L_nm[getIdxMin(n)];
}

int md(int n, int m) { // mod function
    int rlt = n%m;
    if (n < 0) rlt += m;
    return rlt;
}

double getRatio(int n) {
    double mn = (double)getMin(n);
    double mx = (double)getMax(n);
    if (abs((int)mn) <= 0.0000001) return mx;
    return mx/mn;
}

void print(int no, int n) {
    char str[100];
    sprintf(str, "%d.txt", no);
    ofstream out(str);
    for (int i  = 0; i < n; ++i) {
        out << i << ' ' << L_nm[i] << '\n';
    }
    out.close();
    return;
}

void localWave(int p, int hop, int n) {
    lc = 0;

    tmp.resize(L_nm.size());
    copy(L_nm.begin(), L_nm.end(), tmp.begin());

    exNodes = 0;

    while (lc <= hop and exNodes <= hop) {
        int moveLoad = 0;
        int pi = md(p+lc, n);

        if (tmp[pi] > overThres) {
            moveLoad = a_const*(tmp[pi]-thres);
            moveLoad = (int)moveLoad;
        } else if (tmp[pi] > 0 and tmp[pi] > thres) {
            int plc = md(p + lc, n); // index in p+lc
            moveLoad = tmp[plc] - thres;
        }
        tmp[pi] = L_nm[pi] - moveLoad;

        int pii = md(pi+1, n);
        tmp[pii] = L_nm[pii] + moveLoad;

        exNodes = (int)(tmp[pii]/thres - 1);

        exNodes = exNodes < 0 ? 0 : exNodes;

        lc = lc + 1;
    }
    return ;
}

void remoteWave(int p, int exNode, int n) {
    int p_lc_1 = md(p, n); // index in p + lc + 1
    m = getIdxMin(n);

    int tmpLm = L_nm[m];
    int j = 0;

    while (tmpLm <= tmp[p_lc_1] and j <= exNode) {
        int m_j_1 = md(m+j+1, n);

        tmpLm += L_nm[m_j_1];

        j = j + 1;
    }
    rNodes = j;
    return;
}

void nix(int p, int load, int n) {
    int from = md(p, n);      //index of p
    int to   = md(from+1, n); //index of p+1

    if ( L_nm[from] >= load ) {
        L_nm[from] -= load;
        L_nm[to]   += load;

        moveCnt_nm += load * moveCst_nm;                    // count opeation move
    }
    return;
}

void mig(int far_node, int p, int load, int n) {
    // ------- nix (m->m-1, L[m]) ------------
    int neighbor_idx = md(far_node-1, n);
    int tmp_fn = L_nm[far_node];

    L_nm.erase(L_nm.begin()+far_node);                //   <----------------------------------o
                                                    //                                      |
    if (far_node < neighbor_idx) {                  //                                      |
        neighbor_idx -= 1;                          //                                      |
        neighbor_idx = md(neighbor_idx, n);         //                                      +--o
    }                                               //                                      |  |
                                                    //                                      |  V
    L_nm.insert(L_nm.begin()+neighbor_idx, tmp_fn);   // <-----------------------------------<O
    moveCnt_nm += tmp_fn *(moveCst_nm);                   //    count operation (delete + insert) = Move
    
    far_node = neighbor_idx;
    nix(far_node, L_nm[far_node], n);                    // nix(m -> m-1, L[m])

    L_nm.erase(L_nm.begin()+far_node);
    if (far_node < p) {
        p -= 1;
        p = md(p, n);
    }

    // -------  nix(p -> p+1 # p+1 is m that empty, load) -------------
    int next_p = md(p+1, n);
    L_nm.insert(L_nm.begin()+next_p, 0);  //After process m, value in m is empty
    nix(p, load, n);              // nix(p -> p+1 # p+1 is m that empty, load)
    return;
}

void nixmig(int p,int hop, int n) {
    localWave(p, ttl, n);
    if (exNodes > 0) {
        remoteWave(p+lc+1, exNodes, n);
    }

    for (int i = 0; i <= lc; ++i) {
        int pi = md(p+i, n); // index at p + i
        int load = 0;

        if (L_nm[pi] > overThres) {
            load = a_const * (L_nm[pi]-thres);
            load = (int)load;
        } else if (L_nm[pi] > 0 and L_nm[pi] > thres) {
            load = L_nm[pi] - thres;
        }

        nix(pi, load, n);
    }

    if (rNodes > 0) {
        if(DEBUG)printf("rNodes = %d\n", rNodes);
        for (int i = 0; i <= rNodes; ++i) {
            int plc = md(p+lc, n);
            int value = (int)(tmp[plc]/rNodes);

            int m_i_1 = md(m+i+1, n);
            mig(m_i_1, plc, value, n); // m is minimum index
        }
    }

    return;
}

void insertLoad(int idx, int& n) {
    L_nm[idx] += 1;
    addCnt_nm += addCst_nm;                                // count opertaion add

    nixmig(idx, ttl, n);
}

int rt(int n) { // random [0, n-1]
    MTRand m;
    int idx = m.randInt(n-1);
    if (idx <= n-1 && 0 <= idx) return idx;
    else {
        printf("Error index\n");
        assert(idx <= n-1 && 0 <= idx);
    }
    return -1;
}

int main(){
    srand(time(0));
    int T = 1;
    while (T--) {
        int n;
        n = 4;
        //scanf("%d", &n);

        L_nm.clear();

        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            L_nm.push_back(0);
        }

        /*

        // input data in each tank
        for (int i = 0; i < n; ++i) {
            int slt; // select tank
            scanf("%d", &slt);
            insertLoad(slt, n);
        }

        */
        char nameFileRatio[100];
        char nameFileCount[100];

        sprintf(nameFileRatio, "ratio(%d).txt", n);
        sprintf(nameFileCount, "count(%d).txt", n);
        ofstream outRatio(nameFileRatio);
        ofstream outCount(nameFileCount);

        outRatio << "0 " << getRatio(n) << '\n';

        setCnt_nm();
        long long total_nm = 0LL;
        // input data in max load in tank
        for (int i = 0; i < 100; ++i) {
            setCnt_nm(false);
            int idx = rt(n); // random [0, n-1]
            insertLoad(idx, n);

            outRatio << i+1 << ' ' << getRatio(n) << '\n';

            cost_nm = moveCnt_nm + addCnt_nm + deleteCnt_nm;
            outCount << i+1 << ' ' << cost_nm << '\n';

            total_nm += cost_nm;
            //print(i, n);
        }
        outCount << "total_nm = " << total_nm << '\n';
        outCount.close();
        outRatio.close();
    }
    printf("End process\n");
    return 0;
}