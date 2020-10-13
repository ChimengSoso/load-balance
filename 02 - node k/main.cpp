#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include "MersenneTwister.h"
using namespace std;

const bool DEBUG = false;

vector<int> tnk; // tank in algorithm is L
vector<int> tmp; // tmpL in algorithm

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
        if (tnk[i] > tnk[idx])
            idx = i;
    }
    return idx;
}

int getIdxMin(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (tnk[i] < tnk[idx])
            idx = i;
    }
    return idx;
}

int getMax(int n) {
    return tnk[getIdxMax(n)];
}

int getMin(int n) {
    return tnk[getIdxMin(n)];
}

int md(int n, int m) { // mod function
    int rlt = n%m;
    if (n < 0) rlt += m;
    return rlt;
}

void localWave(int p, int hop, int n) {
    lc = 0;

    tmp.resize(tnk.size());
    copy(tnk.begin(), tnk.end(), tmp.begin());

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
        tmp[pi] = tnk[pi] - moveLoad;

        int pii = md(pi+1, n);
        tmp[pii] = tnk[pii] + moveLoad;

        exNodes = (int)(tmp[pii]/thres - 1);

        exNodes = exNodes < 0 ? 0 : exNodes;

        lc = lc + 1;
    }

    if (DEBUG) printf("local Wave pass\n");
    return ;
}

void remoteWave(int p, int exNode, int n) {
    int p_lc_1 = md(p, n); // index in p + lc + 1
    m = getIdxMin(n);

    int tmpLm = tnk[m];
    int j = 0;

    while (tmpLm <= tmp[p_lc_1] and j <= exNode) {
        int m_j_1 = md(m+j+1, n);

        tmpLm += tnk[m_j_1];

        j = j + 1;
    }
    rNodes = j;

    if (DEBUG) printf("remote Wave pass\n");
    return;
}

void nix(int p, int load, int n) {
    int from = md(p, n);      //index of p
    int to   = md(from+1, n); //index of p+1

    if ( tnk[from] >= load ) {
        tnk[from] -= load;
        tnk[to]   += load;

        if (DEBUG) printf("do in if : ");
    }

    if (DEBUG) printf("nix pass\n");
    return;
}

void mig(int far_node, int p, int load, int n) {
    // ------- nix (m->m-1, L[m]) ------------
    int neighbor_idx = md(far_node-1, n);
    int tmp_fn = tnk[far_node];

    tnk.erase(tnk.begin()+far_node);
    if (far_node < neighbor_idx) {
        neighbor_idx -= 1;
        neighbor_idx = md(neighbor_idx, n);
    }

    tnk.insert(tnk.begin()+neighbor_idx, tmp_fn);
    far_node = neighbor_idx;
    nix(far_node, tnk[far_node], n);                    // nix(m -> m-1, L[m])

    tnk.erase(tnk.begin()+far_node);
    if (far_node < p) {
        p -= 1;
        p = md(p, n);
    }

    // -------  nix(p -> p+1 # p+1 is m that empty, load) -------------
    int next_p = md(p+1, n);
    tnk.insert(tnk.begin()+next_p, 0);  //After process m, value in m is empty
    nix(p, load, n);              // nix(p -> p+1 # p+1 is m that empty, load)

    if (DEBUG) printf("mig pass\n");
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

        if (tnk[pi] > overThres) {
            load = a_const * (tnk[pi]-thres);
            load = (int)load;
        } else if (tnk[pi] > 0 and tnk[pi] > thres) {
            load = tnk[pi] - thres;
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

            if(DEBUG)printf("i = %d, m+i+1 = %d, p+lc = %d, value = %d\n", i, md(m+i+1, n), md(p+lc, n), value);
        }
    }

    if (DEBUG) printf("nixmig pass\n");

    return;
}

void insertLoad(int idx, int& n) {
    tnk[idx] += 1;
    nixmig(idx, ttl, n);
    if (DEBUG) printf("insert load pass\n");
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
        out << i << ' ' << tnk[i] << '\n';
    }
    out.close();
    return;
}

int rt(int n) { // random [0, n-1]
    //return 0;
    MTRand m;
    return m.randInt(n-1);
}

int main(){
    srand(time(0));
    int T = 1;
    while (T--) {
        int n;
        n = 1024;
        //scanf("%d", &n);

        tnk.clear();
        
        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            tnk.push_back(0);
        }

        /*

        // input data in each tank
        for (int i = 0; i < n; ++i) {
            int slt; // select tank
            scanf("%d", &slt);
            insertLoad(slt, n);
        }

        */

        int constPoint = rt(n);

        char nameFileRatio[100];

        sprintf(nameFileRatio, "ratio(%d).txt", n);
        ofstream outRatio(nameFileRatio);

        outRatio << "0 " << getRatio(n) << '\n';

        // input data in max load in tank
        for (int i = 0; i < 100000; ++i) {
            int idx = constPoint;
            insertLoad(idx, n);
            outRatio << i+1 << ' ' << getRatio(n) << '\n';

            print(i, n);
        }
        outRatio.close();
    }
    return 0;
}
