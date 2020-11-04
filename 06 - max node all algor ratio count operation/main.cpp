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
using ll = long long;

vector<int> L_jc; // tank of JS => L_jc

vector<int> L_kg; // tank of Karger => L_kg

vector<int> L_nm_1; // tank of nixmig version 1 => L_nm_1
vector<int> tmp_nm_1;  // tmp_L in algorithm nixmig

vector<int> L_nm_2; // tank of nixmig version 2 => L_nm_2
vector<int> tmp_nm_2;  // tmp_L in algorithm nixmig
//
//
//
/** ============== Constant For JC's agorithm ============== */
const double alpha = 5.464;
const double beta = 4.1;
const double eps = 0.023;

int cost_jc;        // use count cost of all operation in algorithm
int moveCnt_jc;     // use count only move operation
int addCnt_jc;      // use count only add operation
int deleteCnt_jc;   // use count only delete operation

const int moveCst_jc   = 1;      // move cost
const int addCst_jc    = 1;      // add cost
const int deleteCst_jc = 1;      // delete cost

// =============================================================
//
//
//
//
//
/** ============= Constant For Karger's algorithm ============ */
int cost_kg;        // use count cost of all operation in algorithm
int moveCnt_kg;     // use count only move operation
int addCnt_kg;      // use count only add operation
int deleteCnt_kg;   // use count only delete operation

const int moveCst_kg   = 1;      // move cost
const int addCst_kg    = 1;      // add cost
const int deleteCst_kg = 1;      // delete cost
// =============================================================
//
//
//
//
//
/** ===== Constant For NigMix's agorithm  [Version 1] ====== */
const int ttl_nm_1   = 5;
const int thres_nm_1 = 60;
const int overthres_nm_1 = 400;
const double a_const_nm_1 = 0.5;

int exNodes_nm_1 = 0;
int rNodes_nm_1 = 0;
int lc_nm_1 = 0;
int m_nm_1; // min index of nixmig algor version 1

int cost_nm_1;        // use count cost of all operation in algorithm
int moveCnt_nm_1;     // use count only move operation
int addCnt_nm_1;      // use count only add operation
int deleteCnt_nm_1;   // use count only delete operation

const int moveCst_nm_1   = 1;      // move cost
const int addCst_nm_1    = 1;      // add cost
const int deleteCst_nm_1 = 1;      // delete cost
// =============================================================
//
//
//
//
//
/** ===== Constant For NigMix's agorithm  [Version 2] ====== */
const int ttl_nm_2   = 5;
const int thres_nm_2 = 60;
const int overthres_nm_2 = 400;
const double a_const_nm_2 = 0.5;

int exNodes_nm_2 = 0;
int rNodes_nm_2 = 0;
int lc_nm_2 = 0;
int m_nm_2; // min index of nixmig algor version 1

int cost_nm_2;        // use count cost of all operation in algorithm
int moveCnt_nm_2;     // use count only move operation
int addCnt_nm_2;      // use count only add operation
int deleteCnt_nm_2;   // use count only delete operation

const int moveCst_nm_2   = 1;      // move cost
const int addCst_nm_2    = 1;      // add cost
const int deleteCst_nm_2 = 1;      // delete cost
// =============================================================
//
//
//
//
//
/* = = ========== Fucntion get Idx Max ========== = = */
int getIdxMax_jc(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_jc[i] > L_jc[idx])
            idx = i;
    }
    return idx;
}

int getIdxMax_kg(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_kg[i] > L_kg[idx])
            idx = i;
    }
    return idx;
}

int getIdxMax_nm_1(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_nm_1[i] > L_nm_1[idx])
            idx = i;
    }
    return idx;
}

int getIdxMax_nm_2(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_nm_2[i] > L_nm_2[idx])
            idx = i;
    }
    return idx;
}
// =====================================================

/* = = ========== Fucntion get Idx Min ========== = = */
int getIdxMin_jc(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_jc[i] < L_jc[idx])
            idx = i;
    }
    return idx;
}

int getIdxMin_kg(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_kg[i] < L_kg[idx])
            idx = i;
    }
    return idx;
}

int getIdxMin_nm_1(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_nm_1[i] < L_nm_1[idx])
            idx = i;
    }
    return idx;
}

int getIdxMin_nm_2(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_nm_2[i] < L_nm_2[idx])
            idx = i;
    }
    return idx;
}
// =====================================================

/* = = ========= Fucntion get Value Max ========= = = */
int getMax_jc(int n) {
    return L_jc[getIdxMax_jc(n)];
}

int getMax_kg(int n) {
    return L_kg[getIdxMax_kg(n)];
}

int getMax_nm_1(int n) {
    return L_nm_1[getIdxMax_nm_1(n)];
}

int getMax_nm_2(int n) {
    return L_nm_2[getIdxMax_nm_2(n)];
}
// =====================================================

/* = = ========= Fucntion get Value Min ========= = = */
int getMin_jc(int n) {
    return L_jc[getIdxMin_jc(n)];
}

int getMin_kg(int n) {
    return L_kg[getIdxMin_kg(n)];
}

int getMin_nm_1(int n) {
    return L_nm_1[getIdxMin_nm_1(n)];
}

int getMin_nm_2(int n) {
    return L_nm_2[getIdxMin_nm_2(n)];
}
// =====================================================

/* = = ========== Fucntion MOD n by m =========== = = */
int md(int n, int m) { 
    int rlt = n%m;
    if (n < 0) rlt += m;
    return rlt;
}
// =====================================================

/* = = ========== Fucntion Good Random ========== = = */
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
// =====================================================

/* ############## Fucntion of Algor'js ############## */

void setCnt_jc(bool all = 1, bool move = 1, bool add = 1, bool deleteC = 1) {
    if (all == true) cost_jc = 0;
    if (move == true) moveCnt_jc = 0;
    if (add == true) addCnt_jc = 0;
    if (deleteC == true) deleteCnt_jc = 0;
}

void minbalance(int u, int& n) {
    int v = getIdxMin_jc(n);
    if (L_jc[u] > alpha*L_jc[v]) {
        int l = md(v-1, n);
        int r = md(v+1, n);
        int z = L_jc[l] < L_jc[r] ? l : r;

        L_jc[z] += L_jc[v];
        moveCnt_jc += L_jc[v]*moveCst_jc;      // count operation move

        #define start L_jc.begin()

        L_jc.erase(start+v);

        if (v < u)
            --u;

        int remain_up = (L_jc[u]+1)/2;
        int remain_down =   L_jc[u]/2;
        L_jc[u] = remain_up;

        moveCnt_jc += remain_up*moveCst_jc;    // count operation move
        L_jc.insert(start+u, remain_down);
    }
    return;
}

void split(int u, int n) {
    int w = getIdxMax_jc(n);
    if (L_jc[u] < L_jc[w]/beta) {
        int l = md(u-1, n);
        int r = md(u+1, n);
        int z = L_jc[l] < L_jc[r] ? l : r;

        if (L_jc[z] <= 2*L_jc[w]/beta) {

            L_jc[z] += L_jc[u];
            L_jc.erase(L_jc.begin() + u);

            if (u < w)
                --w;

            int remain_up = (L_jc[w]+1)/2;
            int remain_down = L_jc[w] / 2;

            L_jc[w] = remain_up;
            L_jc.insert(L_jc.begin() + w, remain_down);
        }else{
            int sum = L_jc[u] + L_jc[z];

            L_jc[u] = (sum+1)/2;
            L_jc[z] = sum /2;
        }
    }
}

double getRatio_jc(int n) {
    double mn = (double)getMin_jc(n);
    double mx = (double)getMax_jc(n);
    if (mn <= 0.0000001) return mx;
    return mx/mn;
}

void print_jc(int no, int n) {
    char str[100];
    sprintf(str, "js_%d.txt", no);
    ofstream out(str);
    for (int i  = 0; i < n; ++i) {
        out << i << ' ' << L_jc[i] << '\n';
    }
    out.close();
    return;
}

int query_jc(int s, int e, int n) {
    int tot = 0;        // tot = total ; it same last pos of interval
    int rlt = 0;        // num of node about [s, e]
    for (int i = 0; i < n; ++i) {
        tot += L_jc[i];
        if (tot > e) {  // last interval
            if (tot-L_jc[i]+1 <= e) ++rlt;
            break;
        }

        // live in Interval [s, e]
        if (s <= tot and tot <= e) {
            ++rlt;
        }
    }
    return rlt;
}

void insertLoad_jc(int idx, int& n) {
    L_jc[idx] += 1;
    addCnt_jc += addCst_jc;
    minbalance(idx, n);
}
// =====================================================

/* ########### Fucntion of Algor'karger ############# */

void setCnt_kg(bool all = 1, bool move = 1, bool add = 1, bool deleteC = 1) {
    if (all == true) cost_kg = 0;
    if (move == true) moveCnt_kg = 0;
    if (add == true) addCnt_kg = 0;
    if (deleteC == true) deleteCnt_kg = 0;
}

void divideEq(int a, int b) {
    moveCnt_kg += ( abs(L_kg[a]-L_kg[b])/2 ) *moveCst_kg;             // count operation move
    int n = (L_kg[a]+L_kg[b]+1)/2;
    int m = (L_kg[a]+L_kg[b])/2;
    L_kg.erase(L_kg.begin()+a);
    if (a < b)
        --b;
    L_kg.erase(L_kg.begin()+b);
    if ((int)L_kg.size() < a)
        a = max(0, (int)L_kg.size()-1);
    L_kg.insert(L_kg.begin()+a, m);
    L_kg.insert(L_kg.begin()+a, n);
    return;
}

void itemBalance(int i, int& n) {
    int j;

    do {
        j = rt(n);
    } while (j == i);


    if (L_kg[i] <= eps*L_kg[j] || L_kg[j] <= eps*L_kg[i]) {
        if (i == j+1) {
            divideEq(i, j);
        }else{ // case 2: i != j+1
            if (L_kg[j+1] > L_kg[i]) {
                divideEq(i, j);
            }else{
                int j_right =  md(j+1, n);
                
                moveCnt_kg += L_kg[j] *moveCst_kg;                // count operation move
                
                L_kg[j_right] += L_kg[j];
                L_kg[j] = 0;

                int X = (L_kg[i]+1)/2;
                int Y = (L_kg[i])/2;

                L_kg.erase(L_kg.begin()+j);

                if (j < i)
                    --i;

                moveCnt_kg += Y * moveCst_kg;                    // count operation add
                L_kg.erase(L_kg.begin()+i);
                L_kg.insert(L_kg.begin()+i, Y);
                L_kg.insert(L_kg.begin()+i, X);
            }
        }
    }
    return;
}

double getRatio_kg(int n) {
    double mn = (double)getMin_kg(n);
    double mx = (double)getMax_kg(n);
    if (mn <= 0.0000001) return mx;
    return mx/mn;
}

void print_kg(int no, int n) {
    char str[100];
    sprintf(str, "kg_%d.txt", no);
    ofstream out(str);
    for (int i  = 0; i < n; ++i) {
        out << i << ' ' << L_kg[i] << '\n';
    }
    out.close();
    return;
}

int query_kg(int s, int e, int n) {
    int tot = 0;        // tot = total ; it same last pos of interval
    int rlt = 0;        // num of node about [s, e]
    for (int i = 0; i < n; ++i) {
        tot += L_kg[i];
        if (tot > e) {  // last interval
            if (tot-L_kg[i]+1 <= e) ++rlt;
            break;
        }

        // live in Interval [s, e]
        if (s <= tot and tot <= e) {
            ++rlt;
        }
    }
    return rlt;
}

void insertLoad_kg(int idx, int& n) {
    L_kg[idx] += 1;
    addCnt_kg += addCst_kg;                                   // count operation add
    itemBalance(idx, n);
}
// =====================================================

/* ########### Fucntion of Algor'Nixmig VERSION 1 ############# */

void setCnt_nm_1(bool all = 1, bool move = 1, bool add = 1, bool deleteC = 1) {
    if (all == true) cost_nm_1 = 0;
    if (move == true) moveCnt_nm_1 = 0;
    if (add == true) addCnt_nm_1 = 0;
    if (deleteC == true) deleteCnt_nm_1 = 0;
}

void localWave_nm_1(int p, int hop, int n) {
    lc_nm_1 = 0;

    tmp_nm_1.resize(L_nm_1.size());
    copy(L_nm_1.begin(), L_nm_1.end(), tmp_nm_1.begin());

    exNodes_nm_1 = 0;

    while (lc_nm_1 <= hop and exNodes_nm_1 <= hop) {
        int moveLoad = 0;
        int pi = md(p+lc_nm_1, n);

        if (tmp_nm_1[pi] > overthres_nm_1) {
            moveLoad = a_const_nm_1*(tmp_nm_1[pi]-thres_nm_1);
            moveLoad = (int)moveLoad;
        } else if (tmp_nm_1[pi] > 0 and tmp_nm_1[pi] > thres_nm_1) {
            int plc_nm_1 = md(p + lc_nm_1, n); // index in p+lc_nm_1
            moveLoad = tmp_nm_1[plc_nm_1] - thres_nm_1;
        }
        tmp_nm_1[pi] = L_nm_1[pi] - moveLoad;

        int pii = md(pi+1, n);
        tmp_nm_1[pii] = L_nm_1[pii] + moveLoad;

        exNodes_nm_1 = (int)(tmp_nm_1[pii]/thres_nm_1 - 1);

        exNodes_nm_1 = exNodes_nm_1 < 0 ? 0 : exNodes_nm_1;

        lc_nm_1 = lc_nm_1 + 1;
    }
    return ;
}

void remoteWave_nm_1(int p, int exNode, int n) {
    int p_lc_nm_1_1 = md(p, n); // index in p + lc_nm_1 + 1
    m_nm_1 = getIdxMin_nm_1(n);

    int tmp_nm_1Lm = L_nm_1[m_nm_1];
    int j = 0;

    while (tmp_nm_1Lm <= tmp_nm_1[p_lc_nm_1_1] and j <= exNode) {
        int m_j_1 = md(m_nm_1+j+1, n);

        tmp_nm_1Lm += L_nm_1[m_j_1];

        j = j + 1;
    }
    rNodes_nm_1 = j;
    return;
}

void nix_nm_1(int p, int load, int n) {
    int from = md(p, n);      //index of p
    int to   = md(from+1, n); //index of p+1

    if ( L_nm_1[from] >= load ) {
        L_nm_1[from] -= load;
        L_nm_1[to]   += load;

        moveCnt_nm_1 += load * moveCst_nm_1;                    // count opeation move
    }
    return;
}

void mig_nm_1(int far_node, int p, int load, int n) {
    // ------- nix (m->m-1, L[m]) ------------
    int neighbor_idx = md(far_node-1, n);
    int tmp_nm_1_fn = L_nm_1[far_node];

    L_nm_1.erase(L_nm_1.begin()+far_node);                //   <----------------------------------o               
                                                    //                                      |
    if (far_node < neighbor_idx) {                  //                                      |
        neighbor_idx -= 1;                          //                                      |
        neighbor_idx = md(neighbor_idx, n);         //                                      +--o
    }                                               //                                      |  |
                                                    //                                      |  V
    L_nm_1.insert(L_nm_1.begin()+neighbor_idx, tmp_nm_1_fn);   // <-----------------------------------<O
    moveCnt_nm_1 += tmp_nm_1_fn *(moveCst_nm_1);                   //    count operation (delete + insert) = Move

    far_node = neighbor_idx;
    nix_nm_1(far_node, L_nm_1[far_node], n);                    // nix(m -> m-1, L[m])

    L_nm_1.erase(L_nm_1.begin()+far_node);
    if (far_node < p) {
        p -= 1;
        p = md(p, n);
    }

    // -------  nix(p -> p+1 # p+1 is m that empty, load) -------------
    int next_p = md(p+1, n);
    L_nm_1.insert(L_nm_1.begin()+next_p, 0);  //After process m, value in m is empty
    nix_nm_1(p, load, n);              // nix(p -> p+1 # p+1 is m that empty, load)
    return;
}

void nixmig_nm_1(int p,int hop, int n) {
    localWave_nm_1(p, ttl_nm_1, n);
    if (exNodes_nm_1 > 0) {
        remoteWave_nm_1(p+lc_nm_1+1, exNodes_nm_1, n);
    }

    for (int i = 0; i <= lc_nm_1; ++i) {
        int pi = md(p+i, n); // index at p + i
        int load = 0;

        if (L_nm_1[pi] > overthres_nm_1) {
            load = a_const_nm_1 * (L_nm_1[pi]-thres_nm_1);
            load = (int)load;
        } else if (L_nm_1[pi] > 0 and L_nm_1[pi] > thres_nm_1) {
            load = L_nm_1[pi] - thres_nm_1;
        }

        nix_nm_1(pi, load, n);
    }

    if (rNodes_nm_1 > 0) {
        for (int i = 0; i <= rNodes_nm_1; ++i) {
            int plc_nm_1 = md(p+lc_nm_1, n);
            int value = (int)(tmp_nm_1[plc_nm_1]/rNodes_nm_1);

            int m_i_1 = md(m_nm_1+i+1, n);
            mig_nm_1(m_i_1, plc_nm_1, value, n); // m is minimum index

        }
    }
    return;
}

double getRatio_nm_1(int n) {
    double mn = (double)getMin_nm_1(n);
    double mx = (double)getMax_nm_1(n);
    if (abs((int)mn) <= 0.0000001) return mx;
    return mx/mn;
}

void print_nm_1(int no, int n) {
    char str[100];
    sprintf(str, "nm_%d.txt", no);
    ofstream out(str);
    for (int i  = 0; i < n; ++i) {
        out << i << ' ' << L_nm_1[i] << '\n';
    }
    out.close();
    return;
}

int query_nm_1(int s, int e, int n) {
    int tot = 0;        // tot = total ; it same last pos of interval
    int rlt = 0;        // num of node about [s, e]
    for (int i = 0; i < n; ++i) {
        tot += L_nm_1[i];
        if (tot > e) {  // last interval
            if (tot-L_nm_1[i]+1 <= e) ++rlt;
            break;
        }

        // live in Interval [s, e]
        if (s <= tot and tot <= e) {
            ++rlt;
        }
    }
    return rlt;
}

void insertLoad_nm_1(int idx, int& n) {
    L_nm_1[idx] += 1;
    addCnt_nm_1 += addCst_nm_1;                                // count opertaion add
    nixmig_nm_1(idx, ttl_nm_1, n);
}
// =====================================================
//
//
/* ########### Fucntion of Algor'Nixmig VERSION 2 ############# */

void setCnt_nm_2(bool all = 1, bool move = 1, bool add = 1, bool deleteC = 1) {
    if (all == true) cost_nm_2 = 0;
    if (move == true) moveCnt_nm_2 = 0;
    if (add == true) addCnt_nm_2 = 0;
    if (deleteC == true) deleteCnt_nm_2 = 0;
}

void localWave_nm_2(int p, int hop, int n) {
    lc_nm_2 = 0;

    tmp_nm_2.resize(L_nm_2.size());
    copy(L_nm_2.begin(), L_nm_2.end(), tmp_nm_2.begin());

    exNodes_nm_2 = 0;

    while (lc_nm_2 <= hop and exNodes_nm_2 <= hop) {
        int moveLoad = 0;
        int pi = md(p+lc_nm_2, n);

        if (tmp_nm_2[pi] > overthres_nm_2) {
            moveLoad = a_const_nm_2*(tmp_nm_2[pi]-thres_nm_2);
            moveLoad = (int)moveLoad;
        } else if (tmp_nm_2[pi] > 0 and tmp_nm_2[pi] > thres_nm_2) {
            int plc_nm_2 = md(p + lc_nm_2, n); // index in p+lc_nm_2
            moveLoad = tmp_nm_2[plc_nm_2] - thres_nm_2;
        }
        tmp_nm_2[pi] = L_nm_2[pi] - moveLoad;

        int pii = md(pi+1, n);
        tmp_nm_2[pii] = L_nm_2[pii] + moveLoad;

        exNodes_nm_2 = (int)(tmp_nm_2[pii]/thres_nm_2 - 1);

        exNodes_nm_2 = exNodes_nm_2 < 0 ? 0 : exNodes_nm_2;

        lc_nm_2 = lc_nm_2 + 1;
    }
    return ;
}

void remoteWave_nm_2(int p, int exNode, int n) {
    int p_lc_nm_2_2 = md(p, n); // index in p + lc_nm_2 + 1
    m_nm_2 = getIdxMin_nm_2(n);

    int tmp_nm_2Lm = L_nm_2[m_nm_2];
    int j = 0;

    while (tmp_nm_2Lm <= tmp_nm_2[p_lc_nm_2_2] and j <= exNode) {
        int m_j_2 = md(m_nm_2+j+1, n);

        tmp_nm_2Lm += L_nm_2[m_j_2];

        j = j + 1;
    }
    rNodes_nm_2 = j;
    return;
}

void nix_nm_2(int p, int load, int n) {
    int from = md(p, n);      //index of p
    int to   = md(from+1, n); //index of p+1

    if ( L_nm_2[from] >= load ) {
        L_nm_2[from] -= load;
        L_nm_2[to]   += load;

        moveCnt_nm_2 += load * moveCst_nm_2;                    // count opeation move
    }
    return;
}

void mig_nm_2(int far_node, int p, int load, int n) {
    // ------- nix (m->m-1, L[m]) ------------
    int neighbor_idx = md(far_node-1, n);
    int tmp_nm_2_fn = L_nm_2[far_node];

    L_nm_2.erase(L_nm_2.begin()+far_node);                //   <----------------------------------o               
                                                    //                                      |
    if (far_node < neighbor_idx) {                  //                                      |
        neighbor_idx -= 1;                          //                                      |
        neighbor_idx = md(neighbor_idx, n);         //                                      +--o
    }                                               //                                      |  |
                                                    //                                      |  V
    L_nm_2.insert(L_nm_2.begin()+neighbor_idx, tmp_nm_2_fn);   // <-----------------------------------<O
    moveCnt_nm_2 += tmp_nm_2_fn *(moveCst_nm_2);                   //    count operation (delete + insert) = Move

    far_node = neighbor_idx;
    nix_nm_2(far_node, L_nm_2[far_node], n);                    // nix(m -> m-1, L[m])

    L_nm_2.erase(L_nm_2.begin()+far_node);
    if (far_node < p) {
        p -= 1;
        p = md(p, n);
    }

    // -------  nix(p -> p+1 # p+1 is m that empty, load) -------------
    int next_p = md(p+1, n);
    L_nm_2.insert(L_nm_2.begin()+next_p, 0);  //After process m, value in m is empty
    nix_nm_2(p, load, n);              // nix(p -> p+1 # p+1 is m that empty, load)
    return;
}

void nixmig_nm_2(int p,int hop, int n) {
    localWave_nm_2(p, ttl_nm_2, n);
    if (exNodes_nm_2 > 0) {
        remoteWave_nm_2(p+lc_nm_2+1, exNodes_nm_2, n);
    }

    for (int i = 0; i <= lc_nm_2; ++i) {
        int pi = md(p+i, n); // index at p + i
        int load = 0;

        if (L_nm_2[pi] > overthres_nm_2) {
            load = a_const_nm_2 * (L_nm_2[pi]-thres_nm_2);
            load = (int)load;
        } else if (L_nm_2[pi] > 0 and L_nm_2[pi] > thres_nm_2) {
            load = L_nm_2[pi] - thres_nm_2;
        }

        nix_nm_2(pi, load, n);
    }

    if (rNodes_nm_2 > 0) {
        for (int i = 0; i <= rNodes_nm_2; ++i) {
            int plc_nm_2 = md(p+lc_nm_2, n);
            int value = (int)(tmp_nm_2[plc_nm_2]/rNodes_nm_2);

            int m_i_2 = md(m_nm_2+i+1, n);
            mig_nm_2(m_i_2, plc_nm_2, value, n); // m is minimum index

        }
    }
    return;
}

double getRatio_nm_2(int n) {
    double mn = (double)getMin_nm_2(n);
    double mx = (double)getMax_nm_2(n);
    if (abs((int)mn) <= 0.0000001) return mx;
    return mx/mn;
}

void print_nm_2(int no, int n) {
    char str[100];
    sprintf(str, "nm_%d.txt", no);
    ofstream out(str);
    for (int i  = 0; i < n; ++i) {
        out << i << ' ' << L_nm_2[i] << '\n';
    }
    out.close();
    return;
}

int query_nm_2(int s, int e, int n) {
    int tot = 0;        // tot = total ; it same last pos of interval
    int rlt = 0;        // num of node about [s, e]
    for (int i = 0; i < n; ++i) {
        tot += L_nm_2[i];
        if (tot > e) {  // last interval
            if (tot-L_nm_2[i]+1 <= e) ++rlt;
            break;
        }

        // live in Interval [s, e]
        if (s <= tot and tot <= e) {
            ++rlt;
        }
    }
    return rlt;
}

void insertLoad_nm_2(int idx, int& n) {
    L_nm_2[idx] += 1;
    addCnt_nm_2 += addCst_nm_2;                                // count opertaion add
    nixmig_nm_2(idx, ttl_nm_2, n);
}
// ==========================================================

void sleep(double sec = 1021) { clock_t s = clock(); while(clock() - s < CLOCKS_PER_SEC * sec); }

/* ########### Fucntion MAIN of Process ############# */
int main(){
    srand(time(0));

    int T = 1;
    while (T--) {

        int n        = 1<<10;        // amount of Tank or L
        int N_thing  = 100000;      // amount of insertion
        // int numQuery = 100;         // amount of Query Interval in Each algorithm

        L_jc.clear();
        L_kg.clear();
        L_nm_1.clear();
        L_nm_2.clear();

        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            L_jc.push_back(0);
            L_kg.push_back(0);
            L_nm_1.push_back(0);
            L_nm_2.push_back(0);
        }

        // ====================== Create RATIO and Ratio File =========================
        char nameFileRatio[100];
        char nameFileCount[100];
        sprintf(nameFileRatio, "ratio(%d).txt", n);
        sprintf(nameFileCount, "count(%d).txt", n);
        ofstream outRatio(nameFileRatio);
        ofstream outCount(nameFileCount);

        ofstream descriptRatio("Description Ratio File.txt");
        descriptRatio << "In File ratio(i).txt has 5 Colume.\n";
        descriptRatio << "Colume 1 : Number of Line.\n";
        descriptRatio << "Colume 2 : Ratio of JC's algorithm.\n";
        descriptRatio << "Colume 3 : Ratio of Karger's algorithm.\n";
        descriptRatio << "Colume 4 : Ratio of NixMig's algorithm [version 1].\n";
        descriptRatio << "Colume 5 : Ratio of NixMig's algorithm [version 2].\n";
        descriptRatio.close();

        ofstream descriptCount("Description Count File.txt");
        descriptCount << "In File count(i).txt has 5 Colume.\n";
        descriptCount << "Colume 1 : Number of Line.\n";
        descriptCount << "Colume 2 : Count of JC's algorithm.\n";
        descriptCount << "Colume 3 : Count of Karger's algorithm.\n";
        descriptCount << "Colume 4 : Count of NixMig's algorithm [version 1].\n";
        descriptCount << "Colume 5 : Count of NixMig's algorithm [version 2].\n";
        descriptCount << "And the last line of this file, it has List's total of Count Operation of all algorithm.\n";
        descriptCount.close();

        double ratio_jc = getRatio_jc(n);
        double ratio_kg = getRatio_kg(n);
        double ratio_nm_1 = getRatio_nm_1(n);
        double ratio_nm_2 = getRatio_nm_2(n);

        outRatio << "0 " << ratio_jc << ' ' << ratio_kg << ' ' << ratio_nm_1 << ' ' << ratio_nm_2 << '\n';
        
        setCnt_jc();
        setCnt_kg();
        setCnt_nm_1();
        setCnt_nm_2();

        ll total_jc = 0LL;
        ll total_kg = 0LL;
        ll totaL_nm_1 = 0LL;
        ll total_nm_2 = 0LL;

        int part_continue = (N_thing+99)/100;

        printf(" 1: \t+ ------------------------------------------------------ +\n");
        printf(" 2: \t|      06 - max node all algor ratio count operation     |\n");
        printf(" 3: \t+ ------------------------------------------------------ +\n 4: \n");

        printf(" 5: n = %d\n", n);
        printf(" 6: number of operation = %d\n 7: \n", N_thing);

        // input data
        for (int i = 0; i < N_thing; ++i) {
            setCnt_jc(0);
            setCnt_kg(0);
            setCnt_nm_1(0);
            setCnt_nm_2(0);

            insertLoad_jc(getIdxMax_jc(n), n);
            insertLoad_kg(getIdxMax_kg(n), n);
            insertLoad_nm_1(getIdxMax_nm_1(n), n);
            insertLoad_nm_2(getIdxMax_nm_2(n), n);

            // ====================== Create RATIO File =========================
            double ratio_jc = getRatio_jc(n);
            double ratio_kg = getRatio_kg(n);
            double ratio_nm_1 = getRatio_nm_1(n);
            double ratio_nm_2 = getRatio_nm_2(n);

            cost_jc = moveCnt_jc + addCnt_jc + deleteCnt_jc;
            cost_kg = moveCnt_kg + addCnt_kg + deleteCnt_kg;
            cost_nm_1 = moveCnt_nm_1 + addCnt_nm_1 + deleteCnt_nm_1;
            cost_nm_2 = moveCnt_nm_2 + addCnt_nm_2 + deleteCnt_nm_2;

            outRatio << i+1 << ' ' << ratio_jc << ' ' << ratio_kg << ' ' << ratio_nm_1 << ' ' << ratio_nm_2 << '\n';
            outCount << i+1 << ' ' <<  cost_jc << ' ' <<  cost_kg << ' ' <<  cost_nm_1 << ' ' <<  cost_nm_2 << '\n';

            total_jc += cost_jc;
            total_kg += cost_kg;
            totaL_nm_1 += cost_nm_1;
            total_nm_2 += cost_nm_2;

            // print_jc(i, n);
            // print_kg(i, n);
            // print_nm_1(i, n);
            // print_nm_2(i, n);

            // see process of program when it runs.
            int now = i+1;
            string tab_status;
            if (now % part_continue == 0) {
                int ten_part = (int)(now/(double)part_continue * 0.25);
                tab_status = "";
                while (ten_part--) tab_status += "P";
                
                printf("\r 8: \tSuccess |%-25s| : %3d%% ...", tab_status.c_str(),  (int)min(100.0/N_thing * now, 100.0));
                if(n==2)sleep(0.0005);
                fflush(stdout);
            }
        }
        printf("\n 9: \n");
        
        // =======check amount of thing in node  =========
        int sum_chk[4] = {};
        for (int i = 0; i < n; ++i) {
            sum_chk[0] += L_jc[i];
            sum_chk[1] += L_kg[i];
            sum_chk[2] += L_nm_1[i];
            sum_chk[3] += L_nm_2[i];
        }

        for (int i = 0; i < 4; ++i) {
            if (sum_chk[i] != N_thing) 
                printf("Program is Bug!\n");
            assert(sum_chk[i] == N_thing);
        }

        printf("10: End Process insertion : closed ratio file...\n");
        
        outCount << "total_jc = " << total_jc << '\n';
        outCount << "total_kg = " << total_kg << '\n';
        outCount << "totaL_nm_1 = " << totaL_nm_1 << '\n';
        outCount << "total_nm_2 = " << total_nm_2 << '\n';

        printf("11: End Process insertion : closed count file...\n");

        // ==== Create RATIO and Count File ====
        outRatio.close();
        outCount.close();
        // =====================================
    }

    printf("12: End all process.\n13: \n");
    printf("14: Pass any key to close program ...");
    getchar();
    printf("15:");
    return 0;
}