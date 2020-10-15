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

vector<int> L_js; // tank of JS => L_js

vector<int> L_kg; // tank of Karger => L_kg

vector<int> L_nm; // tank of nixmig => L_mn
vector<int> tmp;  // tmpL in algorithm nixmig

const double alpha = 5.464;
const double beta = 4.1;
const double eps = 0.023;

const int ttl   = 5;
const int thres = 60;
const int overThres = 400;
const double a_const = 0.5;

int exNodes = 0;
int rNodes = 0;
int lc = 0;
int m; // min index

/* = = ========== Fucntion get Idx Max ========== = = */
int getIdxMax_js(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_js[i] > L_js[idx])
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

int getIdxMax_nm(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_nm[i] > L_nm[idx])
            idx = i;
    }
    return idx;
}
// =====================================================

/* = = ========== Fucntion get Idx Min ========== = = */
int getIdxMin_js(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_js[i] < L_js[idx])
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

int getIdxMin_nm(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_nm[i] < L_nm[idx])
            idx = i;
    }
    return idx;
}
// =====================================================

/* = = ========= Fucntion get Value Max ========= = = */
int getMax_js(int n) {
    return L_js[getIdxMax_js(n)];
}

int getMax_kg(int n) {
    return L_kg[getIdxMax_kg(n)];
}

int getMax_nm(int n) {
    return L_nm[getIdxMax_nm(n)];
}
// =====================================================

/* = = ========= Fucntion get Value Min ========= = = */
int getMin_js(int n) {
    return L_js[getIdxMin_js(n)];
}

int getMin_kg(int n) {
    return L_kg[getIdxMin_kg(n)];
}

int getMin_nm(int n) {
    return L_nm[getIdxMin_nm(n)];
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
void minbalance(int u, int& n) {
    int v = getIdxMin_js(n);
    if (L_js[u] > alpha*L_js[v]) {
        int l = md(v-1, n);
        int r = md(v+1, n);
        int z = L_js[l] < L_js[r] ? l : r;

        L_js[z] += L_js[v];

        #define start L_js.begin()

        L_js.erase(start+v);

        if (v < u)
            --u;

        int remain_up = (L_js[u]+1)/2;
        int remain_down =   L_js[u]/2;
        L_js[u] = remain_up;

        L_js.insert(start+u, remain_down);
    }
    return;
}

void split(int u, int n) {
    int w = getIdxMax_js(n);
    if (L_js[u] < L_js[w]/beta) {
        int l = md(u-1, n);
        int r = md(u+1, n);
        int z = L_js[l] < L_js[r] ? l : r;

        if (L_js[z] <= 2*L_js[w]/beta) {

            L_js[z] += L_js[u];
            L_js.erase(L_js.begin() + u);

            if (u < w)
                --w;

            int remain_up = (L_js[w]+1)/2;
            int remain_down = L_js[w] / 2;

            L_js[w] = remain_up;
            L_js.insert(L_js.begin() + w, remain_down);
        }else{
            int sum = L_js[u] + L_js[z];

            L_js[u] = (sum+1)/2;
            L_js[z] = sum /2;
        }
    }
}

double getRatio_js(int n) {
    double mn = (double)getMin_js(n);
    double mx = (double)getMax_js(n);
    if (mn <= 0.0000001) return mx;
    return mx/mn;
}

void print_js(int no, int n) {
    char str[100];
    sprintf(str, "js_%d.txt", no);
    ofstream out(str);
    for (int i  = 0; i < n; ++i) {
        out << i << ' ' << L_js[i] << '\n';
    }
    out.close();
    return;
}

int query_js(int s, int e, int n) {
    int tot = 0;        // tot = total ; it same last pos of interval
    int rlt = 0;        // num of node about [s, e]
    for (int i = 0; i < n; ++i) {
        tot += L_js[i];
        if (tot > e) {  // last interval
            if (tot-L_js[i]+1 <= e) ++rlt;
            break;
        }

        // live in Interval [s, e]
        if (s <= tot and tot <= e) {
            ++rlt;
        }
    }
    return rlt;
}

void insertLoad_js(int idx, int& n) {
    L_js[idx] += 1;
    minbalance(idx, n);
}
// =====================================================

/* ########### Fucntion of Algor'karger ############# */
void divideEq(int a, int b) {
    int n = (L_kg[a]+L_kg[b]+1)/2;
    int m = (L_kg[a]+L_kg[b])/2;
    L_kg.erase(L_kg.begin()+a);
    if (a < b)
        --b;
    L_kg.erase(L_kg.begin()+b);
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

                int n = (L_kg[i]+1)/2;
                int m = (L_kg[i])/2;

                L_kg[j_right] += L_kg[j];
                L_kg.erase(L_kg.begin()+j);

                if (j < i)
                    --i;

                L_kg.erase(L_kg.begin()+i);
                L_kg.insert(L_kg.begin()+i, m);
                L_kg.insert(L_kg.begin()+i, n);
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
    itemBalance(idx, n);
}
// =====================================================

/* ########### Fucntion of Algor'Nixmig ############# */
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
    m = getIdxMin_nm(n);

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
    }
    return;
}

void mig(int far_node, int p, int load, int n) {
    // ------- nix (m->m-1, L[m]) ------------
    int neighbor_idx = md(far_node-1, n);
    int tmp_fn = L_nm[far_node];

    L_nm.erase(L_nm.begin()+far_node);
    if (far_node < neighbor_idx) {
        neighbor_idx -= 1;
        neighbor_idx = md(neighbor_idx, n);
    }

    L_nm.insert(L_nm.begin()+neighbor_idx, tmp_fn);
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
        for (int i = 0; i <= rNodes; ++i) {
            int plc = md(p+lc, n);
            int value = (int)(tmp[plc]/rNodes);

            int m_i_1 = md(m+i+1, n);
            mig(m_i_1, plc, value, n); // m is minimum index

        }
    }
    return;
}

double getRatio_nm(int n) {
    double mn = (double)getMin_nm(n);
    double mx = (double)getMax_nm(n);
    if (abs((int)mn) <= 0.0000001) return mx;
    return mx/mn;
}

void print_nm(int no, int n) {
    char str[100];
    sprintf(str, "nm_%d.txt", no);
    ofstream out(str);
    for (int i  = 0; i < n; ++i) {
        out << i << ' ' << L_nm[i] << '\n';
    }
    out.close();
    return;
}

int query_nm(int s, int e, int n) {
    int tot = 0;        // tot = total ; it same last pos of interval
    int rlt = 0;        // num of node about [s, e]
    for (int i = 0; i < n; ++i) {
        tot += L_nm[i];
        if (tot > e) {  // last interval
            if (tot-L_nm[i]+1 <= e) ++rlt;
            break;
        }

        // live in Interval [s, e]
        if (s <= tot and tot <= e) {
            ++rlt;
        }
    }
    return rlt;
}

void insertLoad_nm(int idx, int& n) {
    L_nm[idx] += 1;
    nixmig(idx, ttl, n);
}
// =====================================================

/* ########### Fucntion MAIN of Process ############# */
int main(){
    srand(time(0));

    int T = 1;
    while (T--) {

        int n        = 1024;        // amount of Tank or L
        int N_thing  = 100000;      // amount of insertion
        int numQuery = 100;         // amount of Query Interval in Each algorithm

        L_js.clear();
        L_kg.clear();
        L_nm.clear();

        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            L_js.push_back(0);
            L_kg.push_back(0);
            L_nm.push_back(0);
        }
        
        /*
        // ====================== Create RATIO File =========================
        char nameFileRatio[100];
        sprintf(nameFileRatio, "ratio(%d).txt", n);
        ofstream outRatio(nameFileRatio);

        
        double ratio_js = getRatio_js(n)
        double ratio_kg = getRatio_kg(n)
        double ratio_nm = getRatio_nm(n)

        outRatio << "0 " << ratio_js << ' ' << ratio_kg << ' ' << ratio_nm << '\n';
        */

        // input data
        for (int i = 0; i < N_thing; ++i) {
            int idx = rt(n); // random [0, n-1]
            
            insertLoad_js(idx, n);
            insertLoad_kg(idx, n);
            insertLoad_nm(idx, n);

            /*
            // ====================== Create RATIO File =========================
            double ratio_js = getRatio_js(n);
            double ratio_kg = getRatio_kg(n);
            double ratio_nm = getRatio_nm(n);

            outRatio << i+1 << ' ' << ratio_js << ' ' << ratio_kg << ' ' << ratio_nm << '\n';
            
            // print_js(i, n);
            // print_kg(i, n);
            // print_nm(i, n);
            */
        }
        printf("End Process insertion\n");
        
        // ======= check amount of thing in node  =========
        int sum_chk[3] = {};
        for (int i = 0; i < n; ++i) {
            sum_chk[0] += L_js[i];
            sum_chk[1] += L_kg[i];
            sum_chk[2] += L_nm[i];
        }
        for (int i = 0; i < 3; ++i)
            assert(sum_chk[i] == N_thing);

        // ==== Create RATIO File ====
        //outRatio.close();
        // ===========================

        // ===================== Create Query File ===========================
        // itv = [i]n[t]er[v]al
        int itv[] = {10, 20, 50, 100, 200, 500, 1000, 2000}; 
        
        int size_itv = sizeof (itv) / sizeof(int);

        for (int i_range = 0; i_range < size_itv; ++i_range) {
            int range = itv[i_range]; // get range in round i
            
            char nameFileQuery[100];  // make name file query
            sprintf(nameFileQuery, "query(%d).txt", range);
            ofstream outQuery(nameFileQuery);

            for (int t = 1; t <= numQuery; ++t) {
                // random in range [1, N_thing-range+1]; N_think >= 2000
                int pos_start = 1+rt(N_thing-range+1);

                int n_js = query_js(pos_start, pos_start+range-1, n); // x
                int n_kg = query_kg(pos_start, pos_start+range-1, n); // y
                int n_nm = query_nm(pos_start, pos_start+range-1, n); // z
                outQuery << t << ' ' << n_js << ' '  << n_kg << ' ' << n_nm << '\n';
            }

            outQuery.close();
            printf("End Process Query(%d)\n", range);
        }
        printf("End Process Query\n");
    }

    printf("End Process");
    return 0;
}