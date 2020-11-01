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

int cost_jc;        // use count cost of all operation in algorithm
int moveCnt_jc;     // use count only move operation
int addCnt_jc;      // use count only add operation
int deleteCnt_jc;   // use count only delete operation

const int moveCst_jc   = 1;      // move cost
const int addCst_jc    = 1;      // add cost
const int deleteCst_jc = 1;      // delete cost

int cost_kg;        // use count cost of all operation in algorithm
int moveCnt_kg;     // use count only move operation
int addCnt_kg;      // use count only add operation
int deleteCnt_kg;   // use count only delete operation

const int moveCst_kg   = 1;      // move cost
const int addCst_kg    = 1;      // add cost
const int deleteCst_kg = 1;      // delete cost

int cost_nm;        // use count cost of all operation in algorithm
int moveCnt_nm;     // use count only move operation
int addCnt_nm;      // use count only add operation
int deleteCnt_nm;   // use count only delete operation

const int moveCst_nm   = 1;      // move cost
const int addCst_nm    = 1;      // add cost
const int deleteCst_nm = 1;      // delete cost

double alphaZipf = 1.00;

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
int getMax_jc(int n) {
    return L_jc[getIdxMax_jc(n)];
}

int getMax_kg(int n) {
    return L_kg[getIdxMax_kg(n)];
}

int getMax_nm(int n) {
    return L_nm[getIdxMax_nm(n)];
}
// =====================================================

/* = = ========= Fucntion get Value Min ========= = = */
int getMin_jc(int n) {
    return L_jc[getIdxMin_jc(n)];
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
MTRand mtrand1;

int zipf(double alpha, int n) {
    static int first = 1;         // Static first time flag
    static double c = 0;          // Normalization constant
    double z;                     // Uniform random number (0 < z < 1)
    double sum_prob;              // Sum of probabilities
    double zipf_value;            // Computed exponential value to be returned
    int    i;                     // Loop counter

    // Compute normalization constant on first call only
    if (first == 1) {
        for (i=1; i<=n; i++)
            c = c + (1.0 / pow((double) i, alpha));
        c = 1.0 / c;
        first = 0;
    }

    // Pull a uniform random number (0 < z < 1)
    do {
        z = mtrand1.rand();
    } while ((z == 0) || (z == 1));

    // Map z to the value
    sum_prob = 0;
    for (i=1; i<=n; i++) {
        sum_prob = sum_prob + c / pow((double) i, alpha);
        if (sum_prob >= z) {
            zipf_value = i;
            break;
        }
    }

    // Assert that zipf_value is between 1 and N
    assert((zipf_value >=1) && (zipf_value <= n));

    return(zipf_value);
}

int rt(int n) { // random [0, n-1]
    int idx = zipf(alphaZipf, n) - 1;
    return idx;
}
// =====================================================

/* ############## Fucntion of Algor'jc ############## */

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

/* ########### Fucntion of Algor'Nixmig ############# */

void setCnt_nm(bool all = 1, bool move = 1, bool add = 1, bool deleteC = 1) {
    if (all == true) cost_nm = 0;
    if (move == true) moveCnt_nm = 0;
    if (add == true) addCnt_nm = 0;
    if (deleteC == true) deleteCnt_nm = 0;
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
    addCnt_nm += addCst_nm;                                // count opertaion add
    nixmig(idx, ttl, n);
}
// =====================================================
void sleep(double sec = 1021) { clock_t s = clock(); while(clock() - s < CLOCKS_PER_SEC * sec); }
/* ########### Fucntion MAIN of Process ############# */
int main(){
    srand(time(0));

    int T = 1;
    while (T--) {

        int n        = 1024;        // amount of Tank or L
        int N_thing  = 100000;      // amount of insertion
        // int numQuery = 100;         // amount of Query Interval in Each algorithm

        L_jc.clear();
        L_kg.clear();
        L_nm.clear();

        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            L_jc.push_back(0);
            L_kg.push_back(0);
            L_nm.push_back(0);
        }

        // ====================== Create RATIO and Ratio File =========================
        char nameFileRatio[100];
        char nameFileCount[100];
        sprintf(nameFileRatio, "ratio(%d).txt", n);
        sprintf(nameFileCount, "count(%d).txt", n);
        ofstream outRatio(nameFileRatio);
        ofstream outCount(nameFileCount);

        ofstream descriptRatio("Description Ratio File.txt");
        descriptRatio << "In File ratio(i).txt has 4 Colume.\n";
        descriptRatio << "Colume 1 : Number of Line.\n";
        descriptRatio << "Colume 2 : Ratio of JC's algorithm.\n";
        descriptRatio << "Colume 3 : Ratio of Karger's algorithm.\n";
        descriptRatio << "Colume 4 : Ratio of NixMig's algorithm.\n";
        descriptRatio.close();

        ofstream descriptCount("Description Count File.txt");
        descriptCount << "In File count(i).txt has 4 Colume.\n";
        descriptCount << "Colume 1 : Number of Line.\n";
        descriptCount << "Colume 2 : Count of JC's algorithm.\n";
        descriptCount << "Colume 3 : Count of Karger's algorithm.\n";
        descriptCount << "Colume 4 : Count of NixMig's algorithm.\n";
        descriptCount << "And the last line of this file, it has List's total of Count Operation of all algorithm.\n";
        descriptCount.close();

        double ratio_jc = getRatio_jc(n);
        double ratio_kg = getRatio_kg(n);
        double ratio_nm = getRatio_nm(n);

        outRatio << "0 " << ratio_jc << ' ' << ratio_kg << ' ' << ratio_nm << '\n';
        
        setCnt_jc();
        setCnt_kg();
        setCnt_nm();

        ll total_jc = 0LL;
        ll total_kg = 0LL;
        ll total_nm = 0LL;

        int part_continue = (N_thing+99)/100;

        printf(" 1: \t+ ------------------------------------------------------ +\n");
        printf(" 2: \t|        05 - zipf all algor ratio count operation       |\n");
        printf(" 3: \t+ ------------------------------------------------------ +\n 4: \n");

        printf(" 5: n = %d\n", n);
        printf(" 6: number of operation = %d\n 7: \n", N_thing);
        // input data

        
        for (int i = 0; i < N_thing; ++i) {
            int idx = rt(n); // declare by zipf function
            if (idx < 0 or idx > n-1) {
                printf("\t\t\t Random index const Point is wrong.\n");
                assert(0);
                return 0;
            }

            setCnt_jc(0);
            setCnt_kg(0);
            setCnt_nm(0);

            insertLoad_jc(idx, n);
            insertLoad_kg(idx, n);
            insertLoad_nm(idx, n);

            // ====================== Create RATIO File =========================
            double ratio_jc = getRatio_jc(n);
            double ratio_kg = getRatio_kg(n);
            double ratio_nm = getRatio_nm(n);

            cost_jc = moveCnt_jc + addCnt_jc + deleteCnt_jc;
            cost_kg = moveCnt_kg + addCnt_kg + deleteCnt_kg;
            cost_nm = moveCnt_nm + addCnt_nm + deleteCnt_nm;

            outRatio << i+1 << ' ' << ratio_jc << ' ' << ratio_kg << ' ' << ratio_nm << '\n';
            outCount << i+1 << ' ' <<  cost_jc << ' ' <<  cost_kg << ' ' <<  cost_nm << '\n';

            total_jc += cost_jc;
            total_kg += cost_kg;
            total_nm += cost_nm;

            // print_jc(i, n);
            // print_kg(i, n);
            // print_nm(i, n);

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
        int sum_chk[3] = {};
        for (int i = 0; i < n; ++i) {
            sum_chk[0] += L_jc[i];
            sum_chk[1] += L_kg[i];
            sum_chk[2] += L_nm[i];
        }
        for (int i = 0; i < 3; ++i) {
            if (sum_chk[i] != N_thing) {
                printf("Program is Bug!\n");
                printf("i = %d\n", i);
                printf("sum_chk[i] = %d\n", sum_chk[i]);
                printf("N_thing = %d\n", N_thing);
            }
            assert(sum_chk[i] == N_thing);
        }

        printf("10: End Process insertion : closed ratio file...\n");
        
        outCount << "total_jc = " << total_jc << '\n';
        outCount << "total_kg = " << total_kg << '\n';
        outCount << "total_nm = " << total_nm << '\n';

        printf("11: End Process insertion : closed count file...\n");

        // ==== Create RATIO and Count File ====
        outRatio.close();
        outCount.close();
        // ===========================

        // ===================== Create Query File ===========================
        // itv = [i]n[t]er[v]al
        /*int itv[] = {10, 20, 50, 100, 200, 500, 1000, 2000}; 
        
        int size_itv = sizeof (itv) / sizeof(int);

        for (int i_range = 0; i_range < size_itv; ++i_range) {
            int range = itv[i_range]; // get range in round i
            
            char nameFileQuery[100];  // make name file query
            sprintf(nameFileQuery, "query(%d).txt", range);
            ofstream outQuery(nameFileQuery);

            for (int t = 1; t <= numQuery; ++t) {
                // random in range [1, N_thing-range+1]; N_think >= 2000
                int pos_start = 1+rt(N_thing-range+1);

                int n_jc = query_jc(pos_start, pos_start+range-1, n); // x
                int n_kg = query_kg(pos_start, pos_start+range-1, n); // y
                int n_nm = query_nm(pos_start, pos_start+range-1, n); // z
                outQuery << t << ' ' << n_jc << ' '  << n_kg << ' ' << n_nm << '\n';
            }

            outQuery.close();
            printf("End Process Query(%d)\n", range);
        }
        printf("End Process Query\n");*/
    }

    printf("12: End all process.\n13: \n");
    printf("14: Pass any key to close program ...");
    getchar();
    printf("15:");
    return 0;
}