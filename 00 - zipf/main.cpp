#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <set>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include "MersenneTwister.h"

using namespace std;

const int N = 1000;

multiset<int> load;
multiset<int>::iterator l_it; // load iterator

vector<int> tnk; // tank

const double alpha = 5.464;
const double beta = 4.1;

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

void clearVal(int val) {
    l_it = load.lower_bound(val);
    if (l_it == load.end()) return;
    load.erase(l_it);
    return;
}

void minbalance(int u, int& n) {
    int v = getIdxMin(n);
    if (tnk[u] > alpha*tnk[v]) {
        int l = md(v-1, n);
        int r = md(v+1, n);
        int z = tnk[l] < tnk[r] ? l : r;

        clearVal(tnk[u]);
        clearVal(tnk[v]);
        clearVal(tnk[z]);

        tnk[z] += tnk[v];

#define start tnk.begin()

        tnk.erase(start+v);

        if (v < u)
            --u;

        int remain_up = (tnk[u]+1)/2;
        int remain_down =   tnk[u]/2;
        tnk[u] = remain_up;

        tnk.insert(start+u, remain_down);

        load.insert(remain_down);
        load.insert(remain_up);
        load.insert(tnk[z]);
    }
    return;
}

void insertLoad(int idx, int& n) {
    clearVal(tnk[idx]);
    tnk[idx] += 1;
    load.insert(tnk[idx]);
    minbalance(idx, n);
}

void split(int u, int n) {
    int w = getIdxMax(n);
    if (tnk[u] < tnk[w]/beta) {
        int l = md(u-1, n);
        int r = md(u+1, n);
        int z = tnk[l] < tnk[r] ? l : r;

        if (tnk[z] <= 2*tnk[w]/beta) {
            clearVal(tnk[u]);
            clearVal(tnk[z]);
            clearVal(tnk[w]);

            tnk[z] += tnk[u];
            tnk.erase(tnk.begin() + u);

            if (u < w)
                --w;

            int remain_up = (tnk[w]+1)/2;
            int remain_down = tnk[w] / 2;

            tnk[w] = remain_up;
            tnk.insert(tnk.begin() + w, remain_down);

            load.insert(tnk[z]);
            load.insert(tnk[w]);
            load.insert(remain_down);
        } else {
            int sum = tnk[u] + tnk[z];

            clearVal(tnk[u]);
            clearVal(tnk[z]);

            tnk[u] = (sum+1)/2;
            tnk[z] = sum /2;

            load.insert(tnk[u]);
            load.insert(tnk[z]);
        }
    }
}

double getRatio(int n) {
    double mn = (double)getMin(n);
    double mx = (double)getMax(n);
    if (mn <= 0.0000001) return mx;
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

int main() {
    srand(time(0));
    int T = 1;

    double alpha = 1.00;

    while (T--) {
        int n;
        n = 1024;
        //scanf("%d", &n);

        tnk.clear();
        load.clear();

        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            tnk.push_back(0);
            load.insert(0);
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

        sprintf(nameFileRatio, "ratio(%d).txt", n);
        ofstream outRatio(nameFileRatio);

        outRatio << "0 " << getRatio(n) << '\n';

        // input data in max load in tank
        for (int i = 0; i < 100000; ++i) {
            int idx = zipf(alpha, n)-1;
            assert(0 <= idx && idx <= n-1);
            insertLoad(idx, n);
            outRatio << i+1 << ' ' << getRatio(n) << '\n';

            print(i, n);
        }
        outRatio.close();
    }
    return 0;
}
