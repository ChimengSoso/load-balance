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
#include <algorithm>
#include <random>       // default_random_engine
#include <chrono>       // chrono::system_clock

using namespace std;

vector<int> tnk; // tank
vector<int> shuf_tnk;

const double alpha = 5.464;
const double beta = 4.1;
const double eps = 0.023;

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

int rt(int n) { // random [0, n-1]
    MTRand m;
    return m.randInt(n-1);
}

void divideEq(int a, int b) {
    int n = (tnk[a]+tnk[b]+1)/2;
    int m = (tnk[a]+tnk[b])/2;
    tnk.erase(tnk.begin()+a);
    if (a < b)
        --b;
    tnk.erase(tnk.begin()+b);
    tnk.insert(tnk.begin()+a, m);
    tnk.insert(tnk.begin()+a, n);
    return;
}

void itemBalance(int i, int& n) {
    int j;

    do {
        j = rt(n);
    } while (j == i);

    if (tnk[i] <= eps*tnk[j] || tnk[j] <= eps*tnk[i]) {
        if (i == j+1) {
            divideEq(i, j);
        }else{ // case 2: i != j+1
            if (tnk[j+1] > tnk[i]) {
                divideEq(i, j);
            }else{
                int j_right =  md(j+1, n);

                int n = (tnk[i]+1)/2;
                int m = (tnk[i])/2;

                tnk[j_right] += tnk[j];
                tnk.erase(tnk.begin()+j);

                if (j < i)
                    --i;

                tnk.erase(tnk.begin()+i);
                tnk.insert(tnk.begin()+i, m);
                tnk.insert(tnk.begin()+i, n);
            }
        }
    }
    return;
}

void insertLoad(int idx, int& n) {
    tnk[idx] += 1;
    itemBalance(idx, n);
}

double getRatio(int n) {
    double mn = (double)getMin(n);
    double mx = (double)getMax(n);
    if (abs(mn) <= 0.0000001) return mx;
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

    double alphaZipf = 1.00;

    while (T--) {
        int n;
        n = 1024;
        //scanf("%d", &n);

        tnk.clear();
        shuf_tnk.clear();

        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            tnk.push_back(0);
            shuf_tnk.push_back(i);
        }

        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(shuf_tnk.begin(), shuf_tnk.end(), default_random_engine(seed));

        char nameFileRatio[100];

        sprintf(nameFileRatio, "ratio(%d).txt", n);
        ofstream outRatio(nameFileRatio);

        outRatio << "0 " << getRatio(n) << '\n';

        // input data in max load in tank
        for (int i = 0; i < 100000; ++i) {
            int idx = zipf(alphaZipf, n) - 1;
            assert(0 <= idx && idx <= n-1);
            insertLoad(shuf_tnk[idx], n);
            outRatio << i+1 << ' ' << getRatio(n) << '\n';
            print(i, n);
        }
        outRatio.close();
        printf("End Process\n");
    }
    return 0;
}
