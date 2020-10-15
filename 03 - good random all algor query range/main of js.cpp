#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <set>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include "MersenneTwister.h"
using namespace std;

const int N = 1000;

vector<int> L_js; // tank of JS => L_js

const double alpha = 5.464;
const double beta = 4.1;

int getIdxMax_js(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_js[i] > L_js[idx])
            idx = i;
    }
    return idx;
}

int getIdxMin_js(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_js[i] < L_js[idx])
            idx = i;
    }
    return idx;
}

int getMax_js(int n) {
    return L_js[getIdxMax_js(n)];
}

int getMin_js(int n) {
    return L_js[getIdxMin_js(n)];
}

int md(int n, int m) { // mod function
    int rlt = n%m;
    if (n < 0) rlt += m;
    return rlt;
}

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

void insertLoad_js(int idx, int& n) {
    L_js[idx] += 1;
    minbalance(idx, n);
}

double getRatio_js(int n) {
    double mn = (double)getMin_js(n);
    double mx = (double)getMax_js(n);
    if (mn <= 0.0000001) return mx;
    return mx/mn;
}

void print_js(int no, int n) {
    char str[100];
    sprintf(str, "%d.txt", no);
    ofstream out(str);
    for (int i  = 0; i < n; ++i) {
        out << i << ' ' << L_js[i] << '\n';
    }
    out.close();
    return;
}

int rt(int n) { // random [0, n-1]
    MTRand m;
    return m.randInt(n-1);
}

int main(){
    srand(time(0));

    int T = 1;
    while (T--) {

        int n;
        n = 2;

        L_js.clear();

        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            L_js.push_back(0);
        }

        char nameFileRatio[100];

        sprintf(nameFileRatio, "ratio(%d).txt", n);
        ofstream outRatio(nameFileRatio);

        outRatio << "0 " << getRatio_js(n) << '\n';

        // input data in max load in tank
        for (int i = 0; i < 100000; ++i) {
            int idx = rt(n); // random [0, n-1]
            insertLoad_js(idx, n);

            outRatio << i+1 << ' ' << getRatio_js(n) << '\n';
            // print_js(i, n);
        }
        outRatio.close();
    }
    printf("End Process");
    return 0;
}
