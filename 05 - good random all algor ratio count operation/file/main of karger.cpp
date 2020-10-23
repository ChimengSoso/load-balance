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
#include <assert.h>
using namespace std;

vector<int> L_kg; // tank of Karger => L_kg

const double alpha = 5.464;
const double beta = 4.1;
const double eps = 0.023;

int getIdxMax_kg(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_kg[i] > L_kg[idx])
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

int getMax_kg(int n) {
    return L_kg[getIdxMax_kg(n)];
}

int getMin_kg(int n) {
    return L_kg[getIdxMin_kg(n)];
}

int md(int n, int m) { // mod function
    int rlt = n%m;
    if (n < 0) rlt += m;
    return rlt;
}

int rt(int n) { // random [0, n-1]
    MTRand m;
    int idx = m.randInt(n-1);
    if (idx <= n-1 && 0 <= idx) return idx;
    else {
        printf("Error index\n");
        assert(idx <= n-1 && 0 <= idx);
    }
}

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

void insertLoad_kg(int idx, int& n) {
    L_kg[idx] += 1;
    itemBalance(idx, n);
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

int main(){
    srand(time(0));
    int T = 1;
    while (T--) {
        int n;
        n = 1024;
        //scanf("%d", &n);

        L_kg.clear();
        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            L_kg.push_back(0);
        }

        char nameFileRatio[100];

        sprintf(nameFileRatio, "ratio(%d).txt", n);
        ofstream outRatio(nameFileRatio);

        outRatio << "0 " << getRatio_kg(n) << '\n';

        // input data in max load in tank
        for (int i = 0; i < 100000; ++i) {
            int idx = rt(n); // random [0, n-1]
            insertLoad_kg(idx, n);
            outRatio << i+1 << ' ' << getRatio_kg(n) << '\n';
            // print_kg(i, n);
        }
        outRatio.close();
    }
    printf("End Process");
    return 0;
}
