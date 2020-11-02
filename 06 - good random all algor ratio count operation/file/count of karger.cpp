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

vector<int> L_kg; // tank

const double eps = 0.023;

int cost_kg;        // use count cost of all operation in algorithm
int moveCnt_kg;     // use count only move operation
int addCnt_kg;      // use count only add operation
int deleteCnt_kg;   // use count only delete operation

const int moveCst_kg   = 1;      // move cost
const int addCst_kg    = 1;      // add cost
const int deleteCst_kg = 1;      // delete cost

const double alpha = 5.464;
const double beta = 4.1;

void setCnt_kg(bool all = 1, bool move = 1, bool add = 1, bool deleteC = 1) {
    if (all == true) cost_kg = 0;
    if (move == true) moveCnt_kg = 0;
    if (add == true) addCnt_kg = 0;
    if (deleteC == true) deleteCnt_kg = 0;
}

int getIdxMax(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_kg[i] > L_kg[idx])
            idx = i;
    }
    return idx;
}

int getIdxMin(int n) {
    int idx = 0;
    for (int i = 1; i < n; ++i) {
        if (L_kg[i] < L_kg[idx])
            idx = i;
    }
    return idx;
}

int getMax(int n) {
    return L_kg[getIdxMax(n)];
}

int getMin(int n) {
    return L_kg[getIdxMin(n)];
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
        return -1;
    }
}

void divideEq(int a, int b) {
    moveCnt_kg += ( abs(L_kg[a]-L_kg[b])/2 ) *moveCst_kg;             // count operation move
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

                moveCnt_kg += L_kg[j] *moveCst_kg;                // count operation move
                L_kg[j_right] += L_kg[j];
                
                L_kg.erase(L_kg.begin()+j);

                if (j < i)
                    --i;

                moveCnt_kg += m * moveCst_kg;                    // count operation add
                L_kg.erase(L_kg.begin()+i);
                L_kg.insert(L_kg.begin()+i, m);
                L_kg.insert(L_kg.begin()+i, n);
            }
        }
    }
    return;
}

void insertLoad(int idx, int& n) {
    L_kg[idx] += 1;
    addCnt_kg += 1 *addCst_kg;                                   // count operation add
    itemBalance(idx, n);
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
        n = 4;
        //scanf("%d", &n);

        L_kg.clear();

        // assign initial data in tank
        for (int i = 0; i< n; ++i) {
            L_kg.push_back(0);
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

        setCnt_kg();
        long long total_kg = 0LL;
        // input data in max load in tank
        for (int i = 0; i < 100; ++i) {
            setCnt_kg(false);
            int idx = rt(n); // random [0, n-1]
            insertLoad(idx, n);

            outRatio << i+1 << ' ' << getRatio(n) << '\n';

            cost_kg = moveCnt_kg + addCnt_kg + deleteCnt_kg;
            outCount << i+1 << ' ' << cost_kg << '\n';

            total_kg += cost_kg;
            //print(i, n);
        }
        outCount << "total_kg = " << total_kg << '\n';
        outCount.close();
        outRatio.close();
    }
    printf("End process\n");
    return 0;
}
