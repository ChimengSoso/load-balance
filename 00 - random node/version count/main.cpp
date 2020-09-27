#include <stdio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <set>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
using namespace std;

const int N = 1000;

vector<int> tnk; // tank

int cost;        // use count cost of all operation in algorithm
int moveCnt;     // use count only move operation
int addCnt;      // use count only add operation
int deleteCnt;   // use count only delete operation

const int moveCst   = 1;      // move cost
const int addCst    = 1;      // add cost
const int deleteCst = 1;      // delete cost

const double alpha = 5.464;
const double beta = 4.1;

void setCnt(bool all = 1, bool move = 1, bool add = 1, bool deleteC = 1) {
    if (all == true) cost = 0;
    if (move == true) moveCnt = 0;
    if (add == true) addCnt = 0;
    if (deleteC == true) deleteCnt = 0;
}

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

void minbalance(int u, int& n) {
    int v = getIdxMin(n);
    if (tnk[u] > alpha*tnk[v]) {
        int l = md(v-1, n);
        int r = md(v+1, n);
        int z = tnk[l] < tnk[r] ? l : r;

        tnk[z]  += tnk[v];
        moveCnt += tnk[v]*moveCst;      // count operation move

        #define start tnk.begin()

        tnk.erase(start+v);

        if (v < u)
            --u;

        int remain_up = (tnk[u]+1)/2;
        int remain_down =   tnk[u]/2;
        tnk[u] = remain_up;

        moveCnt += remain_up*moveCst;    // count operation move

        tnk.insert(start+u, remain_down);
    }
    return;
}

void insertLoad(int idx, int& n) {
    tnk[idx] += 1;
    addCnt += addCst;

    minbalance(idx, n);
}

void split(int u, int n) {
    int w = getIdxMax(n);
    if (tnk[u] < tnk[w]/beta) {
        int l = md(u-1, n);
        int r = md(u+1, n);
        int z = tnk[l] < tnk[r] ? l : r;

        if (tnk[z] <= 2*tnk[w]/beta) {

            tnk[z] += tnk[u];
            tnk.erase(tnk.begin() + u);

            if (u < w)
                --w;

            int remain_up = (tnk[w]+1)/2;
            int remain_down = tnk[w] / 2;

            tnk[w] = remain_up;
            tnk.insert(tnk.begin() + w, remain_down);

        }else{
            int sum = tnk[u] + tnk[z];

            tnk[u] = (sum+1)/2;
            tnk[z] = sum /2;
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

int rt(int a, int b) { // random to
    return a + rand()%(b-a+1);
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

        char nameFileRatio[100];
        char nameFileCount[100];

        sprintf(nameFileRatio, "ratio(%d).txt", n);
        sprintf(nameFileCount, "count(%d).txt", n);
        ofstream outRatio(nameFileRatio);
        ofstream outCount(nameFileCount);

        outRatio << "0 " << getRatio(n) << '\n';

        setCnt();
        long long total = 0LL;
        // input data in max load in tank
        for (int i = 0; i < 100000; ++i) {
            setCnt(false);
            int idx = rt(0, n-1);
            insertLoad(idx, n);

            outRatio << i+1 << ' ' << getRatio(n) << '\n';

            cost = moveCnt + addCnt + deleteCnt;
            total += cost;
            outCount << i+1 << ' ' << cost << '\n';
            //print(i, n);
        }
        outCount << "total = " << total << '\n';
        outCount.close();
        outRatio.close();
    }
    printf("End process\n");
    return 0;
}
