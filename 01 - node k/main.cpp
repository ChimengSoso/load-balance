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

multiset<int> load;
multiset<int>::iterator l_it; // load iterator

vector<int> tnk; // tank

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

int rt(int a, int b) { // random to
    return a + rand()%(b-a+1);
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
        j = rt(0, n-1);
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

int main(){
    srand(time(0));
    int T = 1;
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

        int constPoint = rt(0, n-1);

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
