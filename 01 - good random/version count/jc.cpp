#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "MersenneTwister.h"

using namespace std;
vector<long> L;
multiset<long> no;
multiset<long>::iterator it;
vector<long>::iterator Li;
//set ค่า alpha beta
float alpha = 5.464;
float beta = 4.1;
double c = 1;
double rohl = 2;

long xxx=0;
long yyy=0;

MTRand mtrand1;
long min_index_node = 0;
long no_insert=1000000;
long random_value[1000000];
double number_of_key_move=0;
float max_ratio=0;
float ratio[1000000];
long call_minbalance;

long findmax(int);
long findmin(int);
long indexmax(int);
long indexmin(int);
long indexminnot0(int p);
long findminnot0(int p);
void randomnumber(int);
int minbalance(long,int);
int split(long,int);
void printvalue(int p);
void printvaluetofile(int x,int p);
double findratio();
long findlevel(long);
void shuffle_randomvalue(int);
void add_hot_spot(double,int);
void add_hot_spots(double,int,int);
char* itoa( int value, char* result, int base ) ;
char buf[15];
char temp[10];
long load_start[1]= {0};
long no_node[4]= {10,1000,10000,100000};
void adjustload_insert(long index,int p);
void adjustload_delete(long index,int p);
void insertnewnode(int);
void deletenode(long index, int p);
void insertkey(long index,int p);
void deletekey(long index,int p);

int main() {
    ofstream out("No_key_move.txt");
    int p=0;
    //for(int p=0;p<4;p++) {
    for(int p=0; p<1; p++) {
        long i;
        //ตั้งชื่อไฟล์
        strcpy(buf,"");
        itoa(no_node[p], temp, 10);
        strcat(buf,"ratio(");
        strcat(buf,temp);
        strcat(buf,")");
        strcat(buf,".txt");
        ofstream out2(buf);
        //เคลียร์ค่า vector ที่แทน load L แทนโหลดของแต่ละโหนดเรียงกันไป ส่วน no เป็น load ที่เรียงงลำดับ เป็น multiset
        L.clear();
        no.clear();
        // init ค่า เพิ่มโหนดและให้ load เป็น 0
        for(i = min_index_node; i < no_node[p]; i++) {
            L.push_back(0);
            no.insert(0);
        }

        long h;
        h=no.size();
        cout<<h;
        number_of_key_move=0;
        long k;
        char x;
        int temp;
        out2<<"0 "<< findratio()<<endl;
        for( i=0; i<100; i++) {
            k=indexmax(p);
            insertkey(k,p);
            out2<<i+1<<" "<< findratio()<<endl;
            printvaluetofile(i,p);
        }
        deletenode(5,p);
        printvaluetofile(100,p);
        deletenode(8,p);
        printvaluetofile(101,p);

        /*
            for( i=10000;i<20000;i++) {
        //    call_minbalance[i]=0;
        //    k=indexmax(p);
              k=indexminnot0(p);
        //    it=no.end();
        //    it--;
              it=no.lower_bound(L[k]);
              no.erase(it);
              L[k]=L[k]-1;
              no.insert(L[k]);
              //adjustload_delete(k,p);
              call_minbalance=split(k,p);
              out2<<i+1<<" "<< findratio()<<endl;
            }
         */
        //  cout<<"adjustload insert = "<< xxx <<endl;
        //  cout<<"adjustload delete = "<< yyy <<endl;
//   xxx=0;    yyy=0;
        out2.close();
        out<<no_node[p]<<" "<< number_of_key_move<<endl;
    }
    out.close();
    return 0;
}
long findmax(int p) {
    it=no.end();
    it--;
    return *it;
}
long indexmax(int p) {
    long answer,index,i;
    answer = L[min_index_node];
    index = min_index_node;
    for(i = min_index_node; i < no_node[p]; i++ ) {
        if(L[i] > answer) {
            answer = L[i];
            index = i;
        }
    }
    return index;
}
long findmin(int p) {
    it=no.begin();
    return *it;
}
long indexmin(int p) {
    long answer,index,i;
    answer = L[min_index_node];
    index = min_index_node;
    for(i = min_index_node; i < no_node[p]; i++ ) {
        if(L[i] < answer) {
            answer = L[i];
            index = i;
        }
    }
    return index;
}
int minbalance(long index,int p) {
    int call=0;
    if(L[index] > (double)(alpha*(double)findmin(p))) {
        long index_minimum;
        index_minimum = indexmin(p);
        //เราจะหาโหนดที่ติดกับโหนดที่น้อยที่สุดที่มี load น้อยสุด แล้วให้รับโหลดจากโหนดน้อยสุด
        if(index_minimum == min_index_node) { //ถ้าเป็นขอบซ้าย
            if(L[min_index_node + 1] > L[no_node[p]-1]) {
                it=no.lower_bound(L[no_node[p]-1]);
                no.erase(it);
                no.insert(L[(no_node[p]-1)] + L[min_index_node]);
                L[(no_node[p]-1)] = L[(no_node[p]-1)] + L[min_index_node];
            } else {
                it=no.lower_bound(L[min_index_node + 1]);
                no.erase(it);
                no.insert(L[min_index_node + 1] + L[min_index_node]);
                L[min_index_node + 1] = L[min_index_node + 1]+L[min_index_node];
            }
        } else if(index_minimum == (no_node[p]-1)) { //ถ้าเป็นขอบขวา
            if(L[min_index_node] > L[(no_node[p]-1) - 1]) {
                it=no.lower_bound(L[(no_node[p]-1) - 1]);
                no.erase(it);
                no.insert(L[(no_node[p]-1) - 1] + L[(no_node[p]-1)]);
                L[(no_node[p]-1) - 1] = L[(no_node[p]-1) - 1] + L[(no_node[p]-1)];
            } else {
                it=no.lower_bound(L[min_index_node]);
                no.erase(it);
                no.insert(L[min_index_node] + L[(no_node[p]-1)]);
                L[min_index_node] = L[min_index_node] + L[(no_node[p]-1)];
            }
        } else {
            if(L[index_minimum - 1] > L[index_minimum + 1]) {
                it=no.lower_bound(L[index_minimum + 1]);
                no.erase(it);
                no.insert(L[index_minimum + 1] + L[index_minimum]);
                L[index_minimum + 1] = L[index_minimum + 1] + L[index_minimum];
            } else {
                it=no.lower_bound(L[index_minimum - 1]);
                no.erase(it);
                no.insert(L[index_minimum - 1] + L[index_minimum]);
                L[index_minimum - 1] = L[index_minimum - 1] + L[index_minimum];
            }
        }//จัดการกับโหนดที่โหลดน้อยสุด
        L[index_minimum] = 0;
        it=no.begin();
        no.erase(it);
        L.erase(L.begin() + index_minimum );
        if(index_minimum < index) {
            index = index - 1;
        }
        long m, n;
        m = (long)floor(L[index] * 0.5);
        n = (long)ceil(L[index] * 0.5);
        it=no.lower_bound(L[index]);
        no.erase(it);
        L[index] = m;
        L.insert(L.begin() + index, n );
        no.insert(m);
        no.insert(n);
        call=1;
        number_of_key_move=number_of_key_move+n;
    }
    xxx++;
    return call;
}
long findlevel(long load) {
    double level=0;
    while(load>(long)pow(rohl,level)) {
        level++;
    }
    return (long)level;
}
void adjustload_insert(long index,int p) {
    int call=0;
    long level = findlevel(L[index]);
    long j=0;
    xxx++;
    if(index==min_index_node) {
        if(L[index+1]>L[no_node[p]-1])
            j=no_node[p]-1;
        else
            j=index+1;
    } else if(index==no_node[p]-1) {
        if(L[min_index_node]>L[index-1])
            j=index-1;
        else
            j=min_index_node;
    } else {
        if(L[index+1]>L[index-1])
            j=index-1;
        else
            j=index+1;
    }
    if(findlevel(L[j])<=level-2) {//NBRADJUST
        long x,y;
        x = (long)ceil((L[j]+L[index])* 0.5);
        y = (long)floor((L[j]+L[index])* 0.5);
        it=no.lower_bound(L[index]);
        no.erase(it);
        it=no.lower_bound(L[j]);
        no.erase(it);
        L[j] = x;
        L[index] = y;
        no.insert(x);
        no.insert(y);
        adjustload_insert(j,p);
        adjustload_insert(index,p);
        //xxx++;
    } else {
        long index_minimum;
        index_minimum = indexmin(p);
        long N;
        if(findlevel(L[index_minimum])<=level-3) {//REORDER
            if(index_minimum == min_index_node) {
                if(L[min_index_node + 1] > L[no_node[p]-1])
                    N=no_node[p]-1;
                else
                    N=min_index_node + 1;
            } else if(index_minimum == (no_node[p]-1)) {
                if(L[min_index_node] > L[(no_node[p]-1) - 1])
                    N=(no_node[p]-1) - 1;
                else
                    N=min_index_node;
            } else {
                if(L[index_minimum - 1] > L[index_minimum + 1])
                    N=index_minimum + 1;
                else
                    N=index_minimum - 1;
            }
            it=no.begin();
            no.erase(it);
            it=no.lower_bound(L[N]);
            no.erase(it);
            no.insert(L[N] + L[index_minimum]);
            L[N] = L[N] + L[index_minimum];
            L.erase(L.begin() + index_minimum );
            if(index_minimum < index) {
                index = index - 1;
            }
            if(index_minimum < N) {
                N = N - 1;
            }
            long x, y;
            x = (long)floor(L[index] * 0.5);
            y = (long)ceil(L[index] * 0.5);
            it=no.lower_bound(L[index]);
            no.erase(it);
            L[index] = x;
            L.insert(L.begin() + index, y );
            no.insert(x);
            no.insert(y);
            adjustload_insert(N,p);
            //xxx++;
        }
    }
}
void adjustload_delete(long index,int p) {
    int call=0;
    long level = findlevel(L[index]);
    long j=0;
    yyy++;
    if(index==min_index_node) {
        if(L[index+1]>L[no_node[p]-1])
            j=no_node[p]-1;
        else
            j=index+1;
    } else if(index==no_node[p]-1) {
        if(L[min_index_node]>L[index-1])
            j=index-1;
        else
            j=min_index_node;
    } else {
        if(L[index+1]>L[index-1])
            j=index-1;
        else
            j=index+1;
    }
    if(findlevel(L[j])>=level+2) {//NBRADJUST
        long x,y;
        x = (long)ceil((L[j]+L[index])* 0.5);
        y = (long)floor((L[j]+L[index])* 0.5);
        it=no.lower_bound(L[index]);
        no.erase(it);
        it=no.lower_bound(L[j]);
        no.erase(it);
        L[j] = y;
        L[index] = x;
        no.insert(x);
        no.insert(y);
        adjustload_delete(j,p);
        adjustload_delete(index,p);
    } else {
        long index_maximum;
        index_maximum = indexmax(p);
        long N;
        if(findlevel(L[index_maximum])>=level+3) {//REORDER
            it=no.lower_bound(L[index]);
            no.erase(it);
            it=no.lower_bound(L[j]);
            no.erase(it);
            no.insert(L[index] + L[j]);
            L[j] = L[j] + L[index];
            L.erase(L.begin() + index );
            if(index_maximum > index) {
                index_maximum = index_maximum - 1;
            }
            long x, y;
            x = (long)floor(L[index_maximum] * 0.5);
            y = (long)ceil(L[index_maximum] * 0.5);
            it=no.lower_bound(L[index_maximum]);
            no.erase(it);
            L[index_maximum] = x;
            L.insert(L.begin() + index_maximum, y );
            no.insert(x);
            no.insert(y);
            adjustload_delete(j,p);
        }
    }
}
void printvalue(int p) {
    long i;
    for(i = min_index_node; i < no_node[p]; i++ ) {
        cout << L[i] <<" ";
    }
}
void printvaluetofile(int x,int p) {
    strcpy(buf,"");
    itoa(x, temp, 10);
    strcat(buf,temp);
    strcat(buf,".txt");
    ofstream out3(buf);
    long i;
    for(i = min_index_node; i < no_node[p]; i++ ) {
        out3<< i<<" "<<L[i] <<endl;
    }
    out3.close();
}
int split(long index,int p) {
    yyy++;
    long maxvalue=findmax(p);
    if(L[index] <= (double)((double)maxvalue/beta)) {
        long index_neighbor;//หา neighbor ของ u ที่ load น้อยสุด
        if(index == min_index_node) {
            if(L[index+1] > L[no_node[p]-1]) {
                index_neighbor = (no_node[p]-1);
            } else {
                index_neighbor = index+1;
            }
        } else if(index == (no_node[p]-1)) {
            if(L[min_index_node] > L[index-1]) {
                index_neighbor = index-1;
            } else {
                index_neighbor = min_index_node;
            }
        } else {
            if(L[index+1] > L[index-1]) {
                index_neighbor = index-1;
            } else {
                index_neighbor = index+1;
            }
        }
        long index_max,loadmax;
        index_max = indexmax(p);
        if(L[index_neighbor] > (double)((double)2.0*maxvalue/beta)) {
            //SPLITNEIGHBOR
            long m,n;
            m = (long)ceil((L[index_neighbor]+L[index])* 0.5);
            n = (long)floor((L[index_neighbor]+L[index])* 0.5);
            number_of_key_move=number_of_key_move+(L[index_neighbor]-m);
            it=no.lower_bound(L[index]);
            no.erase(it);
            it=no.lower_bound(L[index_neighbor]);
            no.erase(it);
            L[index_neighbor] = m;
            L[index] = n;
            no.insert(m);
            no.insert(n);
        } else {
            //SPLITMAX
            it=no.lower_bound(L[index_neighbor]);
            no.erase(it);
            it=no.lower_bound(L[index]);
            no.erase(it);

            no.insert(L[index_neighbor] + L[index]);
            L[index_neighbor] = L[index_neighbor] + L[index];
            L.erase(L.begin() + index );
            long m,n;
            m = (long)ceil(maxvalue* 0.5);
            n = (long)floor(maxvalue* 0.5);
            if(index < index_max) {
                index_max = index_max-1;
            }
            L[index_max] = m;
            it=no.end();
            it--;
            no.erase(it);
            L.insert(L.begin() + index_max, n );
            no.insert(m);
            no.insert(n);
            number_of_key_move=number_of_key_move+n;
        }
    }
}
double findratio() {
    double answer;
    long max,min;
    it=no.end();
    it--;
    max=*it;
    it=no.begin();
    min=*it;
    if(min==0) {
        answer=max;
    } else {
        answer = (double)((double)max/(double)min);
    }
    return answer;
}
void randomnumber(int p) {
    MTRand mtrand1;
    int i;
    for(i=0; i<no_insert; i++) {
        random_value[i] = 0;
    }
    long r;
    for(i=0; i<no_insert; i++) {
        r=mtrand1.randInt( no_node[p]-1);
        random_value[i] = r;
    }
}
void shuffle_randomvalue() {
    random_shuffle(random_value, random_value + (no_insert));
}
void add_hot_spot(double percent,int p) {
    long amount=(long)ceil(no_insert*(double)(percent/100.0));
    MTRand mtrand1;
    long r=mtrand1.randInt( (no_node[p]-1) );
    for(long i=0; i<amount; i++) {
        random_value[i]=r;
    }
    shuffle_randomvalue();
}
void add_hot_spots(double percent,int n,int p) {
    long amount=(long)ceil(no_insert*(double)(percent/100.0));
    MTRand mtrand1;
    long i;
    for(long j=0; j<n; j++) {
        long hotspot_id= mtrand1.randInt( (no_node[p]-1) );
        for(i=(amount*j); i<(amount*(j+1)); i++) {
            random_value[i]=hotspot_id;
        }
    }
    shuffle_randomvalue();
}
char* itoa( int value, char* result, int base ) {
    // check that the base if valid
    if (base < 2 || base > 16) {
        *result = 0;
        return result;
    }
    char* out = result;
    int quotient = value;
    do {
        *out = "0123456789abcdef"[ std::abs( quotient % base ) ];
        ++out;
        quotient /= base;
    } while ( quotient );
    // Only apply negative sign for base 10
    if ( value < 0 && base == 10) *out++ = '-';
    std::reverse( result, out );
    *out = 0;
    return result;
}
long indexminnot0(int p) {
    long answer,index,i;
    answer = findmax(p);
    index = indexmax(p);
    for(i = min_index_node; i < no_node[p]; i++ ) {
        if(L[i] < answer&&L[i]>0) {
            answer = L[i];
            index = i;
        }
    }
    return index;
}
long findminnot0(int p) {
    it=no.lower_bound(1);
    return *it;
}
void insertnewnode(int p) {
    long maxvalue=findmax(p);
    long index_max,loadmax;
    index_max = indexmax(p);
    long m,n;
    m = (long)ceil(maxvalue* 0.5);
    n = (long)floor(maxvalue* 0.5);
    number_of_key_move=number_of_key_move+n;
    L[index_max] = m;
    it=no.end();
    it--;
    no.erase(it);
    L.insert(L.begin() + index_max, n );
    no.insert(m);
    no.insert(n);
    no_node[p]=no_node[p]+1;
}

void deletenode(long index, int p) {
    //ตรวจดูเพื่อนบ้าน ถ้าเพื่อนบ้านรวมกับโหนดที่จะลบแล้วไม่เป็นไร ก็โอนให้เพื่อนบ้าน
    //ถ้าเพื่อนบ้านเยอะกันหมด ไปหาโหนด min แล้ว  insert ให้โหนดเพื่อนบ้าน แล้วไปแทนที่โหนดที่จะลบ
    long index_neighbor,index_minimum,min_neighbor;
    if(index == min_index_node) {
        if(L[index+1] > L[no_node[p]-1]) {
            index_neighbor = (no_node[p]-1);
        } else {
            index_neighbor = index+1;
        }
    } else if(index == (no_node[p]-1)) {
        if(L[min_index_node] > L[index-1]) {
            index_neighbor = index-1;
        } else {
            index_neighbor = min_index_node;
        }
    } else {
        if(L[index+1] > L[index-1]) {
            index_neighbor = index-1;
        } else {
            index_neighbor = index+1;
        }
    }
    long maxvalue = findmax(p);
    if(maxvalue>L[index_neighbor]+L[index]) {
        it=no.lower_bound(L[index]);
        no.erase(it);
        it=no.lower_bound(L[index_neighbor]);
        no.erase(it);
        L[index_neighbor]=L[index_neighbor]+L[index];
        no.insert(L[index_neighbor]);
        L.erase(L.begin() + index );
        no_node[p]=no_node[p]-1;
    } else {

        index_minimum = indexmin(p);
        if(index_minimum == min_index_node) {
            if(L[index_minimum+1] > L[no_node[p]-1]) {
                min_neighbor = (no_node[p]-1);
            } else {
                min_neighbor = index_minimum+1;
            }
        } else if(index_minimum == (no_node[p]-1)) {
            if(L[min_index_node] > L[index_minimum-1]) {
                min_neighbor = index_minimum-1;
            } else {
                min_neighbor = min_index_node;
            }
        } else {
            if(L[index_minimum+1] > L[index_minimum-1]) {
                min_neighbor = index_minimum-1;
            } else {
                min_neighbor = index_minimum+1;
            }
        }

        long temp=L[index_minimum];
        L.erase(L.begin() + index_minimum );
        it=no.lower_bound(L[index_minimum]);
        no.erase(it);
        if(min_neighbor>index_minimum) {
            min_neighbor=min_neighbor-1;
        }
        no_node[p]=no_node[p]-1;
        for(int i=0; i<temp; i++) {
            insertkey(min_neighbor,p);
        }
    }

}
void insertkey(long index,int p) {
    it=no.lower_bound(L[index]);
    no.erase(it);
    L[index]=L[index]+1;
    no.insert(L[index]);
    call_minbalance=minbalance(index,p);
}
void deletekey(long index,int p) {
    it=no.lower_bound(L[index]);
    no.erase(it);
    L[index]=L[index]-1;
    no.insert(L[index]);
    call_minbalance=split(index,p);
}
