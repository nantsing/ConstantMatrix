#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
using namespace std;

const double b = 0.1;
const double  c = 1.0;
const int w = 8;

struct CSD
{
    int data[32];
    int len;
    CSD(int len = 32) : len(len) {
        for (int i = 0; i < 32; i++)
            data[i] = 0;
    }
    void INT2CSD(int a){
        bool neg = false;
        if (a < 0) {
            a = -a;
            neg = true;
        }
        for (int k = 0; k < len; k++) {
            data[k] = (a & (1 << k)) ? 1 : 0;
        }
        for (int i = 0; i < len - 1; i++) {
            if (data[i] == 1 && data[i + 1] == 1) {
                data[i] = -1;
                data[i + 1] = 0;
                data[i + 2] += 1;
                if (i == len - 2)
                    len++;
            }
            int k = i + 2;
            if (k >= len) continue;
            while (data[k] == 2) {
                data[k] = 0;
                data[k + 1] += 1;
                k++;
            }
            len = max(len, k + 1);
        }
        if (neg) {
            for (int i = 0; i < len; i++) {
                data[i] = -data[i];
            }
        }
    }
    friend ostream& operator<<(ostream& os, const CSD& csd) {
        for (int i = csd.len - 1; i >= 0; i--)
            os << csd.data[i];
        return os;
    }
};



struct vec{
    vector<int> data;
    vector<pair<vec* ,pair<int, bool> > > source;
    bool isRoot;
    int width;
    int cost_x, cost_c, cost_b;
    vec() {
        isRoot = false;
        width = 0;
    }
    vec(int n, int k, int w) {
        data.resize(n);
        data[k] = 1;
        isRoot = true;
        width = w;
    }
    vec(vector<pair<vec* ,pair<int, bool> > > source) {
        this->source = source;
        isRoot = false;
        vec tmp;
        tmp.data.resize(source[0].first->data.size());
        for (int i = 0; i < source.size(); i++) {
            int wi = source[i].first->width + source[i].second.first;
            if (source[i].second.second) {
                tmp = tmp - (*(source[i].first) << source[i].second.first);
            } else {
                tmp = tmp + (*(source[i].first) << source[i].second.first);
            }
        }
        data = tmp.data;
        width = 0;
        for (int i = 0; i < data.size(); i++) {
            width += ((1 << w) - 1) * abs(data[i]);
        }
        width = log2(width) + 1;
    }
    double cost() {
        double cost = 0;
        int max_width = 0;
        int min_width = 100000;
        cost_x = 0;
        cost_c = 0;
        cost_b = 0;
        if (isRoot) {
            return 0;
        }
        for (int i = 0; i < source.size(); i++) {
            int wi = source[i].first->width + source[i].second.first;
            if (wi > max_width) {
                max_width = wi;
            }
            if (wi < min_width) {
                min_width = wi;
            }
            cost += wi;
            cost_x += wi;
            if (source[i].second.second) {
                cost += b * wi;
                cost_b += wi;
            }
        }
        cost += c * max_width;
        cost_c = max_width;
        cost -= min_width;
        cost_x -= min_width;
        return cost;
    }
    friend ostream& operator<<(ostream& os, const vec& v) {
        for (int i = 0; i < v.data.size(); i++)
            os << v.data[i] << " ";
        return os;
    }
    vec operator - () {
        vec c;
        for (int i = 0; i < data.size(); i++)
            c.data.push_back(-data[i]);
        return c;
    }
    vec operator << (int k) {
        vec c;
        for (int i = 0; i < data.size(); i++)
            c.data.push_back(data[i]<<k);
        return c;
    }
    vec operator + (const vec& v) {
        vec c;
        assert(data.size() == v.data.size());
        for (int i = 0; i < data.size(); i++)
            c.data.push_back(data[i] + v.data[i]);
        return c;
    }
    vec operator - (const vec& v) {
        vec c;
        assert(data.size() == v.data.size());
        for (int i = 0; i < data.size(); i++)
            c.data.push_back(data[i] - v.data[i]);
        return c;
    }
};


int a[100][100];
CSD csd[100][100];
int main() {
    
    int n, m;
    cin >> n >> m;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            cin >> a[i][j];
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            int l = log2(abs(a[i][j])) + 1;
            l = max(l, 1);
            csd[i][j].len = l;
            csd[i][j].INT2CSD(a[i][j]);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            cout << csd[i][j] << " ";
        cout << endl;
    }

    vector<vec> pool;
    for (int i = 0; i < m; i++) {
        pool.push_back(vec(m, i, w));
    }
    for (int i = 0; i < n; i++) {
        vector<pair<vec*, pair<int, bool> > > source;
        for(int j = 0; j < m; j++) {
            for(int k = 0; k < csd[i][j].len; k++) {
                if (csd[i][j].data[k] == 1) {
                    source.push_back(make_pair(&pool[j], make_pair(k, false)));
                } else if (csd[i][j].data[k] == -1) {
                    source.push_back(make_pair(&pool[j], make_pair(k, true)));
                }
            }
        }
        pool.push_back(vec(source));
    }
    
    double total_cost = 0;
    int total_cost_x = 0;
    int total_cost_b = 0;
    int total_cost_c = 0;
    for(int i = 0; i < pool.size(); i++) {
        cout <<"("<< pool[i] << ") width:" << pool[i].width << " cost:" << pool[i].cost() << endl;
        total_cost += pool[i].cost();
        total_cost_x += pool[i].cost_x;
        total_cost_b += pool[i].cost_b;
        total_cost_c += pool[i].cost_c;
    }
    printf("Total cost: %d + %dc + %db.\n", total_cost_x, total_cost_c, total_cost_b);
    printf("Total cost when c = %f, b = %f: %f\n", c, b, total_cost);
    return 0;
}

