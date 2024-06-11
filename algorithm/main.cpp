#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
using namespace std;

const double b = 0.1;
const double  c = 2.0;
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
bool IsReplace[30][30][33];
int Replace[30][30][33];
int ShiftReplace[30][30][33];
CSD csd[100][100];
CSD csd_matching[100][100];
vec pool_baseline[100000];
vec pool[100000];

int _string_mathcing(vector<int >& str, int i, int j, int l, int n, int num)
{
    int sum = 0;
    // 尝试匹配
    for (int u = i + 1; u < n; ++u) {
        // 从高位开始匹配
        for (int s = csd[u][j].len - l; s >= 0; --s) {
            bool IsMatching = true;
            for (int k = 0; k < l; ++k) {
                if (str[k] != csd_matching[u][j].data[s+k]) {
                    IsMatching = false;
                    break;
                }
            }
            if (IsMatching){
                ++sum;
                bool Rep = false;
                for (int k = 0; k < l; ++k) {
                    if (csd_matching[u][j].data[s+k] != 0) {
                        csd_matching[u][j].data[s+k] = 0;
                        IsReplace[u][j][s+k] = true;
                        if (Rep == false){
                            Replace[u][j][s+k] = num;
                            ShiftReplace[u][j][s+k] = s;
                            Rep = true;
                        }
                        else {
                            Replace[u][j][s+k] = -1;
                            ShiftReplace[u][j][s+k] = -1;
                        }
                    }
                }
            }
        }
    }
    return sum;
}

int string_matching(int n, int m, int l_min = 3, int l_max = 5)
{
    int num = 0;
    vector<int >str;
    vector<pair<vec*, pair<int, bool> > > source;
    // 枚举匹配字符串的长度
    for (int l = l_max; l >= l_min; --l) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                int sum = 0;
                // 枚举字符串起点
                for (int s = csd[i][j].len - l; s >= 0; --s) {
                    str.clear();
                    if (csd_matching[i][j].data[s + l - 1] == 0 || csd_matching[i][j].data[s] == 0) continue;
                    for (int k = 0; k < l; ++k)
                        str.push_back(csd_matching[i][j].data[s+k]);


                    sum = _string_mathcing(str, i, j, l, n, num);

                    if (sum > 0){
                        source.clear();
                        bool Rep = false;
                        for (int k = 0; k < l; ++k) {
                            if (str[k] != 0) {
                                if (str[k] == 1) {
                                    source.push_back(make_pair(&pool_baseline[j], make_pair(k, false)));
                                } else if (str[k] == -1) {
                                    source.push_back(make_pair(&pool_baseline[j], make_pair(k, true)));
                                }
                                csd_matching[i][j].data[s+k] = 0;
                                IsReplace[i][j][s+k] = true;
                                if (Rep == false){
                                    Replace[i][j][s+k] = num;
                                    ShiftReplace[i][j][s+k] = s;
                                    Rep = true;
                                }
                                else {
                                    Replace[i][j][s+k] = -1;
                                    ShiftReplace[i][j][s+k] = -1;
                                }
                            }
                        }
                        // 增加加法器
                        if (Rep) {
                            pool[num++] = (vec(source));
                        }
                    }
                }
            }
        }
    }

    return num;
}

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
            csd_matching[i][j].len = l;
            csd_matching[i][j].INT2CSD(a[i][j]);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            cout << csd[i][j] << " ";
        cout << endl;
    }

    int num_baseline = 0;
    for (int i = 0; i < m; i++) {
        pool_baseline[num_baseline++] = (vec(m, i, w));
    }
    for (int i = 0; i < n; i++) {
        vector<pair<vec*, pair<int, bool> > > source;
        for(int j = 0; j < m; j++) {
            for(int k = 0; k < csd[i][j].len; k++) {
                if (csd[i][j].data[k] == 1) {
                    source.push_back(make_pair(&pool_baseline[j], make_pair(k, false)));
                } else if (csd[i][j].data[k] == -1) {
                    source.push_back(make_pair(&pool_baseline[j], make_pair(k, true)));
                }
            }
        }
        pool_baseline[num_baseline++] = (vec(source));
    }

    double total_cost = 0;
    int total_cost_x = 0;
    int total_cost_b = 0;
    int total_cost_c = 0;
    for(int i = 0; i < num_baseline; i++) {
        cout <<"("<< pool_baseline[i] << ") width:" << pool_baseline[i].width << " cost:" << pool_baseline[i].cost() << endl;
        total_cost += pool_baseline[i].cost();
        total_cost_x += pool_baseline[i].cost_x;
        total_cost_b += pool_baseline[i].cost_b;
        total_cost_c += pool_baseline[i].cost_c;
    }
    printf("Total cost: %d + %dc + %db.\n", total_cost_x, total_cost_c, total_cost_b);
    printf("Total cost when c = %f, b = %f: %f\n", c, b, total_cost);
    
    printf("\n");

    int num = string_matching(n, m, 3, 3);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++)
            cout << csd_matching[i][j] << " ";
        cout << endl;
    }

    for (int i = 0; i < n; i++) {
        vector<pair<vec*, pair<int, bool> > > source;
        for(int j = 0; j < m; j++) {
            for(int k = 0; k < csd[i][j].len; k++) {
                if (IsReplace[i][j][k] == true) {
                    if (Replace[i][j][k] != -1){
                        source.push_back(make_pair(&pool[Replace[i][j][k]], make_pair(ShiftReplace[i][j][k], false)));
                    }
                    continue;
                }

                if (csd[i][j].data[k] == 1) {
                    source.push_back(make_pair(&pool_baseline[j], make_pair(k, false)));
                } else if (csd[i][j].data[k] == -1) {
                    source.push_back(make_pair(&pool_baseline[j], make_pair(k, true)));
                }
            }
        }
        pool[num++] = (vec(source));
    }

    total_cost = 0;
    total_cost_x = 0;
    total_cost_b = 0;
    total_cost_c = 0;
    for(int i = 0; i < num; i++) {
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

