#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>
#include <random>
#include <string.h>
using namespace std;

const double b = 0.0;
const double  c = 1.2;
const int w = 8;
const int ROUND = 1000;

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
    CSD& operator=(const CSD& other) {
        if (this == &other) {
            return *this;
        }
        len = other.len;
        for (int i = 0; i < len; i++) {
            data[i] = other.data[i];
        }
        return *this;
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
    vec& operator=(const vec& other) {
        if (this == &other) {
            return *this;
        }
        data = other.data;
        source = other.source;
        isRoot = other.isRoot;
        width = other.width;
        cost_x = other.cost_x;
        cost_c = other.cost_c;
        cost_b = other.cost_b;
        return *this;
    }
};

int a[100][100];
bool IsReplace[30][30][33];
bool NegReplace[30][30][33];
int Replace[30][30][33];
int ShiftReplace[30][30][33];
CSD csd[100][100];
CSD csd_matching[100][100];
vec pool_baseline[100000];
vec pool[100000];
vec pool_tmp[100000];
int current_cost_x, current_cost_c, current_cost_b;

bool _local_search(int n, int m, int num, double p, mt19937 &gen)
{
    // 检查是否是赚的，若是，则以p的概率返回true，否则以1-p的概率返回true

    // return true;
    for (int i = 0; i < n; i++) {
        vector<pair<vec*, pair<int, bool> > > source;
        for(int j = 0; j < m; j++) {
            for(int k = 0; k < csd[i][j].len; k++) {
                if (IsReplace[i][j][k] == true) {
                    if (Replace[i][j][k] != -1){
                        source.push_back(make_pair(&pool_tmp[Replace[i][j][k]], \
                             make_pair(ShiftReplace[i][j][k], NegReplace[i][j][k])));
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
        pool_tmp[num++] = (vec(source));
    }

    double total_cost = 0;
    int total_cost_x = 0;
    int total_cost_b = 0;
    int total_cost_c = 0;
    for(int i = 0; i < num; i++) {
        total_cost += pool_tmp[i].cost();
        total_cost_x += pool_tmp[i].cost_x;
        total_cost_b += pool_tmp[i].cost_b;
        total_cost_c += pool_tmp[i].cost_c;
    }

    double current_cost1 = (double)current_cost_x + (double)current_cost_c * c + (double)current_cost_b * b;
    double cost1 = (double)total_cost_x + (double)total_cost_c * c + (double)total_cost_b * b;
    // double current_cost2 = (double)current_cost_x + (double)current_cost_c * 2.0 + (double)current_cost_b * b;
    // double cost2 = (double)total_cost_x + (double)total_cost_c * 2.0 + (double)total_cost_b * b;

    // if ((cost1 < current_cost1 && cost2 <= current_cost2) || \
    //         (cost1 <= current_cost1 && cost2 < current_cost2)){
    //             current_cost_x = total_cost_x;
    //             current_cost_b = total_cost_b;
    //             current_cost_c = total_cost_c;
    //             return true;
    //         }
    
    uniform_real_distribution<double> distr(0, 1);
    if ((cost1 < current_cost1 && distr(gen) <= p) ||
        (cost1 >= current_cost1 && distr(gen) > p)){
        current_cost_x = total_cost_x;
        current_cost_b = total_cost_b;
        current_cost_c = total_cost_c;
        return true;
    }

    return false;
}

// 检查csd[u][j]的[s:s+l]在乘上符号后是否与str相同
bool _Is_matching(vector<int >& str, int l, int u, int j, int s, int sig)
{
    for (int k = 0; k < l; ++k) {
        if (str[k] != sig * csd_matching[u][j].data[s+k]) {
            return false;
            break;
        }
    }
    return true;
}

void _recover(vector<int >& str, int i, int j, int l, int n)
{
    for (int u = i; u < n; ++u) {
        for (int s = csd[u][j].len - l; s >= 0; --s) {
            bool IsMatching = _Is_matching(str, l, u, j, s, 1);
            if (IsMatching){
                for (int k = 0; k < l; ++k) 
                    if (csd_matching[u][j].data[s+k] != 0) 
                        IsReplace[u][j][s+k] = false;
            continue;
            }
            
            IsMatching = _Is_matching(str, l, u, j, s, -1);
            if (IsMatching){
                for (int k = 0; k < l; ++k) 
                    if (csd_matching[u][j].data[s+k] != 0) 
                        IsReplace[u][j][s+k] = false;
            }

        }
    }

}

void _set0(int l, int u, int j, int s)
{
    for (int k = 0; k < l; ++k) {
        if (csd_matching[u][j].data[s+k] != 0) 
            csd_matching[u][j].data[s+k] = 0;
    }
}

void _string_matching_set0(vector<int >& str, int i, int j, int l, int n)
{
    for (int u = i; u < n; ++u) {
        for (int s = csd[u][j].len - l; s >= 0; --s) {
            bool IsMatching = _Is_matching(str, l, u, j, s, 1);
            if (IsMatching){
                _set0(l, u, j, s);
                continue;
            }

            // 匹配负数
            IsMatching = _Is_matching(str, l, u, j, s, -1);
            if (IsMatching){
                _set0(l, u, j, s);
            }

        }
    }
}

void _modify_matching(int l, int u, int j, int s, int num, bool neg)
{
    bool Rep = false;
    for (int k = 0; k < l; ++k) {
        if (csd_matching[u][j].data[s+k] != 0) {
            IsReplace[u][j][s+k] = true;
            if (Rep == false){
                Replace[u][j][s+k] = num;
                ShiftReplace[u][j][s+k] = s;
                NegReplace[u][j][s+k] = neg;
                Rep = true;
            }
            else {
                Replace[u][j][s+k] = -1;
                ShiftReplace[u][j][s+k] = -1;
            }
        }
    }
}

int _string_matching(vector<int >& str, int i, int j, int l, int n, int num)
{
    int sum = 0;
    // 尝试匹配
    for (int u = i; u < n; ++u) {
        // 从高位开始匹配
        for (int s = csd[u][j].len - l; s >= 0; --s) {
            bool IsMatching = _Is_matching(str, l, u, j, s, 1);
            if (IsMatching){
                ++sum;
                _modify_matching(l, u, j, s, num, false);
                s = s - l;
                continue;
            }
            
            // 匹配负数
            IsMatching = _Is_matching(str, l, u, j, s, -1);
            if (IsMatching){
                _modify_matching(l, u, j, s, num, true);
                s = s - l;
                ++sum;
            }
        }
    }
    return sum;
}

int string_matching(int n, int m, int l_min = 3, int l_max = 5, double p = 1.0)
{
    int num = 0;
    vector<int >str;
    vector<pair<vec*, pair<int, bool> > > source;

    // 随机数生成
    random_device rd;
    mt19937 gen(rd());

    // 枚举匹配字符串的长度
    vector<int> l_list;
    for (int l = l_min; l <= l_max; ++l) {
        l_list.push_back(l);
    }
    if (p < 1.0 - 1e-6) {
        shuffle(l_list.begin(), l_list.end(), gen);
    }
    for (int l : l_list) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                int sum = 0;
                // 枚举字符串起点
                for (int s = csd[i][j].len - l; s >= 0; --s) {
                    str.clear();
                    if (csd_matching[i][j].data[s + l - 1] == 0 || csd_matching[i][j].data[s] == 0) continue;
                    for (int k = 0; k < l; ++k)
                        str.push_back(csd_matching[i][j].data[s+k]);

                    sum = _string_matching(str, i, j, l, n, num);

                    if (sum > 0){
                        source.clear();
                        for (int k = 0; k < l; ++k) {
                            if (str[k] != 0) {
                                if (str[k] == 1) {
                                    source.push_back(make_pair(&pool_baseline[j], make_pair(k, false)));
                                } else if (str[k] == -1) {
                                    source.push_back(make_pair(&pool_baseline[j], make_pair(k, true)));
                                }
                            }
                        }
                        // 尝试增加加法器
                        pool_tmp[num++] = (vec(source));
                        if (_local_search(n, m, num, p, gen)) {
                            pool[num - 1] = (vec(source));
                            _string_matching_set0(str, i, j, l, n);
                        }
                        else {
                            _recover(str, i, j, l, n);
                            num--;
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
        for (int j = 0; j < m; j++){
            cin >> a[i][j];
        }
        // 由于位移没有代价，所以可以把公用的末尾0去掉
        while(1) {
            bool flag = true;
            bool no_zero = true;
            for (int j = 0; j < m; j++) {
                if (a[i][j] != 0) {
                    no_zero = false;
                }
                if (a[i][j] % 2 != 0) {
                    flag = false;
                    break;
                }
            }
            if(!flag || no_zero) break;
            for (int j = 0; j < m; j++) {
                a[i][j] /= 2;
            }
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            int l = log2(abs(a[i][j])) + 1;
            l = max(l, 1);
            csd[i][j].len = l;
            csd[i][j].INT2CSD(a[i][j]);
        }
    }

    cout << "translated goal:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++){
            cout << a[i][j] << " ";
        }
        cout << endl;
    }

    cout << "CSD representation:" <<endl;

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
    int base_cost_x = total_cost_x;
    int base_cost_b = total_cost_b;
    int base_cost_c = total_cost_c;
    double base_cost = total_cost;
    double no_random_cost;

    printf("Total cost: %d + %dc + %db.\n", total_cost_x, total_cost_c, total_cost_b);
    printf("Total cost when c = %f, b = %f: %f\n", c, b, total_cost);
    
    printf("\n");

    double best_cost = -1;
    int num = 0;

    for (double p = 1.0; p >= 0.0; p -= 1.0/ROUND) {
        current_cost_x = base_cost_x;
        current_cost_b = base_cost_b;
        current_cost_c = base_cost_c;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                csd_matching[i][j] = csd[i][j];
            }
        }

        memset(IsReplace, 0, sizeof(IsReplace));
        memset(pool, 0, sizeof(pool));
        memset(pool_tmp, 0, sizeof(pool_tmp));
                    
        int num = string_matching(n, m, 3, 10, p);

        for (int i = 0; i < n; i++) {
            vector<pair<vec*, pair<int, bool> > > source;
            for(int j = 0; j < m; j++) {
                for(int k = 0; k < csd[i][j].len; k++) {
                    if (IsReplace[i][j][k] == true) {
                        if (Replace[i][j][k] != -1){
                            source.push_back(make_pair(&pool[Replace[i][j][k]], make_pair(ShiftReplace[i][j][k], NegReplace[i][j][k])));
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
            total_cost += pool[i].cost();
            total_cost_x += pool[i].cost_x;
            total_cost_b += pool[i].cost_b;
            total_cost_c += pool[i].cost_c;
        }
        // printf("Total cost: %d + %dc + %db.\n", total_cost_x, total_cost_c, total_cost_b);
        // printf("Total cost when c = %f, b = %f: %f\n", c, b, total_cost);

        if (best_cost == -1 || total_cost < best_cost) {
            if (best_cost == -1) no_random_cost = total_cost;
            best_cost = total_cost;

            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++)
                    cout << csd_matching[i][j] << " ";
                cout << endl;
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
        }
    }
    cout<< "Best cost: " << best_cost << endl;
    cout<< "No random cost: " << no_random_cost << endl;
    cout<< "Baseline cost: " << base_cost << endl;
    return 0;
}

