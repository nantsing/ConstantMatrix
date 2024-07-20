#include <set>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <memory>
#include <cstring>
#include <cmath>
typedef unsigned int uint;

const double B = 0.0;
const double  C = 1.2;
const int w = 8;

template <uint width>
class CSD {
    public:
    int bits[width + 5];
    uint len;
    int num;
    CSD() {
        for (int i = 0; i < width + 5; i++) {
            bits[i] = 0;
        }
        len = 0;
        num = 0;
    }
    CSD(int x) {
        memset(bits, 0, sizeof(bits));
        num = x;
        len = 0;
        bool neg = false;
        if (x < 0) {
            x = -x;
            neg = true;
        }
        while (x) {
            bits[len++] = x & 1;
            x >>= 1;
        }
        for (int i = 0; i < (int)len - 1; i++) {
            if (bits[i] == 1 && bits[i + 1] == 1) {
                bits[i] = -1;
                bits[i + 1] = 0;
                bits[i + 2] += 1;
                if (i == len - 2)
                    len++;
            }
            uint k = i + 2;
            if (k >= len) continue;
            while (bits[k] == 2) {
                bits[k] = 0;
                bits[k + 1] += 1;
                k++;
            }
            len = std::max(len, k + 1);
        }
        if (neg) {
            for (int i = 0; i < len; i++) {
                bits[i] = -bits[i];
            }
        }
    }
    friend std::ostream& operator<<(std::ostream& os, const CSD& num) {
        for (int i = num.len - 1; i >= 0; i--)
            os << num.bits[i];
        return os;
    }
    CSD& operator=(const CSD& other) {
        if (this == &other) {
            return *this;
        }
        len = other.len;
        num = other.num;
        for (int i = 0; i < len; i++) {
            bits[i] = other.bits[i];
        }
        return *this;
    }
    bool operator==(const CSD& other) {
        if (len != other.len) return false;
        for (int i = 0; i < len; i++) {
            if (bits[i] != other.bits[i]) return false;
        }
        return true;
    }
    bool operator!=(const CSD& other) {
        return !(*this == other);
    }
    CSD operator -() {
        CSD<width> res;
        for (int i = 0; i < len; i++) {
            res.bits[i] = -bits[i];
        }
        res.len = len;
        res.num = -num;
        return res;
    }
    // get the sub word from [l, r)
    CSD get_sub(int l, int r) {
        CSD<width> res;
        for (int i = l; i < r; i++) {
            res.bits[i - l] = bits[i];
        }
        res.len = r - l;
        res.num = 0;
        for (int i = res.len - 1; i >= 0; i--) {
            res.num = res.num * 2 + res.bits[i];
        }
        return res;
    }
};

struct Cost {
    double x, c, b;
    Cost() : x(0), c(0), b(0) {}
    void operator+=(const Cost& other) {
        x += other.x;
        c += other.c;
        b += other.b;
    }
    void operator-=(const Cost& other) {
        x -= other.x;
        c -= other.c;
        b -= other.b;
    }
    Cost operator+(const Cost& other) {
        Cost res;
        res.x = x + other.x;
        res.c = c + other.c;
        res.b = b + other.b;
        return res;
    }
    Cost operator-() {
        Cost res;
        res.x = -x;
        res.c = -c;
        res.b = -b;
        return res;
    }
    double sum () const {
        return x + c * C + b * B;
    }
    bool operator<(const Cost& other) {
        return sum() < other.sum();
    }
};

template <uint width>
class Adder : public std::enable_shared_from_this< Adder<width> > {
    public:

    uint elements, bit_width;
    std::vector< CSD<width> > data; // total [$elements] of elements
    uint** srcid; // srcid[i][j] = the unique id of the src of the j-th bit of the i-th element
    bool** hasUpdate; // hasUpdate[i][j] = whether the j-th bit of the i-th element has been updated
    
    std::unordered_map<unsigned long long, std::pair< std::shared_ptr<Adder>, std::pair<int, bool> > > src;
    std::unordered_map<unsigned long long, std::pair< std::shared_ptr<Adder>, std::pair<int, bool> > > dst;
    bool isRoot;

    // alloc mem for srcid and hasUpdate
    void alloc_mem() {
        srcid = new uint*[elements];
        hasUpdate = new bool*[elements];
        for (int i = 0; i < elements; i++) {
            srcid[i] = new uint[width];
            hasUpdate[i] = new bool[width];
            for (int j = 0; j < width; j++) {
                srcid[i][j] = 0;
                hasUpdate[i][j] = false;
            }
        }
    }

    // empty adder
    Adder(uint elements): isRoot(false), elements(elements), bit_width(width) {
        alloc_mem();
        data.resize(elements);
    }

    // root elements
    Adder(int k, uint elements) : isRoot(true), elements(k), bit_width(width) {
        alloc_mem();
        data.resize(elements);
        data[k] = CSD<width>(1);
    }

    // update output width of the adder
    void update_width() {
        uint largest = (1 << width) - 1;
        uint sum = 0;
        for (int i = 0; i < elements; i++) {
            sum += std::abs(data[i].num) * largest;
        }
        bit_width = width;
        while ((1 << bit_width) - 1 < sum) {
            bit_width++;
        }
    }
    
    ~Adder() {
        for (int i = 0; i < elements; i++) {
            delete[] srcid[i];
            delete[] hasUpdate[i];
        }
        delete[] srcid;
        delete[] hasUpdate;
    }

    // the cost of the adder
    Cost cost() {
        Cost ret;
        int max_width = 0;
        int min_width = 100000;

        if (isRoot) {
            return ret;
        }
        for (auto it = src.begin(); it != src.end(); it++){
            auto t = it->second;
            int wi = t.first->bit_width + t.second.first;
            if (wi > max_width) {
                max_width = wi;
            }
            if (wi < min_width) {
                min_width = wi;
            }
            ret.x += wi;
            if (t.second.second) {
                ret.b += wi;
            }
        }
        ret.c = max_width;
        ret.x -= min_width;
        return ret;
    }
};

template <uint width>
class Circuit {
    public:
    unsigned long long count = 0; // counter for unique src/dst id

    uint n, m; // n-rows m-cols
    std::vector<std::shared_ptr<Adder<width> > > input;
    std::vector<std::shared_ptr<Adder<width> > > output;
    std::vector< std::vector< std::shared_ptr<Adder<width> > > > layers;

    Circuit(uint n, uint m, int* matrix) : n(n), m(m) {
        input.resize(m);
        output.resize(n);
        layers.resize(0);
        for (int i = 0; i < m; i++) {
            input[i] = std::make_shared<Adder<width> >(i, m);
        }
        for (int i = 0; i < n; i++, matrix += m) {
            output[i] = std::make_shared<Adder<width> >(m);
            naive_connect(output[i], matrix);
        }
    }

    Cost cost() {
        Cost ret;
        for (int i = 0; i < m; i++) {
            ret += input[i]->cost();
        }
        for (int i = 0; i < layers.size(); i++) {
            for (int j = 0; j < layers[i].size(); j++) {
                ret += layers[i][j]->cost();
            }
        }
        for (int i = 0; i < n; i++) {
            ret += output[i]->cost();
        }
        return ret;
    }
    
    // the starting connection of the adder
    void naive_connect(std::shared_ptr<Adder<width> >& cur, int coefficient[]) {
        for (int i = 0; i < m; i++) {
            cur->data[i] = CSD<width>(coefficient[i]);
            std::cout << cur->data[i] << " ";
            for (int j = 0; j < cur->data[i].len; j++) {
                 if (cur->data[i].bits[j] == 1) {
                    cur->srcid[i][j] = count++;
                    cur->src[cur->srcid[i][j]] = std::make_pair(input[i], std::make_pair(j, false));
                    input[i]->dst[count] = std::make_pair(cur, std::make_pair(j, false));
                } else if (cur->data[i].bits[j] == -1) {
                    cur->srcid[i][j] = count++;
                    cur->src[cur->srcid[i][j]] = std::make_pair(input[i], std::make_pair(j, true));
                    input[i]->dst[count] = std::make_pair(cur, std::make_pair(j, true));
                }
            }
        }
        std::cout << std::endl;
        cur->update_width();
    }

    void print_test_info() {
        for (int i = 0; i < m; i++) {
            std::cout << "input " << i << " : ";
            for (int j = 0; j < m; j++) {
                std::cout << input[i]->data[j].num << " ";
            }
            std::cout << " width:" << input[i]->bit_width << " cost:" << input[i]->cost().sum() << std::endl;
        }
        for (int i = 0; i < n; i++) {
            std::cout << "output " << i << " : ";
            for (int j = 0; j < m; j++) {
                std::cout << output[i]->data[j].num << " ";
            }
            std::cout << " width:" << output[i]->bit_width << " cost:" << output[i]->cost().sum() << std::endl;
            // for (auto it = output[i]->src.begin(); it != output[i]->src.end(); it++) {
            //     auto t = it->second;
            //     std::cout << "src " << it->first << " : " << " " << t.second.first << " " << t.second.second << std::endl;
            // }
        }
    }
};

int main() {
    int n, m;
    int *matrix;
    std::cin >> n >> m;
    matrix = new int[n * m];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            std::cin >> matrix[i * m + j];
        }
        // 由于位移没有代价，所以可以把公用的末尾0去掉
        while(1) {
            bool flag = true;
            bool no_zero = true;
            for (int j = 0; j < m; j++) {
                if (matrix[i * m + j] != 0) {
                    no_zero = false;
                }
                if (matrix[i * m + j] % 2 != 0) {
                    flag = false;
                    break;
                }
            }
            if(!flag || no_zero) break;
            for (int j = 0; j < m; j++) {
                matrix[i * m + j] /= 2;
            }
        }
    }

    std::cout << "translated goal:" << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            std::cout << matrix[i * m + j] << " ";
        }
        std::cout << std::endl;
    }

    Circuit<w> circuit(n, m, matrix);

    
    circuit.print_test_info();
    Cost cost = circuit.cost();
    std::cout << cost.sum() << std::endl;

    delete [] matrix;
    return 0;
}