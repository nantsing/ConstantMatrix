#pragma once

extern const double B, C;
#include "csd.hpp"
#include <vector>
#include <memory>
#include <map>

class Cost {
    public:
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

class Index {
    public:
    int layer, node, bitshift;
    Index(int layer, int node, int bit) : layer(layer), node(node), bitshift(bit) {}
    Index() : layer(0), node(0), bitshift(0) {}
    bool operator < (const Index& other) const {
        if (layer != other.layer) return layer < other.layer;
        if (node != other.node) return node < other.node;
        return bitshift < other.bitshift;
    }
};

template <uint width>
class Adder : public std::enable_shared_from_this< Adder<width> > {
    public:
    uint counter = 0;
    uint nodeid, layerid; // the id of the node in the layer
    uint elements, bit_width;
    std::vector< CSD<width> > data; // total [$elements] of elements
    int** srcid; // srcid[i][j] = the unique id of the src of the j-th bit of the i-th element
    bool** hasUpdate; // hasUpdate[i][j] = whether the j-th bit of the i-th element has been updated
    
    std::map<uint , std::pair< std::shared_ptr<Adder>, std::pair<int, bool> > > src;
    std::map<Index, std::pair< std::shared_ptr<Adder>, std::pair<int, bool> > > dst;
    bool isRoot;

    bool check_unupdated(int j, int s, int e); // check the s-th to e-th bits of the j-th element are unupdated

    void print_srcid(); // print the srcid matrix

    // alloc mem for srcid and hasUpdate
    void alloc_mem();

    // empty adder
    Adder(uint elements, uint nodeid, uint layerid);

    // root elements
    Adder(int k, uint elements, uint nodeid, uint layerid);

    // copy constructor
    Adder(const Adder& other);

    // update output width of the adder
    void update_width();
    
    ~Adder();

    // the cost of the adder
    Cost cost();
};

template <uint width>
bool Adder<width>::check_unupdated(int j, int s, int e) {
    for (int k = s; k < e; k++) {
        if (hasUpdate[j][k]) {
            return false;
        }
    }
    return true;
}

template <uint width>
void Adder<width>::print_srcid() {
    for (int i = 0; i < elements; i++) {
        for (int j = 0; j < width + EXTRA_BITS; j++) {
            std::cout << srcid[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

template <uint width>
void Adder<width>::alloc_mem() {
    srcid = new int*[elements];
    hasUpdate = new bool*[elements];
    for (int i = 0; i < elements; i++) {
        srcid[i] = new int[width + EXTRA_BITS];
        hasUpdate[i] = new bool[width + EXTRA_BITS];
        for (int j = 0; j < width + EXTRA_BITS; j++) {
            srcid[i][j] = -1;
            hasUpdate[i][j] = false;
        }
    }
}

template <uint width>
Adder<width>::Adder(uint elements, uint nodeid, uint layerid): isRoot(false), elements(elements), bit_width(width), nodeid(nodeid), layerid(layerid) {
    alloc_mem();
    data.resize(elements);
}

template <uint width>
Adder<width>::Adder(int k, uint elements, uint nodeid, uint layerid) : isRoot(true), elements(elements), bit_width(width), nodeid(nodeid), layerid(layerid) {
    alloc_mem();
    data.resize(elements);
    data[k] = CSD<width>(1);
}

template <uint width>
Adder<width>::Adder(const Adder& other) {
    counter = other.counter;
    elements = other.elements;
    bit_width = other.bit_width;
    nodeid = other.nodeid;
    layerid = other.layerid;
    isRoot = other.isRoot;
    alloc_mem();
    data = other.data;
    for (int i = 0; i < elements; i++) {
        for (int j = 0; j < width + EXTRA_BITS; j++) {
            srcid[i][j] = other.srcid[i][j];
            hasUpdate[i][j] = other.hasUpdate[i][j];
        }
    }
    src = other.src;
    dst = other.dst;
}

template <uint width>
void Adder<width>::update_width() {
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

template <uint width>
Adder<width>::~Adder() {
    for (int i = 0; i < elements; i++) {
        delete[] srcid[i];
        delete[] hasUpdate[i];
    }
    delete[] srcid;
    delete[] hasUpdate;
}

template <uint width>
Cost Adder<width>::cost() {
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