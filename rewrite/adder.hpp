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
    uint** srcid; // srcid[i][j] = the unique id of the src of the j-th bit of the i-th element
    bool** hasUpdate; // hasUpdate[i][j] = whether the j-th bit of the i-th element has been updated
    
    std::map<uint , std::pair< std::shared_ptr<Adder>, std::pair<int, bool> > > src;
    std::map<Index, std::pair< std::shared_ptr<Adder>, std::pair<int, bool> > > dst;
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
    Adder(uint elements, uint nodeid, uint layerid): isRoot(false), elements(elements), bit_width(width), nodeid(nodeid), layerid(layerid) {
        alloc_mem();
        data.resize(elements);
    }

    // root elements
    Adder(int k, uint elements, uint nodeid, uint layerid) : isRoot(true), elements(elements), bit_width(width), nodeid(nodeid), layerid(layerid) {
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