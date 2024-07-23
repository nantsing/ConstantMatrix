#pragma once
#include "adder.hpp"
#include <vector>
#include <memory>
#include <unordered_map>

template <uint width>
class Circuit {
    public:
    uint n, m; // n-rows m-cols
    std::vector<std::shared_ptr<Adder<width> > > input;
    std::vector<std::shared_ptr<Adder<width> > > output;
    std::vector< std::vector< std::shared_ptr<Adder<width> > > > layers;

    Circuit(uint n, uint m, int* matrix); // constructor

    Cost cost(); // return the total cost of the circuit
    
    // the starting connection of the adder
    void naive_connect(std::shared_ptr<Adder<width> >& cur, int coefficient[]);

    void print_test_info(); // print the test info
};

template <uint width>
Circuit<width>::Circuit(uint n, uint m, int* matrix) : n(n), m(m) {
    input.resize(m);
    output.resize(n);
    layers.resize(0);
    for (int i = 0; i < m; i++) {
        input[i] = std::make_shared<Adder<width> >(i, m, i, 0xffffffff);
    }
    for (int i = 0; i < n; i++, matrix += m) {
        output[i] = std::make_shared<Adder<width> >(m, i, 0);
        naive_connect(output[i], matrix);
    }
}

template <uint width>
Cost Circuit<width>::cost() {
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

template <uint width>
void Circuit<width>::naive_connect(std::shared_ptr<Adder<width> >& cur, int coefficient[]) {
    for (int i = 0; i < m; i++) {
        cur->data[i] = CSD<width>(coefficient[i]);
        std::cout << cur->data[i] << " ";
        for (int j = 0; j < cur->data[i].len; j++) {
            Index dst_info(cur->layerid, cur->nodeid, j);
            // std::cout << dst_info.node << " " << dst_info.bitshift << " " << dst_info.layer << " " << std::endl;
             if (cur->data[i].bits[j] == 1) {
                cur->srcid[i][j] = cur->counter++;
                cur->src[cur->srcid[i][j]] = std::make_pair(input[i], std::make_pair(j, false));
                input[i]->dst[dst_info] = std::make_pair(cur, std::make_pair(j, false));
            } else if (cur->data[i].bits[j] == -1) {
                cur->srcid[i][j] = cur->counter++;
                cur->src[cur->srcid[i][j]] = std::make_pair(input[i], std::make_pair(j, true));
                input[i]->dst[dst_info] = std::make_pair(cur, std::make_pair(j, true));
            }
        }
    }
    std::cout << std::endl;
    cur->update_width();
}

template <uint width>
void Circuit<width>::print_test_info() {
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
    }
}