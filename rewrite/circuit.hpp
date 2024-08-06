#pragma once
#include "adder.hpp"
#include <utility>
#include <vector>
#include <memory>
#include <unordered_map>
#include <random>
#include <cassert>
#include <algorithm>
#include <ranges>

struct difference{
    std::pair<int, int> index; // layer, node
    std::pair<int, int> shift; // shift1, shift2 for the target number
    std::pair<bool, bool> neg; // neg1, neg2 for the target number
    std::pair<int, bool> differ; // Δ shift, same/different sign
    friend bool operator < (const difference &a, const difference &b) {
        if (a.differ.first != b.differ.first) return a.differ.first < b.differ.first;
        else if (a.differ.second == true && b.differ.second == false) return false;
        return true;
    }
};

template <uint width>
class Circuit {
    // return : {nodeid, {shift, neg}}
    std::vector< std::pair<int, std::pair<int, bool> > > _sub_num_match(CSD<width> goal, int layer, int col);

    bool _check_col_update_nice(CSD<width> goal, int layer, int col, std::vector< std::pair<int, std::pair<int, bool> > >& match, int new_layer_id, int new_node_id); 

    bool _check_point_merge_update_nice(std::shared_ptr<Adder<width> > new_node, std::vector<std::shared_ptr<Adder<width>>>& layer, std::vector<difference>& differ, int new_layer_id, int new_node_id, bool keep_adderi, bool keep_adderj, std::pair<int, int> src_node, std::pair<int, int> new_shift, std::pair<bool, bool> new_neg);

    public:
    std::mt19937 gen;
    uint n, m; // n-rows m-cols
    double p; // probability factor
    std::vector<std::shared_ptr<Adder<width> > > input;
    std::vector< std::vector< std::shared_ptr<Adder<width> > > > layers;

    Circuit(uint n, uint m, int* matrix, double p = 1.0); // constructor

    Circuit(const Circuit& other); // copy constructor

    Cost cost(); // return the total cost of the circuit
    
    // the starting connection of the adder
    void naive_connect(std::shared_ptr<Adder<width> >& cur, int coefficient[]);

    void col_wise_optimization(int l_min = 3, int l_max = 5, int layer = 0); // optimize the circuit in row-wise

    void print_test_info(); // print the test info

    void print_csd_form(); // print the csd form of the circuit

    void print_dst(); // print the dst of the points

    void point_merge (int layer, bool is_input = false); // merge the points in the same layer
};

template <uint width>
void Circuit<width>::print_dst() {
    for (int i = 0; i < m; i++) {
        std::cout << "input " << i << " : ";
        for (auto& dst : input[i]->dst) {
            std::cout << "(" << dst.first.layer << " " << dst.first.node << " " << dst.first.bitshift << ") ";
        }
        std::cout << std::endl;
    }
    for (int i = layers.size() - 1; i > 0; i--) {
        for (int j = 0; j < layers[i].size(); j++) {
            std::cout << "layer " << i << " node " << j << " : ";
            for (auto& dst : layers[i][j]->dst) {
                std::cout << "(" << dst.first.layer << " " << dst.first.node << " " << dst.first.bitshift << ") ";
            }
            std::cout << std::endl;
        }
    }
}

template <uint width>
void Circuit<width>::print_csd_form() {
    for (int i = 0; i < m; i++) {
        std::cout << "input " << i << " : ";
        for (int j = 0; j < m; j++) {
            std::cout << input[i]->data[j] << "\t";
        }
        std::cout << std::endl;
    }
    for (int i = layers.size() - 1; i >= 0; i--) {
        for (int j = 0; j < layers[i].size(); j++) {
            std::cout << "layer " << i << " node " << j << " : ";
            for (int k = 0; k < m; k++) {
                std::cout << layers[i][j]->data[k] << "\t";
            }
            std::cout << std::endl;
        }
    }
}

template <uint width>
Circuit<width>::Circuit(uint n, uint m, int* matrix, double p) : n(n), m(m), p(p) {
    std::random_device rd;
    gen = std::mt19937(rd());
    input.resize(m);
    layers.resize(1);
    layers[0].resize(n);
    for (int i = 0; i < m; i++) {
        input[i] = std::make_shared<Adder<width> >(i, m, i, 0xffffffff);
    }
    for (int i = 0; i < n; i++, matrix += m) {
        layers[0][i] = std::make_shared<Adder<width> >(m, i, 0);
        naive_connect(layers[0][i], matrix);
    }
}

template <uint width>
Circuit<width>::Circuit(const Circuit& other) : n(other.n), m(other.m), p(other.p) {
    std::random_device rd;
    gen = std::mt19937(rd());
    
    // copy the size
    input.resize(other.input.size());
    layers.resize(other.layers.size());
    for (int i = 0; i < layers.size(); i++) {
        layers[i].resize(other.layers[i].size());
    }

    // copy the data
    for (int i = 0; i < input.size(); i++) {
        input[i] = std::make_shared<Adder<width> >(*other.input[i]);
    }
    for (int i = 0; i < layers.size(); i++) {
        for (int j = 0; j < layers[i].size(); j++) {
            layers[i][j] = std::make_shared<Adder<width> >(*other.layers[i][j]);
        }
    }

    // revise the src and dst
    for (int i = 0; i < input.size(); i++) {
        for (auto& dst : input[i]->dst) {
            auto index = dst.first;
            auto src = other.input[i]->dst[index].first;
            dst.second.first = layers[src.lock()->layerid][src.lock()->nodeid];
        }
    }
    for (int i = 0; i < layers.size(); i++) {
        for (int j = 0; j < layers[i].size(); j++) {
            for (auto& src : layers[i][j]->src) {
                auto index = src.first;
                auto dst = other.layers[i][j]->src[index];
                if (dst.first->layerid == 0xffffffff) {
                    src.second.first = input[dst.first->nodeid];
                } else {
                    src.second.first = layers[dst.first->layerid][dst.first->nodeid];
                }
            }
            for (auto& dst : layers[i][j]->dst) {
                auto index = dst.first;
                auto src = other.layers[i][j]->dst[index].first;
                dst.second.first = layers[src.lock()->layerid][src.lock()->nodeid];
            }
        }
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
            if (i == 0 || layers[i][j]->dst.size() > 0)
                ret += layers[i][j]->cost();
        }
    }
    return ret;
}

template <uint width>
void Circuit<width>::naive_connect(std::shared_ptr<Adder<width> >& cur, int coefficient[]) {
    for (int i = 0; i < m; i++) {
        cur->data[i] = CSD<width>(coefficient[i]);
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
    for (int i = layers.size() - 1; i >= 0; i--) {
        for (int j = 0; j < layers[i].size(); j++) {
            std::cout << "layer " << i << " node " << j << " : ";
            for (int k = 0; k < m; k++) {
                std::cout << layers[i][j]->data[k].num << " ";
            }
            std::cout << " width:" << layers[i][j]->bit_width << " cost:" << layers[i][j]->cost().sum() << std::endl;
        }
    }
}

template <uint width>
void Circuit<width>::col_wise_optimization(int l_min, int l_max, int layer) {
    // 现在只能处理最靠近输入的一层
    
    int new_layer_id = layers.size();
    layers.resize(layers.size() + 1);
    
    // 枚举匹配字符串的长度
    std::vector<int> l_list;
    for (int l = l_min; l <= l_max; ++l) {
        l_list.push_back(l);
    }
    if (p < 1.0 - 1e-6) {
        std::shuffle(l_list.begin(), l_list.end(), gen);
    }

    for (int l : l_list) {
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < layers[layer].size(); i++) {
                // 枚举字符串起点
                for (int s = layers[layer][i]->data[j].len - l; s >= 0; s--) {
                    int e = s + l;
                    assert(e <= layers[layer][i]->data[j].len);
                    if (layers[layer][i]->data[j].bits[s] == 0 || layers[layer][i]->data[j].bits[e - 1] == 0) continue;
                    if (!layers[layer][i]->check_unupdated(j, s, e)) continue;
                    auto sub_num = layers[layer][i]->data[j].get_sub(s, e);
                    auto match = _sub_num_match(sub_num, layer, j);
                    assert(match.size() > 0);
                    if (match.size() == 1) continue;
                    // 如果需要取负的较多，就把goal改为-goal
                    int pos = 0, neg = 0;
                    for (auto& m : match) {
                        if (m.second.second) neg++;
                        else pos++;
                    }
                    if (neg > pos) {
                        sub_num = -sub_num;
                        for (auto& m : match) {
                            m.second.second = !m.second.second;
                        }
                    }
                    _check_col_update_nice(sub_num, layer, j, match, new_layer_id, layers[new_layer_id].size());
                }
            }
        }
    }
    if (layers[new_layer_id].size() == 0) {
        layers.resize(layers.size() - 1);
    }
}

template <uint width>
std::vector< std::pair<int, std::pair<int, bool> > > Circuit<width>::_sub_num_match(CSD<width> goal, int layer, int col) {
    std::vector< std::pair<int, std::pair<int, bool> > > ret;
    for (int i = 0; i < layers[layer].size(); i++) {
        for (int s = layers[layer][i]->data[col].len - goal.len; s >= 0; s--) {
            int e = s + goal.len;
            assert(e <= layers[layer][i]->data[col].len);
            if (!layers[layer][i]->check_unupdated(col, s, e)) continue;
            auto sub_num = layers[layer][i]->data[col].get_sub(s, e);
            if (sub_num == goal) {
                ret.push_back(std::make_pair(i, std::make_pair(s, false)));
                s -= goal.len - 1;
            } else if (-sub_num == goal) {
                ret.push_back(std::make_pair(i, std::make_pair(s, true)));
                s -= goal.len - 1;
            }
        }
    }
    return ret;
}

template <uint width>
bool Circuit<width>::_check_col_update_nice(CSD<width> goal, int layer, int col, std::vector< std::pair<int, std::pair<int, bool> > >& match, int new_layer_id, int new_node_id) {
    std::shared_ptr< Adder<width> > new_node = std::make_shared<Adder<width> >(m, new_node_id, new_layer_id);
    int co[m];
    for (int i = 0; i < m; i++) {
        co[i] = 0;
    }
    co[col] = goal.num;
    naive_connect(new_node, co);
    
    Cost old_cost;
    int last = -1;
    Cost new_cost = new_node->cost();
    std::shared_ptr< Adder<width> > cur;
    for (auto& m : match) {
        if (m.first != last) {
            if (last != -1) {
                new_cost += cur->cost();
                old_cost += layers[layer][last]->cost();
            }
            cur = std::make_shared< Adder<width> >(*layers[layer][m.first]);
            cur->set_temp();
        }
        last = m.first;
        int id = cur->counter++;
        cur->src[id] = std::make_pair(new_node, m.second);
        for (int i = m.second.first; i < m.second.first + goal.len; i++) {
            cur->hasUpdate[col][i] = true;
            if (cur->srcid[col][i] != -1) {
                cur->src.erase(cur->srcid[col][i]);
            }
            cur->srcid[col][i] = id;
        }
    }
    new_cost += cur->cost();
    old_cost += layers[layer][last]->cost();

    std::uniform_real_distribution<double> distr(0, 1);

    if ((new_cost.sum() < old_cost.sum() && distr(gen) <= p) || 
        (new_cost.sum() >= old_cost.sum() && distr(gen) < 1 - p)) {
        for (auto& m : match) {
            auto cur = layers[layer][m.first];
            int new_id = cur->counter++;
            for (int i = m.second.first; i < m.second.first + goal.len; i++) {
                cur->hasUpdate[col][i] = true;
                if (cur->srcid[col][i] != -1) {
                    Index idx(cur->layerid, cur->nodeid, i);
                    if (cur->src.find(cur->srcid[col][i]) != cur->src.end()) {
                        cur->src[cur->srcid[col][i]].first->dst.erase(idx);
                        cur->src.erase(cur->srcid[col][i]);
                    }
                }
                cur->srcid[col][i] = new_id;
            }
            cur->src[new_id] = std::make_pair(new_node, m.second);
            new_node->dst[Index(cur->layerid, cur->nodeid, m.second.first)] = std::make_pair(cur, m.second);
        }
        layers[new_layer_id].push_back(new_node);
        return true;
    }

    return false;
}

template <uint width>
void Circuit<width>::point_merge(int layer, bool is_input) {
    int new_layer_id = layers.size();
    layers.resize(layers.size() + 1);

    std::vector<std::shared_ptr<Adder<width>>>* layer_p;

    if (is_input) layer_p = &input;
    else layer_p = &layers[layer];

    // 枚举当前层中的两个点看能否合并
    for (int i = 0; i < (*layer_p).size(); i++) {
        for (int j = 0; j < i; j++) {
            std::vector<difference> diff; // 存储同一目标点中发生的位移和取负的差别
            for (auto iter_i = (*layer_p)[i]->dst.begin(); iter_i != (*layer_p)[i]->dst.end(); iter_i++) {
                assert(iter_i->first.bitshift == iter_i->second.second.first);
                auto node_info = std::make_pair(iter_i->first.layer, iter_i->first.node);

                Index start(node_info.first, node_info.second, INT32_MIN);
                Index end(node_info.first, node_info.second, INT32_MAX);
                auto it_start = (*layer_p)[j]->dst.lower_bound(start);
                auto it_end = (*layer_p)[j]->dst.upper_bound(end);
                for (auto iter_j = it_start; iter_j != it_end; iter_j++) {
                    assert(node_info.first == iter_j->first.layer);
                    assert(node_info.second == iter_j->first.node);
                    assert(iter_j->first.bitshift == iter_j->second.second.first);

                    auto delta = std::make_pair(iter_i->first.bitshift - iter_j->first.bitshift, iter_i->second.second.second == iter_j->second.second.second);
                    auto shift = std::make_pair(iter_i->first.bitshift, iter_j->first.bitshift);
                    auto neg = std::make_pair(iter_i->second.second.second, iter_j->second.second.second);
                    diff.push_back(difference{node_info, shift, neg, delta});
                }
            }
            
            if (diff.size() == 0) continue;

            // find the differ with most number of use
            std::sort(diff.begin(), diff.end());
            
            std::vector<difference> best_diff;
            std::pair<int, bool> last;
            int len = 0;
            for (int k = 0; k <= diff.size(); k++) {
                if (k != 0 && (k == diff.size() || last != diff[k].differ)) {
                    bool better = false;
                    if (len > best_diff.size()) {
                        better = true;
                    } else if (len == best_diff.size()) {
                        if (abs(last.first) < abs(best_diff[0].differ.first)) {
                            better = true;
                        } else if (abs(last.first) == abs(best_diff[0].differ.first) && last.second == false) {
                            better = true;
                        }
                    }

                    if (better) {
                        best_diff.clear();
                        for (int l = k - len; l < k; l++) {
                            best_diff.push_back(diff[l]);
                        }
                    }
                    len = 0;
                    if (k == diff.size()) break;
                }
                last = diff[k].differ;
                len++;
            }
            std::cout << i << " " << j << std::endl;
            for (auto& it: best_diff) {
                std::cout << it.index.second << std::endl;
            }

            // make the new adder
            bool keep_adderi = true, keep_adderj = true;
            if (!is_input) {
                if ((*layer_p)[i]->dst.size() == best_diff.size()) {
                    keep_adderi = false;
                }
                if ((*layer_p)[j]->dst.size() == best_diff.size()) {
                    keep_adderj = false;
                }
            }
            int new_node_id = layers[new_layer_id].size();
            auto new_node = std::make_shared<Adder<width>>(m, new_node_id, new_layer_id);
            
            // get the shift number and neg for the new_node
            std::pair<int, int> new_shift;
            std::pair<bool, bool> new_neg;
            if (best_diff[0].differ.first > 0) {
                new_shift = std::make_pair(best_diff[0].differ.first, 0);
            } else {
                new_shift = std::make_pair(0, -best_diff[0].differ.first);
            }

            // check pos or neg is better
            uint sum_pos = 0, sum_neg = 0;
            for (auto& it: best_diff) {
                if (it.neg.first) {
                    sum_neg++;
                } else {
                    sum_pos++;
                }
            }
            if (best_diff[0].differ.second) {
                new_neg = sum_neg > sum_pos ? std::make_pair(true, true) : std::make_pair(false, false);
            } else {
                new_neg = sum_pos > sum_neg ? std::make_pair(false, true) : std::make_pair(true, false);
            }

            // for (auto& it: best_diff) {
            //     it.shift.first -= new_shift.first;
            //     it.shift.second -= new_shift.second;
            //     it.neg.first = (it.neg.first != new_neg.first);
            //     it.neg.second = (it.neg.second != new_neg.second);
            // }

            if (keep_adderi && keep_adderj) {
                new_node->connect((*layer_p)[i], new_shift.first, new_neg.first);
                new_node->connect((*layer_p)[j], new_shift.second, new_neg.second);
            } else if (!keep_adderi &&!keep_adderj) {
                int co[m];
                for (int k = 0; k < m; k++) {
                    co[k] = (*layer_p)[i]->data[k].num + (*layer_p)[j]->data[k].num;
                }
                naive_connect(new_node, co);
            } else if (keep_adderi &&!keep_adderj) {
                int co[m];
                for (int k = 0; k < m; k++) {
                    co[k] = (*layer_p)[j]->data[k].num;
                }
                naive_connect(new_node, co);
                new_node->connect((*layer_p)[i], new_shift.first, new_neg.first);
            } else if(!keep_adderi && keep_adderj) {
                int co[m];
                for (int k = 0; k < m; k++) {
                    co[k] = (*layer_p)[i]->data[k].num;
                }
                naive_connect(new_node, co);
                new_node->connect((*layer_p)[j], new_shift.second, new_neg.second);                
            }
            new_node->update_width();
            // std::cout << "New Node Information:" << std::endl;
            // std::cout << "ID: " << new_node_id << std::endl;
            // std::cout << "Layer: " << new_layer_id << std::endl;
            // std::cout << "Shift: " << new_shift.first << ", " << new_shift.second << std::endl;
            // std::cout << "Neg: " << new_neg.first << ", " << new_neg.second << std::endl;
            _check_point_merge_update_nice(new_node, *layer_p, best_diff, new_layer_id, new_node_id, keep_adderi, keep_adderj, std::make_pair(i, j), new_shift, new_neg);
        }
    }
    if (layers[new_layer_id].size() == 0) {
        layers.resize(layers.size() - 1);
    }
}

template <uint width>
bool Circuit<width>::_check_point_merge_update_nice(std::shared_ptr< Adder<width> > new_node, 
std::vector<std::shared_ptr<Adder<width>>>& layer, std::vector<difference>& differ, 
int new_layer_id, int new_node_id, bool keep_adderi, bool keep_adderj, std::pair<int, int> src_node,
std::pair<int, int> new_shift, std::pair<bool, bool> new_neg) {

    Cost old_cost;
    if (!keep_adderi) {
        old_cost += layer[src_node.first]->cost();
    }
    if (!keep_adderj) {
        old_cost += layer[src_node.second]->cost();
    }
    // 需要考虑同一个dst中用到两次src的情况
    std::ranges::sort(differ, {}, &difference::index);
    auto last = std::make_pair(-1, -1);
    
    Cost new_cost = new_node->cost();
    std::shared_ptr< Adder<width> > cur;
    for (auto& it : differ) {
        if (it.index != last) {
            if (last != std::make_pair(-1, -1)) {
                new_cost += cur->cost();
                old_cost += layers[last.first][last.second]->cost();
            }
            cur = std::make_shared< Adder<width> >(*layers[it.index.first][it.index.second]);
            cur->set_temp();
        }
        last = it.index;
        int id = cur->counter++;

        int shift = it.shift.first - new_shift.first;
        bool neg = (it.neg.first != new_neg.first);

        cur->src[id] = std::make_pair(new_node, std::make_pair(shift, neg));

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < layer[src_node.first]->data[i].len; j++) {
                assert (j + it.shift.first < cur->data[i].len);
                if (cur->srcid[i][j + it.shift.first] != -1) {
                    cur->src.erase(cur->srcid[i][j + it.shift.first]);
                }
                cur->srcid[i][j + it.shift.first] = id;
            }
        }
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < layer[src_node.second]->data[i].len; j++) {
                assert (j + it.shift.second < cur->data[i].len);
                if (cur->srcid[i][j + it.shift.second] != -1) {
                    cur->src.erase(cur->srcid[i][j + it.shift.second]);
                }
                cur->srcid[i][j + it.shift.second] = id;
            }
        }
    }
    new_cost += cur->cost();
    old_cost += layers[last.first][last.second]->cost();

    std::uniform_real_distribution<double> distr(0, 1);

    if ((new_cost.sum() < old_cost.sum() && distr(gen) <= p) || 
        (new_cost.sum() >= old_cost.sum() && distr(gen) < 1 - p)) {
        for (auto& it : differ) {

            auto cur = layers[it.index.first][it.index.second];
            int id = cur->counter++;

            int shift = it.shift.first - new_shift.first;
            bool neg = (it.neg.first != new_neg.first);

            cur->src[id] = std::make_pair(new_node, std::make_pair(shift, neg));
            new_node->dst[Index(cur->layerid, cur->nodeid, shift)] = std::make_pair(cur, std::make_pair(shift, neg));

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < layer[src_node.first]->data[i].len; j++) {
                    assert (j + it.shift.first < cur->data[i].len);
                    if (cur->srcid[i][j + it.shift.first] != -1) {
                        Index idx(cur->layerid, cur->nodeid, j + it.shift.first);
                        if (cur->src.find(cur->srcid[i][j + it.shift.first]) != cur->src.end()) {
                            cur->src[cur->srcid[i][j + it.shift.first]].first->dst.erase(idx);
                            cur->src.erase(cur->srcid[i][j + it.shift.first]);
                        }
                    }
                    cur->srcid[i][j + it.shift.first] = id;
                }
            }
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < layer[src_node.second]->data[i].len; j++) {
                    assert (j + it.shift.second < cur->data[i].len);
                    if (cur->srcid[i][j + it.shift.second] != -1) {
                        Index idx(cur->layerid, cur->nodeid, j + it.shift.second);
                        if (cur->src.find(cur->srcid[i][j + it.shift.second]) != cur->src.end()) {
                            cur->src[cur->srcid[i][j + it.shift.second]].first->dst.erase(idx);
                            cur->src.erase(cur->srcid[i][j + it.shift.second]);
                        }
                    }
                    cur->srcid[i][j + it.shift.second] = id;
                }
            }
        }
        layers[new_layer_id].push_back(new_node);

        return true;
    }

    return false;
}
