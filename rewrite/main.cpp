#include <iostream>
#include <cstring>
#include <cmath>
// #define DEBUG
#include "circuit.hpp"

typedef unsigned int uint;

const double B = 0.0;
const double  C = 1.2;
const int w = 16;
#define ROUND 10000


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

    Circuit<w> best_circuit(n, m, matrix);
    Circuit<w> best_circuit_copy(best_circuit);

    Cost best_cost = best_circuit.cost();
    std::cout << best_cost.sum() << std::endl;


    for (double p = 1.0; p >= 0.0; p -= 1.0 / ROUND) {
        Circuit<w> new_circuit(best_circuit_copy);
        new_circuit.p = p;
        new_circuit.col_wise_optimization(3, 10);
        new_circuit.point_merge(1);
        new_circuit.point_merge(-1, true);
        Cost cost = new_circuit.cost();
        if (cost < best_cost) {
            best_cost = cost;
            best_circuit = new_circuit;
        }
    }
    best_circuit.check();
    best_circuit.print_csd_form();
    best_circuit.print_adder_value();
    best_circuit.print_dst();
    best_circuit.print_src();
    best_cost = best_circuit.cost();
    std::cout << best_cost.sum() << std::endl;

    delete [] matrix;
    return 0;
}