#include <iostream>
#include <cstring>
#include <cmath>
#include "circuit.hpp"

typedef unsigned int uint;

const double B = 0.0;
const double  C = 1.2;
const int w = 8;


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

    circuit.col_wise_optimization();
    circuit.print_test_info();
    cost = circuit.cost();
    std::cout << cost.sum() << std::endl;
    
    delete [] matrix;
    return 0;
}