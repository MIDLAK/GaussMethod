//Решение СЛАУ методом Гаусса с выбором главного элемента
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <math.h>

void print_matrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("\t%.2lf", matrix[i][j]);
        }
        printf("\t| %.2lf", matrix[i][cols]);
        std::cout << "\n";
    }
    std::cout << "\n";
}

void matrix_filling(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols+1; j++) {
            matrix[i][j] = 100 + rand() % (900 - 100 + 1);
        }
    }
}

int main(int arg, char *argv[]) {
    int cols;
    std::cout << "N> ";
    std::cin >> cols;

    int rows = cols;

    double** matrix = new double* [rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols+1]; //+1 т.к. вектро b хранится в отдельном столбце
    }

    matrix_filling(matrix, rows, cols);

    print_matrix(matrix, rows, cols);

    double max_el;
    int max_el_row;
    int permutations; //количество перестановок строк

    for (int j = 0; j < cols; j++) {
        max_el = matrix[j][j];
        max_el_row = j;
        for(int i = j; i < rows; i++) {
            if (fabs(matrix[i][j]) > max_el) {
                max_el = matrix[i][j];
                max_el_row = i;
                permutations++;
            }
        }
        //обмен строк местами
        double* swp_address = matrix[max_el_row];
        matrix[max_el_row] = matrix[j]; //j совпадает с номером нужной строки
        matrix[j] = swp_address;

        //прямой ход (нули ниже главной диагонали)
        for (int k = 1 + j; k < rows; k++) {
            double c = matrix[k][j]/matrix[j][j];
            for (int m = 0; m <= rows; m++) {
                if (fabs(c) == 0) {
                    break;
                }
                double c2 = (matrix[j][m]*c);
                double t = matrix[k][m] - c2;
                matrix[k][m] = t;
            }
        }
        print_matrix(matrix, rows, cols);
    }

    double x[cols+1];
    //обратный ход
    for (int i = rows-1; i >= 0; i--) {
        double sum = 0;
        for (int j = cols-1; j > i; j--) {
            sum += x[j]*matrix[i][j];
        }
        x[i] = (matrix[i][cols] - sum)/matrix[i][i];
        if (x[i] == 0) {
            x[i] = fabs(x[i]);
        }
    }

    double determinant = 1.0;
    for (int i = 0; i < rows; i++) {
        determinant *= matrix[i][i];
    }
    determinant *= pow(-1, permutations);

    //печать результатов
    printf("x = [");
    for (int i = 0; i < cols; i++) {
        printf("\t %.2lf", x[i]);
    }
    printf("\t]");
    printf("\ndet(A) = %.2lf", determinant);
}