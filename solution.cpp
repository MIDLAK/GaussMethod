//Решение СЛАУ методом Гаусса с выбором главного элемента
#include <cstdlib>
#include <fstream>
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

int matrix_filling(double** matrix, int rows, int cols) {
    std::ifstream f;
    f.open("file.txt");
    if (f.fail()) {
        return -1;
    }
    double tmp;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j <= cols; j++) {
            f >> matrix[i][j];
        }
    }
    f.close();
    return 0;
}

//количество строк в файле (размерность матрицы)
int file_length(char* path) {
    std::ifstream file;
    int rows = 0;
    std::string line;

    file.open(path);
    if (file.fail()) {
        return -1;
    }

    while (getline(file, line)) {
        rows++;
    }
    file.close();

    return rows;
}

double** inverse_matrix(double** matrix, int rows, int cols) {

    double** matrix_inverse = new double* [rows];

    //генерация единичной матрицы
    for (int i = 0; i < rows; i++) {
        matrix_inverse[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            if (i == j) {
                matrix_inverse[i][j] = 1;
            } else {
                matrix_inverse[i][j] = 0;
            }
        }
    }
    
    double cft;
    for(int t = 0; t < rows; t++) {
        cft = matrix[t][t];
        //деление текущей строки на диагональный элемент
        for (int j = 0; j < cols; j++) {
            matrix[t][j] /= cft;
            matrix_inverse[t][j] /= cft;
        }

        //вычитание из следующей строки предудующей с домножением
        for (int i = t + 1; i < rows; i++) {
            cft = matrix[i][t]; //первый элемент текущей строки
            for (int j = 0; j < cols; j++) {
                matrix[i][j] -= matrix[t][j] * cft;
                matrix_inverse[i][j] -= matrix_inverse[t][j] * cft;
            }
        }
    }

    for (int t = rows - 1; t > 0; t--) {
        for (int i = t - 1; i >= 0; i--) {
            cft = matrix[i][t]; //последний элемент текущей строки
            for (int j = 0; j < cols; j++) {
                matrix[i][j] -= matrix[t][j] * cft;
                matrix_inverse[i][j] -= matrix_inverse[t][j] * cft;
            }
        }
    }

    return matrix_inverse;
}

int main(int arg, char *argv[]) {

    int rows = file_length(argv[1]);
    int cols = rows;

    double** matrix = new double* [rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols+1]; //+1 т.к. вектро b хранится в отдельном столбце
    }

    matrix_filling(matrix, rows, cols);

    double** matrix_copy = new double* [rows];
    double** matrix_copy_2 = new double* [rows];
    for (int i = 0; i < rows; i++) {
        matrix_copy[i] = new double[cols+1];
        matrix_copy_2[i] = new double[cols+1];
        for (int j = 0; j < cols+1; j++) {
            matrix_copy[i][j] = matrix[i][j];
            matrix_copy_2[i][j] = matrix[i][j];
        }
    }

    double max_el;
    int max_el_row;
    int permutations; //количество перестановок строк

    for (int j = 0; j < cols; j++) {

        //поиск максимального по модулю элемента столбца
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
            for (int m = 0; m <= cols; m++) {
                if (fabs(c) == 0) {
                    break;
                }
                double c2 = (matrix[j][m]*c);
                double t = matrix[k][m] - c2;
                matrix[k][m] = t;
            }
        }
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

    //определитель
    double determinant = 1.0;
    for (int i = 0; i < rows; i++) {
        determinant *= matrix[i][i];
    }
    determinant *= pow(-1, permutations);
    if (determinant == 0) {
        std::cout << "No solutions" << std::endl;
        return 0;
    }

    double** matrix_inverse = inverse_matrix(matrix_copy_2, rows, cols);

    //определение невязок
    double residuals[rows];
    for (int i = 0; i < rows; i++) {
        double sum = 0;
        for (int j = 0; j < cols; j++) {
           sum += matrix_copy[i][j] * x[j];
        }
        residuals[i] = matrix_copy[i][cols] - sum;
    }

    //печать результатов
    print_matrix(matrix_copy, rows, cols);
    print_matrix(matrix, rows, cols);

    //вывод обратной матрицы
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("\t%.2lf", matrix_inverse[i][j]);
        }
        printf("\n");
    }

    printf("x = [");
    for (int i = 0; i < cols; i++) {
        printf("\t %.2lf", x[i]);
    }
    printf("\t]");

    printf("\ndet(A) = %.2lf", determinant);

    printf("\nr = [");
    for (int i = 0; i < rows; i++) {
        printf("\t %.20lf", residuals[i]);
    }
    printf("\t]");

    //запись результатов в файл
    std::fstream f;
    f.open(argv[2],  std::fstream::out);
    for (int i = 0; i < rows; i++) {
        f << x[i] << " ";
    }
    f << std::endl;

    for (int i = 0; i < rows; i++) {
        f << residuals[i] << " ";
    }
    f << std::endl;
    f << determinant;
    f.close();
}
