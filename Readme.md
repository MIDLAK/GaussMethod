## О программе
Программа решает СЛАУ вида Ax=b.

## Компиляция программы
### Windows
cl /EHsc \<filename\>.cpp

### Linux
g++ \<filename\>.cpp

## Пример входных и выходных данных
В качестве аргументов командной строки программа принимает 2 параметра, а именно имена входного и выходного файлов. Заранее выходной файл создавать не нужно.

gauss.exe input.txt output.txt

### input.txt
Последний столбец матрицы содержит вектор b.

5 7 6 5 23
7 10 8 7 32
6 8 10 9 33
5 7 9 10 31

### output.txt
x = 1 1 1 1 
r = 3.55271e-15 3.55271e-15 7.10543e-15 3.55271e-15 
d = -1
inverse: 
	68	-41	-17	10
	-41	25	10	-6
	-17	10	5	-3
	10	-6	-3	2
