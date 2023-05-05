# Лабораторная работа №1
Численное решение уравнения переноса с использованием технологии MPI

Изначально в папке находятся файлы:

func.c - описание функций и шагов

main.c - программа с распараллеливанием

Makefile - файл для сборки

plot.py - построение графиков в Python3

Запуск:

$ make build

$ cd build

/build$ mpirun -n (number of processes) -oversubscribe task

В папке build появляются 2 .csv файла для графиков. Чтобы графики Эффективности и Ускорения были презентабельными, стоит последнюю команду ввести несколько раз подряд, меняя (number of processes) последовательно от 1 до (сколько угодно). 

Для графика ускорения в файле time.csv вручную нужно добавить строчку в самое начало: n,time

Запуск построения графиков:

$ python3 plot.py

Удаление результатов вычисления: $ make clean
