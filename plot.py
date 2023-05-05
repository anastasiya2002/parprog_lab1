import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Загрузка данных из CSV файла
dataframe = pd.read_csv('build/output.csv')
x = dataframe['x']
t = dataframe['t']
u = dataframe['u']

# Создание 3D графика
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(x, t, u, cmap='viridis')

# Настройка осей
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u')
ax.set_title('График функции u(x,t)')

# Отображение графика
plt.show()




data = pd.read_csv('build/time.csv')
n = data['n']
time = data['time']

plt.plot(n,time[0] / (time))
plt.xlabel('n')
plt.ylabel('S')
plt.title('Ускорение S(n)')
plt.show()

data = pd.read_csv('build/time.csv')
n = data['n']
time = data['time']

plt.plot(n,time[0] / (n * time))
plt.xlabel('n')
plt.ylabel('S')
plt.title('Эффективность E(n)')
plt.show()
