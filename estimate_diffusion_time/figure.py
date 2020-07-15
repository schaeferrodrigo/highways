# -*- coding: utf-8 -*-
#===============================================================================
import numpy as np
import pandas as pd
#from main_file import *
#from data_frame import *
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D



data = pd.read_csv('diffusion_time.csv')
x = data['angle']
y = data['eps']
print x
print y


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for index , c , m  in [('t_1' , 'blue','o'), ('t_2' , 'red','^'),('t_3' ,'green','s'),('t_4' , 'k','d') ]:
    ax.scatter(x, y, data[index], c=c, marker=m)
for index , c , m  in [('t_1' , 'blue','o'), ('t_2' , 'red','^'),('t_3' ,'green','s'),('t_4' , 'k','d') ]:
    ax.scatter(0.5, y, (y**(-2))*2, c='red', marker=m)
ax.set_xlabel(r'$\alpha$')
ax.set_ylabel(r'$\varepsilon$')
plt.title('Total time of diffusion ')






plt.show()

# for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
#     xs = randrange(n, 23, 32)
#     ys = randrange(n, 0, 100)
#     zs = randrange(n, zlow, zhigh)
#     ax.scatter(xs, ys, zs, c=c, marker=m)
#
#
# plt.show()
