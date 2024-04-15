import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import scipy.fft
import torch.nn as nn
import torch
from sklearn.neighbors import NearestNeighbors
from matplotlib.animation import FuncAnimation
import os 
import sys
import imageio

ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) # This is your Project Root
sys.path.append(ROOT_DIR)

def xyz_to_dataframe(adress):
    with open(str(adress), 'r') as file:
        lines = file.readlines()
    df = pd.DataFrame([line.split() for line in lines], columns=["Element", "X", "Y", "Z"])
    

    # Получение индексов разделительных строк
    separator_indices = df.index[df['X'].isnull()].tolist()

    # Разбиение DataFrame на отдельные части (NumPy матрицы) с помощью цикла
    matrix_chunks = [df.iloc[separator_indices[i-1]+2:separator_indices[i]] for i in range(1, len(separator_indices))]

    # Преобразование каждой части в отдельную NumPy матрицу
    return [chunk[['Element', 'X', 'Y', 'Z']].to_numpy() for chunk in matrix_chunks]




def make_table(file_name):
    return pd.read_csv(file_name, sep=',',header=None, names=None, index_col=None).values


def create_histogram(frame):
    fig, ax = plt.subplots()
    numpy_matrices=xyz_to_dataframe(ROOT_DIR +'/logs.xyz')
    result_matrix=numpy_matrices[frame*5]
    hydrogen=result_matrix[result_matrix[:,0]=='AtomType2'][:,1:4].astype(float)
    carbon=result_matrix[result_matrix[:,0]=='AtomType1'][:,1:4].astype(float)
    carbon_low=carbon[carbon[:,2]<0]
    carbon_high=carbon[carbon[:,2]>0]
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(carbon_low)
    distances, nearest_neighbor_index = nbrs.kneighbors(carbon_high)
    data=distances[:,0]  # новые данные

    ax.hist(data, bins=100, density=True)
    ax.set_title('Boundary length hystogram. Straight algoruthm. Frame {}'.format(frame*5))
    
    # Сохранение изображения
    filename = 'frame_{}.png'.format(frame)
    plt.savefig(filename)
    plt.close()




# frames = 17 # Количество кадров в анимации

# # Генерация гистограмм для каждого кадра и сохранение их как изображения
# for frame in range(frames):
#     create_histogram(frame)

# # Создание гиф-анимации из сохраненных изображений
# images = []
# for frame in range(frames):
#     filename = 'frame_{}.png'.format(frame)
#     print("Frame "+str(frame)+" created")
#     images.append(imageio.imread(filename))

# # Сохранение гиф-анимации
# imageio.mimsave('hist_evo2.gif', images, duration=100)

# # Удаление временных файлов изображений
# import os
# for frame in range(frames):
#     filename = 'frame_{}.png'.format(frame)
#     os.remove(filename)
    








plt.figure(figsize=(16, 9))
# plt.xticks(np.arange(-1, 21, step=0.2), minor = True)
# plt.yticks(np.arange(-1, 5, step=0.02), minor = True)
# plt.xticks(np.arange(-1, 21, step=1))
# plt.yticks(np.arange(-1, 5, step=0.1))
# plt.grid(which='minor', alpha=0.2)
# plt.grid(which='major', alpha=1)

color=["green","blue","red", "purple",  "orange","yellow","grey", "black"]
subtext=str("M Ce(SO4)2, KBrO3, CH2(COOH)2:   " )
label=["0.001, 0.06, 0.03","0.001, 0.06, 0.3","0.001, 0.06, 1.2","0.001, 0.06, 0.8", "0.001, 0.06, 1.5","0.001, 0.1, 1.2","0.001, 0.01, 1.2","8"]


actual_part=1
excication_part=0.02
data3=np.arange(100, dtype=float)

# numpy_matrices=xyz_to_dataframe(ROOT_DIR +'/logs.xyz')[::4]



for j in range(1):

    

    data1=xyz_to_dataframe('logs_m_'+str(j)+'.xyz')
    # data2=xyz_to_dataframe('logs_s_'+str(j)+'.xyz')
    data1_1=make_table('logs_m_'+str(j)+'.csv')
    # data2_1=make_table('logs_s_'+str(j)+'.csv')
    
    
    E_h=-1.12157786283786


    E_bg=-252214
    E_t=data1_1[:,0]
    x=np.zeros((150))
    y=np.zeros((150))
    for i in range(150):
        result_matrix=data1[i]
        hydrogen=result_matrix[result_matrix[:,0]=='AtomType2'][:,1:4].astype(float)
        # carbon=result_matrix[result_matrix[:,0]=='AtomType1'][:,1:4].astype(float)
        # carbon_low=carbon[carbon[:,2]<0]
        # carbon_high=carbon[carbon[:,2]>0]
        # nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(carbon_low)
        # distances, nearest_neighbor_index = nbrs.kneighbors(carbon_high)
        n=np.shape(hydrogen)[0]
        x[i]=n
        y[i]=(E_t[i]-E_bg-n*E_h)/n

        

        print("dot "+str(i)+ " plotted")


    print("plotted")
    plt.plot(x, y, color='red')

    # E_bg=-252214
    # E_t=data2_1[::3,0]
    # x=np.zeros((100))
    # y=np.zeros((100))
    # for i in range(1,100):
    #     result_matrix=data2[4*i]
    #     hydrogen=result_matrix[result_matrix[:,0]=='AtomType2'][:,1:4].astype(float)
    #     # carbon=result_matrix[result_matrix[:,0]=='AtomType1'][:,1:4].astype(float)
    #     # carbon_low=carbon[carbon[:,2]<0]
    #     # carbon_high=carbon[carbon[:,2]>0]
    #     # nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(carbon_low)
    #     # distances, nearest_neighbor_index = nbrs.kneighbors(carbon_high)
    #     n=np.shape(hydrogen)[0]
    #     x[i]=n
    #     y[i]=(E_t[i]-E_bg-n*E_h)/n

    #     dx = x[i] - x[i-1]
    #     dy = y[i] - y[i-1]

    #     print("dot "+str(i)+ " plotted")



    #     plt.arrow( x[i-1], y[i-1],  dx, dy,  fc='blue', ec='blue', color = "blue")
    # print("plotted")
    # plt.scatter(x, y, marker='o', linestyle='None', color='blue')



plt.legend()
plt.ylabel(r"$specific E_ $", fontsize = 25)
plt.xlabel(r"$N$, number of hydrogen ", fontsize = 25)

# plt.suptitle("Evolution of system with different hydrogen seeding algorythms", fontsize=30, y=1.2)

plt.show() 



