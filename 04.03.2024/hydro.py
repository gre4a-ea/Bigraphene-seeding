import os 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import sys
import subprocess
import math
import csv
from datetime import datetime


def layer_creator(border, alpha, basis, a, mixshift, height):
    A=np.array([[math.cos(math.pi/180*alpha), -math.sin(math.pi/180*alpha)],[math.sin(math.pi/180*alpha), math.cos(math.pi/180*alpha)]])
    coordinates=np.vstack(((np.arange(-border, border, 1).reshape(-1,1)*np.ones(2*border)).T.reshape(1, 4*border**2), (np.arange(-border, border, 1).reshape(-1,1)*np.ones(2*border)).reshape(1, 4*border**2)))
    shift=np.array([a/math.sqrt(3),0])
    data_3=np.hstack(   (np.dot(A,(np.dot(basis, coordinates).T + shift+mixshift).T).T,  height*np.ones(((2*border)**2, 1)))    ) 
    data_4=np.hstack(   (np.dot(A,(np.dot(basis, coordinates).T+mixshift).T).T,  height*np.ones(((2*border)**2, 1)))    )  
    return  np.vstack((  data_3, data_4))

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

def rename_existing_file(filename):
    if os.path.exists(filename):
        # Получаем текущую дату и время
        now = datetime.now()
        # Создаем строку с текущим временем и датой в формате ГГГГММДД_ЧЧММСС
        timestamp = now.strftime("%Y%m%d_%H%M%S")
        # Получаем расширение файла
        file_extension = os.path.splitext(filename)[1]
        # Генерируем новое имя файла с добавлением времени и даты
        new_filename = f"{timestamp}{file_extension}"
        # Переименовываем существующий файл
        new_filename = os.path.join(ROOT_DIR, new_filename)

        os.rename(filename, new_filename)
        print(f"Файл {filename} переименован в {new_filename}")
    else:
        print(f"Файл {filename} не найден")

class bigraphene:
    def __init__(self, alpha, border, a, mixshift, thickness):
        self.alpha = alpha
        self.border = border
        self.a = a
        self.mixshift = mixshift
        self.thickness = thickness
        radius = border*a*math.sqrt(3)/2
        self.radius = radius

        basis=np.array([[ math.sqrt(3)/2*a, 0],[-a/2, a]])
        carbon_low= layer_creator(border, 0, basis, a, 0, -thickness/2)
        carbon_high= layer_creator(border, alpha, basis, a, 0, thickness/2)
        
        
        carbon_high=carbon_high[carbon_high[:,0]**2+ (carbon_high[:,1])**2 <= (radius)**2 ]+mixshift
        carbon_low=carbon_low[carbon_low[:,0]**2+ carbon_low[:,1]**2 <= (radius)**2 ]
 
        self.carbon_high = carbon_high
        self.carbon_low = carbon_low
        self.activated_index_high = np.full(np.shape(carbon_high)[0], False, dtype=bool).reshape(-1,1)
        self.activated_index_low = np.full(np.shape(carbon_low)[0], False, dtype=bool).reshape(-1,1)
        self.hydrogen_high=np.empty((0, 3))
        self.hydrogen_low=np.empty((0, 3))
        self.energy=None
        
    def AddHydrogenNearest(self, bording_distance):


        
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_low)
        distances_high, nearest_neighbor_index = nbrs.kneighbors(self.carbon_high)

        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_high)
        distances_low, nearest_neighbor_index = nbrs.kneighbors(self.carbon_low)



    
        adding_indexes_low =np.logical_and( (distances_low<=bording_distance), (self.carbon_low[:,0]**2+ self.carbon_low[:,1]**2 <= (self.radius)**2 ).reshape(-1,1))
        adding_indexes_high=np.logical_and( (distances_high<=bording_distance), (self.carbon_high[:,0]**2+ self.carbon_high[:,1]**2 <= (self.radius)**2 ).reshape(-1,1))

        if np.any(adding_indexes_low) or np.any(adding_indexes_high):

            nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(self.carbon_low)
            distances, nearest_neighbor_index = nbrs.kneighbors(self.carbon_low[adding_indexes_low[:,0]])
            mask_low = np.zeros(np.shape(self.carbon_low)[0], dtype=bool).reshape(-1,1)
            mask_low[nearest_neighbor_index[:,1]]=True
            mask_low[nearest_neighbor_index[:,2]]=True
            mask_low[nearest_neighbor_index[:,3]]=True

            nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(self.carbon_high)
            distances, nearest_neighbor_index = nbrs.kneighbors(self.carbon_high[adding_indexes_high[:,0]])
            mask_high = np.zeros(np.shape(self.carbon_high)[0], dtype=bool).reshape(-1,1)
            mask_high[nearest_neighbor_index[:,1]]=True
            mask_high[nearest_neighbor_index[:,2]]=True
            mask_high[nearest_neighbor_index[:,3]]=True
        else:
            return self

        


        mask_low=np.logical_and( mask_low, ~self.activated_index_low, ~adding_indexes_low)
        mask_high=np.logical_and(mask_high, ~self.activated_index_high, ~adding_indexes_high)
              


        hydrogen_low=self.carbon_low[np.any(mask_low, axis=1)]-np.array([0,0,1.5])
        hydrogen_high=self.carbon_high[np.any(mask_high, axis=1)]+np.array([0,0,1.5])

        self.hydrogen_low=np.vstack((self.hydrogen_low, hydrogen_low))
        self.hydrogen_high=np.vstack((self.hydrogen_high, hydrogen_high))

    
        self.carbon_low[np.any(mask_low, axis=1)]-=np.array([0,0,0.1])
        self.carbon_high[np.any(mask_high, axis=1)]+=np.array([0,0,0.1])

        self.activated_index_high=np.logical_or(mask_high, self.activated_index_high)
        self.activated_index_low=np.logical_or(mask_low, self.activated_index_low)

        self.activated_index_high=np.logical_or(adding_indexes_high, self.activated_index_high)
        self.activated_index_low=np.logical_or(adding_indexes_low, self.activated_index_low)
    
    def AddHydrogenPartly(self, part, border):


        
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_low)
        distances_high, nearest_neighbor_index = nbrs.kneighbors(self.carbon_high)

        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_high)
        distances_low, nearest_neighbor_index = nbrs.kneighbors(self.carbon_low)

    

    
        adding_indexes_low = (self.carbon_low[:,0]**2+ self.carbon_low[:,1]**2 <= (self.radius*border)**2 ).reshape(-1,1)
        adding_indexes_high= (self.carbon_high[:,0]**2+ self.carbon_high[:,1]**2 <= (self.radius*border)**2 ).reshape(-1,1)

        mask_low=np.logical_and( ~self.activated_index_low, ~adding_indexes_low)
        mask_high=np.logical_and( ~self.activated_index_high, ~adding_indexes_high)

        distances_high_ = np.sort(distances_high[mask_high])
        distances_low_ = np.sort(distances_low[mask_low])

        print(np.shape(distances_high_)[0])

        bording_distance=max(distances_high_[int(np.shape(distances_high_)[0]*part)], distances_low_[int(np.shape(distances_low_)[0]*part)])
        

        adding_indexes_low =np.logical_and( (distances_low<=bording_distance), (self.carbon_low[:,0]**2+ self.carbon_low[:,1]**2 <= (self.radius*border)**2 ).reshape(-1,1))
        adding_indexes_high=np.logical_and( (distances_high<=bording_distance), (self.carbon_high[:,0]**2+ self.carbon_high[:,1]**2 <= (self.radius*border)**2 ).reshape(-1,1))

        if np.any(adding_indexes_low) or np.any(adding_indexes_high):

            nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(self.carbon_low)
            distances, nearest_neighbor_index = nbrs.kneighbors(self.carbon_low[adding_indexes_low[:,0]])
            mask_low = np.zeros(np.shape(self.carbon_low)[0], dtype=bool).reshape(-1,1)
            mask_low[nearest_neighbor_index[:,1]]=True
            mask_low[nearest_neighbor_index[:,2]]=True
            mask_low[nearest_neighbor_index[:,3]]=True

            nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(self.carbon_high)
            distances, nearest_neighbor_index = nbrs.kneighbors(self.carbon_high[adding_indexes_high[:,0]])
            mask_high = np.zeros(np.shape(self.carbon_high)[0], dtype=bool).reshape(-1,1)
            mask_high[nearest_neighbor_index[:,1]]=True
            mask_high[nearest_neighbor_index[:,2]]=True
            mask_high[nearest_neighbor_index[:,3]]=True
        else:
            return self

        


        mask_low=np.logical_and( mask_low,  ~adding_indexes_low)
        mask_high=np.logical_and(mask_high, ~adding_indexes_high)

        mask_low=np.logical_and( mask_low, ~self.activated_index_low)
        mask_high=np.logical_and(mask_high, ~self.activated_index_high)
              


        hydrogen_low=self.carbon_low[np.any(mask_low, axis=1)]-np.array([0,0,1.1])
        hydrogen_high=self.carbon_high[np.any(mask_high, axis=1)]+np.array([0,0,1.1])

        self.hydrogen_low=np.vstack((self.hydrogen_low, hydrogen_low))
        self.hydrogen_high=np.vstack((self.hydrogen_high, hydrogen_high))

    
        self.carbon_low[np.any(mask_low, axis=1)]-=np.array([0,0,0.1])
        self.carbon_high[np.any(mask_high, axis=1)]+=np.array([0,0,0.1])

        self.activated_index_high=np.logical_or(mask_high, self.activated_index_high)
        self.activated_index_low=np.logical_or(mask_low, self.activated_index_low)

        self.activated_index_high=np.logical_or(adding_indexes_high, self.activated_index_high)
        self.activated_index_low=np.logical_or(adding_indexes_low, self.activated_index_low)

    
        
        # nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_low)
        # distances_high, nearest_neighbor_index = nbrs.kneighbors(self.carbon_high)

        # nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_high)
        # distances_low, nearest_neighbor_index = nbrs.kneighbors(self.carbon_low)


        # mask_low =np.logical_and( (distances_low>=bording_distance), (self.carbon_low[:,0]**2+ self.carbon_low[:,1]**2 <= (self.radius*border)**2 ).reshape(-1,1))
        # mask_high = np.logical_and( (distances_high>=bording_distance), (self.carbon_high[:,0]**2+ self.carbon_high[:,1]**2 <= (self.radius*border)**2 ).reshape(-1,1))

        # mask_low=np.array(np.logical_and( mask_low, ~self.activated_index_low))
        # mask_high=np.logical_and(mask_high, ~self.activated_index_high)

        # hydrogen_low=self.carbon_low[np.any(mask_low, axis=1)]-np.array([0,0,1.1])
        # hydrogen_high=self.carbon_high[np.any(mask_high, axis=1)]+np.array([0,0,1.1])

        # self.hydrogen_low=np.vstack((self.hydrogen_low, hydrogen_low))
        # self.hydrogen_high=np.vstack((self.hydrogen_high, hydrogen_high))

        # self.carbon_low[np.any(mask_low, axis=1)]-=np.array([0,0,0.1])
        # self.carbon_high[np.any(mask_high, axis=1)]+=np.array([0,0,0.1])



    def RemoveH2(self, border):
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(self.hydrogen_high)
        distances_high, index_high = nbrs.kneighbors(self.hydrogen_high)
        mask_high = (index_high[:,1]== index_high[index_high[index_high[:,1]][:,1]][:,1])
        self.hydrogen_high=self.hydrogen_high[~np.logical_and(distances_high[:,1]<border, mask_high)]
        
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(self.hydrogen_low)
        distances_low, index_low = nbrs.kneighbors(self.hydrogen_low)
        mask_low = (index_low[:,1]== index_low[index_low[index_low[:,1]][:,1]][:,1])
        self.hydrogen_low=self.hydrogen_low[~np.logical_and(distances_low[:,1]<border, mask_low)]

    def MoveToCsv(self, file_adress):
        new_c = pd.DataFrame(np.hstack((self.carbon_low, self.carbon_high)), columns=['x_c_l', 'y_c_l', 'z_c_l', 'x_c_h', 'y_c_h', 'z_c_h'])
        new_h = pd.DataFrame(np.hstack((self.hydrogen_low, self.hydrogen_high)), columns=['x_h_l', 'y_h_l', 'z_h_l', 'x_h_h', 'y_h_h', 'z_h_h'])
        

        if os.path.exists(file_adress):
            with open(file_adress, 'a') as f:
                new_c.to_csv(f)
                new_h.to_csv(f)
        else:
            new_c.to_csv(file_adress)
            new_h.to_csv(file_adress)

    def MoveToCSVLogs(self, file_adress, info):
        new_c = self.energy
        if os.path.exists(file_adress):
            with open(file_adress, mode='a', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([self.energy, (np.shape(self.hydrogen_low)[0]+np.shape(self.hydrogen_high)[0])/(np.shape(self.carbon_low)[0]+np.shape(self.carbon_high)[0]), info])
        else:
            with open(file_adress, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow([self.energy, (np.shape(self.hydrogen_low)[0]+np.shape(self.hydrogen_high)[0])/(np.shape(self.carbon_low)[0]+np.shape(self.carbon_high)[0]), info])

    def DumpToXYZ(self, file_adress):
        atom1_positions=np.vstack((self.carbon_high,self.carbon_low))
        atom2_positions=np.vstack((self.hydrogen_high,self.hydrogen_low))

        if os.path.exists(file_adress):
            num_atoms = len(atom1_positions) + len(atom2_positions)
            with open(file_adress, 'a') as f:
                f.write(f"{num_atoms}\n")
                f.write("Generated XYZ file\n")

                for pos in atom1_positions:
                    f.write(f"AtomType1 {pos[0]} {pos[1]} {pos[2]}\n")

                for pos in atom2_positions:
                    f.write(f"AtomType2 {pos[0]} {pos[1]} {pos[2]}\n")
        else:
            num_atoms = len(atom1_positions) + len(atom2_positions)

            with open(file_adress, 'w') as f:
                f.write(f"{num_atoms}\n")
                f.write("Generated XYZ file\n")

                for pos in atom1_positions:
                    f.write(f"AtomType1 {pos[0]} {pos[1]} {pos[2]}\n")

                for pos in atom2_positions:
                    f.write(f"AtomType2 {pos[0]} {pos[1]} {pos[2]}\n")

    def ReadEnergy(self, file_adress):
        with open(file_adress, 'r') as file:
            lines = file.readlines()
        self.energy = float(lines[-2].strip())

    def MoveToCoo(self, root_address):
        carbon= np.vstack((self.carbon_low, self.carbon_high))
        hydrogen=np.vstack((self.hydrogen_low, self.hydrogen_high))

        ones_carbon = 1*np.ones((carbon.shape[0], 1), dtype='int64')
        carbon = np.append(ones_carbon, carbon,  axis=1)

        ones_hydrogen = 2*np.ones((hydrogen.shape[0], 1), dtype='int64')
        hydrogen = np.append(ones_hydrogen ,hydrogen,  axis=1)

        positions = np.vstack((hydrogen, carbon))
        natoms=np.shape(positions)[0]
        row_numbers = np.arange(positions.shape[0]).reshape(-1, 1)
        positions = np.hstack((row_numbers, positions))
                
        with open(root_address+'/coo','w') as fdata:
            # First line is a comment line 
            fdata.write('Mixed graphene (30 degrees)\n')
            #--- Header ---#
            # Specify number of atoms and atom types 
            fdata.write('{} atoms\n'.format(natoms))
            fdata.write('{} atom types\n\n'.format(2))
            # Specify box dimensions
            fdata.write('{} {} xlo xhi\n'.format(-self.radius, self.radius))
            fdata.write('{} {} ylo yhi\n'.format(-self.radius, self.radius))
            fdata.write('{} {} zlo zhi\n'.format(-self.radius, self.radius))
            fdata.write('\n')
            
            # Atoms section
            fdata.write('Masses\n\n')
            fdata.write('1 12.0107\n')
            fdata.write('2 1.001\n')
            fdata.write('\n')
            fdata.write('Atoms\n\n')
            # Write each position 
            for i, atom_type ,pos_x, pos_y, pos_z in positions:
                fdata.write('{} {} {} {} {}\n'.format(int(i+1), int(atom_type), pos_x, pos_y, pos_z))
    
    def ReadXYZ(self, file_address):
        numpy_matrices=xyz_to_dataframe(file_address)
        result_matrix=numpy_matrices[-1]
        hydrogen=result_matrix[result_matrix[:,0]=='H'][:,1:4].astype(float)
        carbon=result_matrix[result_matrix[:,0]=='C'][:,1:4].astype(float)


        self.carbon_low=carbon[carbon[:,2]<0]
        self.carbon_high=carbon[carbon[:,2]>0]

        self.hydrogen_low=hydrogen[hydrogen[:,2]<0]
        self.hydrogen_high=hydrogen[hydrogen[:,2]>0]

        self.activated_index_high = np.full(np.shape(self.carbon_high)[0], False, dtype=bool).reshape(-1,1)
        self.activated_index_low = np.full(np.shape(self.carbon_low)[0], False, dtype=bool).reshape(-1,1)
        
        if np.shape(self.hydrogen_low)[0] != 0:
            nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_low)
            distances_high, nearest_neighbor_index = nbrs.kneighbors(self.hydrogen_low)
            self.activated_index_low[nearest_neighbor_index[:,0]] = True
        else:
            print("No self.hydrogen_low")
        if np.shape(self.hydrogen_high)[0] != 0:
            nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_high)
            distances_low, nearest_neighbor_index = nbrs.kneighbors(self.hydrogen_high)
            self.activated_index_high[nearest_neighbor_index[:,0]] = True
        else:
            print("No self.hydrogen_high")
        
        
        nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(self.carbon_low[:,0:2])
        distances_high, nearest_neighbor_index = nbrs.kneighbors(self.carbon_low[self.activated_index_low[:,0]][:,0:2])
                
        nni=np.hstack((nearest_neighbor_index[:,1],nearest_neighbor_index[:,2],nearest_neighbor_index[:,3]))
        self.activated_index_low[np.where(np.bincount(nni)>=2)[0]] = True

        nbrs = NearestNeighbors(n_neighbors=4, algorithm='ball_tree').fit(self.carbon_high[:,0:2])
        distances_high, nearest_neighbor_index = nbrs.kneighbors(self.carbon_high[self.activated_index_high[:,0]][:,0:2])
                
        nni=np.hstack((nearest_neighbor_index[:,1],nearest_neighbor_index[:,2],nearest_neighbor_index[:,3]))
        self.activated_index_high[np.where(np.bincount(nni)>=2)[0]] = True

    def RemoveHydrogen(self, border):
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_high)
        distances_high, index_high = nbrs.kneighbors(self.hydrogen_high)
        self.hydrogen_high=self.hydrogen_high[distances_high[:,0]<border]
        
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(self.carbon_low)
        distances_low, index_high = nbrs.kneighbors(self.hydrogen_low)
        self.hydrogen_low=self.hydrogen_low[distances_low[:,0]<border]

    def ShowPlot(self):
        fig = plt.figure(figsize = (10, 10))
        ax = plt.axes(projection ="3d")
        
        ax.scatter3D(self.hydrogen_low[0], self.hydrogen_low[1], self.hydrogen_low[2], s=7, alpha = 0.5, color = "red")
        

        plt.title("simple 3D scatter plot")
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.legend()
        plt.show() 

alpha=30
border=50
a=2.46
mixshift=[0.1,0.1,0]
thickness=3.17
actual_part=1
excication_part=0.02


ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) # This is your Project Root
sys.path.append(ROOT_DIR)



    


x=bigraphene(alpha, border, a, mixshift, thickness)
x.carbon_high=None
x.carbon_low=None
x.hydrogen_high=None
x.hydrogen_low=[0,0,0]
print(x)
x.ShowPlot()
x.MoveToCoo(ROOT_DIR)
return_code = subprocess.call("mpirun -np 6 --use-hwthread-cpus ~/soft/exe/lmp_mlip < input | tee out", shell=True)
if return_code!=0:
    print("Модель сколлапсировала на шаге "+str(2*i))
    exit()
else: 
    print("Модель оптимизирована на шаге "+str(2*i))        
x.ReadEnergy(ROOT_DIR+"/out")
        
        