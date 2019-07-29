# Based on pdb_query

# Import all of the files
# Get mean chain length
# Get max chain length
# Remove outliers if needed
# Observe the dictionary
# Choose 40 families and num select points.
# Choose num_select points


# Imports
from Bio import PDB
import os
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import random



# Gets list with [x,y,z] and num select amount of points.
# If amount of points is smaller then the list size -> randomly choose num select points from the list
# If amount of points is larger than the list size -> Take all of the points and 0 pad the difference.
# points array is list
# num_select is integer


def scatter_plot(x_vals,y_vals,z_vals):

    fig = pyplot.figure()
    ax = Axes3D(fig)
    ax.scatter(x_vals,y_vals,z_vals)
    pyplot.show()

    return


def randomly_choose_points(points_array , num_select):
    pdb_points = np.zeros((num_select,3))
    len_points = len(points_array)
    np_points_array = np.array(points_array)
    index = np.random.choice(range(len_points), np.min([len_points,num_select]), replace=False)
    np_points_chosen = np_points_array[index]
    pdb_points[0:len(np_points_chosen)] = np_points_chosen
    return np.expand_dims(pdb_points,axis=0)


train_example_num = 400
family_num = 40

# Need to do something for all of the families together
# Create a dict where's every key in the dict is PDB file and the values are x,y,z coordinates

# Families folder located in tryout
families_folder = 'chosen families'
# Get all dirs in the folder
entries = os.listdir(families_folder)
#Initiate dict
dict_families = {}
# For all of the families - files in the folder
for family_name in entries:
    print(family_name)
    dict_pdbs = {}
    # Get all PDBS from the folder
    pdbs = os.listdir(families_folder +'/' + family_name)
    pdb_index = 0
    pdbs_available = len(pdbs)
    indices = np.random.choice(range(pdbs_available), train_example_num, replace=False)
# Run over all of the PDBS
    for ind in indices:
        pdb = pdbs[ind]
#        print(pdb)
        # Get x,y,z coordinate of the atom
        parser = PDB.PDBParser()
        io = PDB.PDBIO()
        struct = parser.get_structure(pdb,families_folder +'/' + family_name + '/' + pdb)
        dict_pdbs[pdb] = []
        for model in struct:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        x,y,z = atom.get_coord()
                        dict_pdbs[pdb].append([x,y,z])
                        # print(x,y,z)
    dict_families[family_name] = dict_pdbs


# Having the dict, select amount of points and get coordinates
num_select = 512
X_train = np.zeros([train_example_num*family_num,num_select,3])
X_labels = np.zeros(train_example_num*family_num)


none_zero_dict = {}

# Set label initial value
label = 0
for family_key in dict_families.keys():
    non_zero_count = 0
    non_zero_mean = 0 
    pdb_xyz_concat = np.zeros([1,num_select,3])
    family_dict = dict_families[family_key]
    for key in dict_families[family_key].keys():
        len_pdb = len(family_dict[key])
        if  len_pdb < num_select:
            non_zero_count += 1
            non_zero_mean += len_pdb
#            print(key , len_pdb)
        pdb_xyz = randomly_choose_points(family_dict[key], num_select)
        pdb_xyz_concat = np.concatenate((pdb_xyz_concat, pdb_xyz), axis=0)
        # scatter_plot(pdb_xyz[0][:,0],pdb_xyz[0][:,1],pdb_xyz[0][:,2])
    family_data  = pdb_xyz_concat[1:,:,:]
    family_label = np.ones([family_data.shape[0]]) * label
    print(family_key)
    print(['None zero count',non_zero_count])
    if non_zero_count == 0:
        Mean_length = 512
    if non_zero_count > 0:
        Mean_length = (non_zero_mean*non_zero_count/train_example_num) + 512*(train_example_num-non_zero_count)/train_example_num
        
    print(['Mean length',Mean_length ])
    print()
#    print(families_folder +'/' + family_key + '/' + family_key + '.npy')

    X_train[label*train_example_num:(label+1)*train_example_num,:,:] = family_data
    X_labels[label*train_example_num:(label+1)*train_example_num] = family_label
    label = label + 1
    none_zero_dict[family_key] = [non_zero_count,Mean_length]



np.save('X_train'  + '.npy', X_train) # save
np.save('X_label'  + '.npy', X_labels) # save

f = open("none_zero_dict.txt","w")
f.write( str(none_zero_dict) )
f.close()