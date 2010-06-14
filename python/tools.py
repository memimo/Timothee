'''This module contains some of the frequently used functions'''

from config import *
from numpy import *

def get_onoff_map(idx, common):
    '''Read the recorded layer lgn map for npy files'''
   
    tmp = common['on_off_map'][idx]
    
    on_map = maximum(0, tmp)
    off_map = - minimum(0, tmp)
    
    size_var = shape(tmp)
    tmp2 = zeros((size_var[0], size_var[1], 2), float64)
    tmp2[:, :, 0] = on_map
    tmp2[:, :, 1] = off_map

    return tmp2


def get_s1_map(idx, common):
    '''Read the recorded layer s1 map for npy files'''

    global_idx = mod(idx, common['n_map']) + 1
    block = ceil(global_idx / (params['c']['block']))

   
    if (not common.has_key('current_block')) or block != \
    common['current_block']:
        s1_map = load(paths['map'] + 's1.' + str(int(block - 1)) + '.' + \
                      params['s']['type'] + '.npy')
        common['s1_map'] = s1_map
        common['current_block'] = block 

    idx_in_block = global_idx - (block - 1) * params['c']['block'] -1
 
    return common['s1_map'][:, :, :, idx_in_block], common


def dot_product(vect1, vect2):
    '''custom defined dot product, which just consider the first dimenstion'''

    n_num = min(len(vect1), len(vect2))
    return dot(vect1[:n_num], vect2[:n_num])

def dot_product2(vect1, vect2):
    '''normal dot product'''

    return dot(vect1[:, 0, 0], vect2[:, 0, 0])


def multi_dimensional_max(in_mat):
    '''Find the max and it's index in each dimenstion of a multi-dimenstional
    matrix'''

    size_var = shape(in_mat)
    n_dim = len(size_var)
    idx = zeros(n_dim)

    modulo = zeros((n_dim), float64)
    modulo[0] = 1

    for d_ind in range(n_dim -1):
        modulo[d_ind + 1] = modulo[d_ind] * size_var[d_ind]

    val = in_mat.max()
    rest = in_mat.flatten(1).argmax() + 1

    for d_ind in range(n_dim-1, -1, -1):
        
        idx[d_ind] = floor((rest-1) / (modulo[d_ind])) + 1
        rest = rest - (idx[d_ind] -1) * modulo[d_ind]

    return val, idx - 1


def permute(num):
    '''custom permute for random seeding reason'''

    return random.random(num).argsort()


