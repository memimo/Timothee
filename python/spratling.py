from numpy import *
from config import *
from tools import *

def spratling2005(n_iter, n_feat, n_node, patch_size, weight, param, common):
    ''' implementation of Spratling 2005'''

    if params['test_mode']:
        random.mtrand.seed(5489)

    param['learning_rate'] = 0.001

    if len(weight) == 0:
        weight = zeros((patch_size[0],patch_size[1], n_feat, n_node), float64)
        for n_ind in range(n_node):
            weight[:,:,:,n_ind] = get_initial_weight(patch_size, n_feat)

        former_y = ones((n_node), float64)
        mean_former_y = former_y.mean()

        n_firing_output = zeros((n_node), float64)
        n_firing_input = zeros((patch_size[0], patch_size[1], n_feat), float64)

        evol = zeros((param['n'], floor(param['n_iter'] / 100.0)), float64)

    print 'Spratling 2005...'

    offset = floor(random.random())

    for i_ind in range(n_iter):
        x_var, common = get_s1_map(i_ind + offset, common)

        node_max = ((weight.max(0)).max(0)).max(0)
        node_max = node_max.flatten(1)
        input_max = weight.max(3)
        iput_max = input_max +  ~input_max.astype(bool)

        y_var = zeros((n_node), float64)

        for n_ind in range(n_node):
            if node_max[n_ind] == 0:
                raise NameError('node_max[n]==0')
            if any(input_max == 0):
                raise NameError('any(input_max == 0)')

            z_var = x_var  * (weight[:,:,:,n_ind]**2.0) / input_max / node_max[n_ind]
            y_var[n_ind] = z_var.max()
            winner = z_var.flatten(1).argmax()

            if former_y[n_ind] > mean_former_y:
                weight[:,:,:,n_ind] = weight[:,:,:,n_ind] - param['learning_rate'] * (former_y[n_ind] - mean_former_y) / (n_node * mean_former_y) * x_var
                winner_f = floor((winner) / patch_size[0] / patch_size[1]) 
                winner_j = floor((winner + 1 - (winner_f) * patch_size[0] * patch_size[1] - 1) / patch_size[0]) 
                winner_i = winner + 1 - (winner_f)*patch_size[0] * patch_size[1] - (winner_j) * patch_size[1] - 1
                weight[winner_i, winner_j, winner_f, n_ind] = weight[winner_i, winner_j, winner_f, n_ind] + 2.0 * param['learning_rate'] * x_var[winner_i, winner_j, winner_f] * (former_y[n_ind] - mean_former_y) / (n_node * mean_former_y)

                n_firing_input[winner_i, winner_j, winner_f] = n_firing_input[winner_i, winner_j, winner_f] + 1
                n_firing_output[n_ind] = n_firing_output[n_ind] + 1


        weight = maximum(weight, 0)
        norm_factor = sum(weight, 3)
        for n_ind in range(n_node):
            weight[:,:,:,n_ind] = weight[:,:,:,n_ind] / norm_factor

        
        if mod(i_ind + 1, 100) == 0:
            evol[ :, int((i_ind + 1)/100) - 1] = weight.sum()

        former_y = y_var
        mean_former_y = former_y.mean()

        if mod(i_ind, 10000) == 0:
            print '.',

    common['weight'] = weight 
    common['evol'] = evol 
    common['n_firing_input'] = n_firing_input
    common['n_firing_output'] = n_firing_output 
    
    return common

def get_initial_weight(patch_size, n_feat):
    weight = 1 + .01*random.random((patch_size[0], patch_size[1], n_feat))
    weight = maximum(0, weight)
    weight = weight / weight.sum()

    return weight

