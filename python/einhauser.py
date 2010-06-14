from config import *
from numpy import *
import pdb
from tools import *

def einhauser(n_iter, n_feat, n_node, patch_size, weight, param, common):
    ''' implemention of Einhauser 2005'''

    if params['test_mode']:
        random.mtrand.seed(5489)

    param['nu'] = inf
    param['learning_rate'] = 5 * 10 ** -4
    param['thr_decay'] = 2 ** -15

    if len(weight) == 0:
        weight = zeros((patch_size[0], patch_size[1], n_feat, n_node), float64)
        for n_ind in range(n_node):
            #format: i x j x feat x node
            weight[:, :, :, n_ind] = get_initial_weight(patch_size, n_feat)

        #init middle layer
        if param['nu'] < inf:
            tr_am = 1.0 / param['nu'] * ones((patch_size[0], patch_size[1],
                                             n_feat), float64)
        else:
            tr_am = ones((patch_size[0], patch_size[1], n_feat), float64)

        threshold_am = .001 * ones((patch_size[0], patch_size[1], n_feat),
                                   float64)
        previous_max_am = 0.0
        previous_winner_mi = 0.0
        previous_winner_mj = 0.0
        previous_winner_mf = 0.0

        #init top layer
        if param['nu'] < inf:
            tr_at = 1.0 / param['nu'] * ones(n_node)
        else:
            tr_at = ones((n_node), float64)


        previous_max_at = -1
        previous_winner_t = -1

        #reporting
        n_firing_output = zeros((n_node), float64)
        n_firing_input = zeros((patch_size[0], patch_size[1], n_feat), float64)

        evol = zeros((param['n'], floor(param['n_iter'] / 100)), float64)

    print 'Einhauser 2002...'

    offset = floor(random.random())

    #Main iteration
    for i_ind in range(n_iter):
        x_var, common = get_s1_map(i_ind + offset, common)

        am_var = x_var / tr_am

        max_am = am_var.max()
        winner_m = float((am_var.flatten(1)).argmax())

        winner_mf = floor((winner_m) / patch_size[0] / patch_size[1]) 
        winner_mj = floor((winner_m + 1 - (winner_mf) * patch_size[0] * \
                          patch_size[1] -1.0) / patch_size[1])
        winner_mi = winner_m + 1 - (winner_mf) * patch_size[0] * patch_size[1] \
            - (winner_mj) * patch_size[0] - 1

        at_var = zeros((n_node), float64)
        for n_ind in range(n_node):
            at_var[n_ind] = (weight[:, :, :, n_ind] * am_var).max() / \
            tr_at[n_ind]

        #find the winner
        max_at = at_var.max()
        winner_t = at_var.argmax()

        #learning
        if i_ind > 0 and max_am > threshold_am[winner_mi, winner_mj, winner_mf]:
            threshold_am[winner_mi, winner_mj, winner_mf] = max_am
            weight[:, :, :, previous_winner_t] = (1.0 - param['learning_rate'])\
            * weight[:, :, :, previous_winner_t]
            weight[winner_mi, winner_mj, winner_mf, previous_winner_t] = \
            weight[winner_mi, winner_mj, winner_mf, previous_winner_t] + \
                param['learning_rate']

            n_firing_input[previous_winner_mi, previous_winner_mj,
                previous_winner_mf] = n_firing_input[previous_winner_mi,
                previous_winner_mj, previous_winner_mf] + 1
            n_firing_output[winner_t] = n_firing_output[winner_t] + 1


        if mod(i_ind + 1, 100) == 0:
            evol[:, int((i_ind + 1) / 100) - 1] = weight.sum()

        #update trace
        tr_am = am_var / param['nu'] + (1 - 1.0 / param['nu']) * tr_am
        tr_at = at_var / param['nu'] + (1 - 1.0 / param['nu']) * tr_at

        #thr decay
        threshold_am = threshold_am * (1-param['thr_decay'])

        #previous = curretn
        previus_max_am = max_am
        previous_winner_mi = winner_mi
        previous_winner_mj = winner_mj
        previous_winner_mf = winner_mf
        previous_winner_t = winner_t
        previous_max_at = max_at

        if mod(i_ind, 10000) == 0:
            print '.',

    common['weight'] = weight 
    common['evol'] = evol 
    common['thr'] = threshold_am 
    common['n_firing_input'] = n_firing_input
    common['n_firing_output'] = n_firing_output 
    
    return common

def get_initial_weight(patch_size, n_feat):
    '''Return random weight matrix'''
    
    weight = 0.75 * (1.0 + 0.0 * (random.random((patch_size[0], patch_size[1],
                     n_feat))- 0.5))
    return weight 

