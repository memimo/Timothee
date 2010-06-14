from numpy import *
from config import *
from tools import *
import pdb

def foldiak1991(weight, param, common):
    '''Implementation of Foldiak 1991 '''

    if params['test_mode']:
        random.mtrand.seed(5489)

    param['nu'] = 5.0
    param['learning_rate'] = 0.02

    if len(weight) == 0:
        weight = zeros((param['rf_size'][0],param['rf_size'][1], param['n_feat'], param['n']), float64)
        for n_ind in range(param['n']):
            weight[:,:,:,n_ind] = get_initial_weight(param['rf_size'], param['n_feat'])
            
        tr_y = 1.0 / param['nu']*ones(param['n'], float64)
        evol = zeros((param['n'], floor(param['n_iter'] / 100.0)), float64)
    print 'Foldiak 1991...'

    offset = floor(random.random())
    for i_ind in range(param['n_iter']):
        x_var, common = get_s1_map(i_ind + offset, common)

        y_var = zeros(param['n'], float64)
        for n_ind in range(param['n']):
            y_var[n_ind] = (weight[:,:,:,n_ind]*x_var).sum()

        max_y = y_var.max()
        winner = y_var.argmax()

        weight[:,:,:,winner] = weight[:,:,:,winner] + param['learning_rate'] * tr_y[winner] * (x_var - weight[:,:,:,winner])

        if mod(i_ind + 1, 100) == 0:
            evol[ :, int((i_ind + 1)/100) - 1] = weight.sum()

        tr_y = (array(range(param['n'])) == winner).astype(float64) / param['nu'] + (1.0 - 1.0 / param['nu']) * tr_y

        if mod(i_ind, 1000) == 0:
            print '.',


    common['weight'] = weight 
    common['evol'] = evol 
    
    return common

def get_initial_weight(patch_size, n_feat):
    return 0.1 * array(random.random((patch_size[0], patch_size[1], n_feat)))
    
