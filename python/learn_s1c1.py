'''Main module to run the whole program'''


import os.path

from config import *
import layer_c1
import layer_lgn
import layer_s1


def l_lgn():
    '''Compute the lgn layer'''

    #Check if it is already computed
    exists = zeros((params['c']['rf_size'][0], params['c']['rf_size'][1]), bool)
    for i_ind in range(params['c']['rf_size'][0]):
        for j_ind in range(params['c']['rf_size'][1]):
            if os.path.isfile(paths['map'] + 'onoff.' + str(i_ind) + '.' + \
                              str(j_ind) + '.npy'):
                exists[i_ind, j_ind] = True

    #If not then compute
    if not exists.all():
        layer_lgn.compute_layer()


def l_s1():
    '''Computer Layer S1 of the network'''

    tmp_lgn = load(paths['map'] + 'onoff.0.0.npy')
    count = ceil(float(len(tmp_lgn)) / params['c']['block'])
    if count == 0:
        count = 1
    exists = zeros((count), bool)
    
    for block in range(count):
        if os.path.isfile(paths['map'] + 's1.' + str(block) + '.' + \
                          params['s']['type'] + '.npy'):
            exists[block] = True        
            
    #If not then compute
    if not exists.all():
        layer_s1.layer_s1()



def l_c1():
    '''Computer Layer C1 of the network'''

    #Timothee    
    if not os.path.isfile(paths['map'] + 'common.c_timothee.' + \
                          params['comp_name'] + '.' + \
                          params['s']['type'] + '.pck'):
        layer_c1.learn_c1('timothee', {})

    #Foldiak
    if not os.path.isfile(paths['map'] + 'common.c_foldiak.' + \
                          params['comp_name'] + '.' + \
                          params['s']['type'] + '.pck'):
        layer_c1.learn_c1('foldiak', {})

    #Einhasuer
    if not os.path.isfile(paths['map'] + 'common.c_einhauser.' + \
                          params['comp_name'] + '.' + \
                          params['s']['type'] + '.pck'):
        layer_c1.learn_c1('einhasuer', {})

    #Spratling
    if not os.path.isfile(paths['map'] + 'common.c_spratling.' + \
                          params['comp_name'] + '.' + \
                          params['s']['type'] + '.pck'):
        layer_c1.learn_c1('spratling', {})



def run_all():
    '''Run the three layers consequently'''
    l_lgn()
    l_s1()
    l_c1()


if __name__ == "__main__":
    run_all()
