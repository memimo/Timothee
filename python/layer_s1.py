from config import *
from learn_unlearn import *
from numpy import *
import pickle
from tools import *


def record_response(range_var, weight_sharig_mode, center, common):
    '''
    Timothee Masquelier timothee.masquelier@alum.mit.edu August 2006
    n_map: number of input maps
    get_map: handle to a function i -> map
    weight_sharing_mode: if true convolution mode + max response, if false simple local dot product
    center: structure that must contains at least center.patch
    '''
    n_center = shape(center)[2]
    
    patch = []
    for c_ind in range(shape(center)[2]):
        tmp = center[:, :, c_ind, :]
        patch.append(tmp.flatten(1))
          
    response = zeros((n_center, len(range_var)), float64)
    print 'Recording responces...'
    if weight_sharig_mode:
        raise NameError('this does not work yet')
    else:
        for m_ind in range_var:
            img = get_onoff_map(m_ind, common)
            img = img / linalg.norm(img.flatten(1))
            
            for cent_ind in range(n_center): 
                response[cent_ind, m_ind - range_var[0]] = dot(img.flatten(1),
                                                               patch[cent_ind])
      
    return response
    


def layer_s1():

    common = {}
    common['center'] = zeros((params['c']['rf_size'][0], 
                            params['c']['rf_size'][1],
                            params['s']['rf_size'][0],
                            params['s']['rf_size'][1], params['s']['n'], 2))
    for i_ind in range(params['c']['rf_size'][0]):
        for j_ind in range(params['c']['rf_size'][1]):
            print 'S1 MAP (' + str(i_ind) + ', ' + str(j_ind) + ')'
            
            common['on_off_map'] = load(paths['map'] + 'onoff.' + str(i_ind) + \
                                        '.' + str(j_ind) + '.npy')
            #learning
            center = learn_unlearn(len(common['on_off_map']),
                                   params['s']['n_iter'], [], [],
                                   params['rec_filter'], params['s'],
                                   common)
           
    	    tmp_cent = zeros((params['s']['rf_size'][0],
                             params['s']['rf_size'][1], len(center), 2))
            for cent_ind in range(len(center)):
                tmp_cent[:, :, cent_ind, :] = center[cent_ind]['patch']
           
            
            common['center'][i_ind, j_ind, :, :, :, :] = tmp_cent
            
            #record responce
            for block in range(ceil(float(len(common['on_off_map'])) / \
                               params['c']['block'])):
                name = 's1.' + str(block) + '.' + params['s']['type'] + '.npy'
           
                if i_ind == 0 and j_ind == 0:
                    s1_map = zeros((params['c']['rf_size'][0],
                                   params['c']['rf_size'][1], params['s']['n'],
                                   params['c']['block']))
                else:
                    s1_map = load(paths['map'] + name)

                range_var = range(block * params['c']['block'],
                                  minimum(len(common['on_off_map']),
                                  (block + 1) * params['c']['block']))
                s1_map[i_ind, j_ind, :, 0:len(range_var)] = \
                    record_response(range_var, False,
                                    common['center'][i_ind, j_ind, :, :, :, :],
                                    common)
                save(paths['map'] + name, s1_map)

            common['n_map'] = len(common['on_off_map'])
            #clear s1
            del common['on_off_map'] 

    cfile = open(paths['map'] + 'common.' + params['comp_name'] + '.' + \
                 params['s']['type'] + '.pck', 'w')
    pickle.dump(common, cfile)
    cfile.close()


if __name__ == "__main__":
    #params['test_mode'] = True
    layer_s1()

