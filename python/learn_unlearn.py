from config import *
import numpy
from numpy import *
from tools import *

def learn_unlearn(n_map, n_iter, center, patch, in_param, common):
    '''
    Timothee Masquelier timothee.masquelier@alum.mit.edu Nov 2006
    n_map: number of input maps
    center: eventually previous centers (to go on a computation), otherwise
    leave empty: []
    patch: eventually a patch to begin with (instead of random weights),
           otherwise leave empty: []
    rec_filter: filter to use for reconstructions
    PARAM structure:
    PARAM.n: number of cells
    PARAM.RFSize: RF size of cells
    PARAM.kWTA: limit the number of cells that can learn with a given image
    '''
    
    #Set random seed for test mode
    if params['test_mode']:
        random.mtrand.seed(5489)

    params['s']['normalize'] = False
    params['s']['learning_period'] = 200.0
    tmp = get_onoff_map(0, common)

    n_filter = shape(tmp)[2]

    if center == []:
        has_initial_path = (patch != [])
        #instantiate main structure
        for cent_ind in range(in_param['n']):
            if not has_initial_path:
                patch.append(get_random_centers(1, in_param['rf_size'],
                             n_filter))
                center.append({})
                center[cent_ind]['original_patch'] = patch[cent_ind]
                patch[cent_ind] = patch[cent_ind].flatten()
            else:
                raise NameError('Error')


            center[cent_ind]['activity'] = numpy.float64(-1.0)
            center[cent_ind]['n_firing'] = 0.0
            center[cent_ind]['evol_theta_original'] = zeros(n_iter,
                                                            numpy.float64)
            center[cent_ind]['evol_threshold'] = zeros(n_iter, numpy.float64)
            center[cent_ind]['activity_threshold'] = numpy.float64(0.0)
            center[cent_ind]['learning_rate'] = in_param['learning_rate_max']


    trace = None
    if in_param['nu'] < inf:
        trace = 1.0 / in_param['nu'] * ones(in_param['n'], numpy.float64)
    else: #mode with no trace
        trace = ones(in_param['n'], numpy.float64)


    neutral = ones((in_param['rf_size'][0], in_param['rf_size'][1], n_filter),
                   numpy.float64)
    neutral = neutral.flatten(1) / linalg.norm(neutral.flatten(1))

    im_order = permute(n_map) 
    current_im_idx = 0

    print 'lenarnUnlearn: size ', str(in_param['rf_size']), '...'

    
    #Iteration-------------------------------------------------
    for iteration in range(n_iter):
        #if iteration == 38115:
        #    pdb.set_trace()

        onoff_map = get_onoff_map(im_order[current_im_idx], common)
        for cent_ind in range(in_param['n']):
            img = onoff_map.flatten(1)
            #this normalization is very useful (Tim 04/2007)
            img = img / linalg.norm(img)

            center[cent_ind]['raw_activity'] = dot(img, patch[cent_ind])
            center[cent_ind]['activity'] = center[cent_ind]['raw_activity'] / \
                trace[cent_ind]
            
        #inhibition
        #apply inhibition to raw inputs
        for cent_ind in range(in_param['n']): #affect inhibited activites
            center[cent_ind]['activity'] = \
            solve_nonselect_inhib_circuit(center[cent_ind]['activity'],
                                              in_param['inhib'])
            
         
        max_activity = float64(0.0)
        most_active = 0
        for cent_ind in range(in_param['n']):
            if center[cent_ind]['activity'] > max_activity:
                max_activity = center[cent_ind]['activity']
                most_active = cent_ind

       
        #Check for fire---       
        if max_activity >= center[most_active]['activity_threshold']:
            
            center[most_active], patch[most_active] = fire(center[most_active],
                                                           patch[most_active],
                                                           img, iteration,
                                                           in_param)
            
            if in_param['push_away']:
                for cent_ind in range(in_param['n']):
                    if cent_ind == most_active:
                        continue
                    if center[cent_ind]['activity'] >= \
                        center[cent_ind]['activity_threshold']:
                        patch[cent_ind] = patch[cent_ind] + 1.0 * \
                        center[cent_ind]['learnng_rate'] * 1.0 * \
                        (neutral - patch[cent_ind])
                        if in_param['normalize']:
                            patch[cent_ind] = patch[cent_ind] / \
                            linalg.norm(patch[cent_ind])


        #Update traces
        for cent_ind in range(in_param['n']):
            trace[cent_ind] = center[cent_ind]['raw_activity'] / in_param['nu']\
            + float64(1.0 - 1.0 / in_param['nu']) * trace[cent_ind]

        #Update threshold
        for cent_ind in range(in_param['n']):
            center[cent_ind]['activity_threshold'] = \
            center[cent_ind]['activity_threshold'] * \
            float64(1.0 - in_param['thr_decay'])

        if current_im_idx == n_map -1: #last image
            current_im_idx = 0
            im_order = permute(n_map) #repermute
        else:
            current_im_idx += 1 #move on to next image

        if mod(iteration, round(n_iter / 100.0)) == 0:
            print '.',

        if iteration == 599:
            tsum = 0            
            for cent_ind in range(in_param['n']):
                tsum += center[cent_ind]['n_firing']
            if tsum == 0:
                print 'Neurons do not fire'
                break

        

    #-End of oiteration ---------------------------


    for cent_ind in range(in_param['n']):
        center[cent_ind]['patch'] = reshape(patch[cent_ind],
                                            (in_param['rf_size'][0],
                                            in_param['rf_size'][1], n_filter),
                                            order='F')

    return center



#---Fire Function--------------------------------------------------
def fire(center, patch, best_patch, iteration, in_param):
    '''Fire the neuron'''

    #inc firing_rate
    center['n_firing'] += 1

    #Hebbian modification: dW = a.y.(X-W)
    y_var = float64(1.0)
    
    new_center_patch = patch + (center['learning_rate'] * y_var * \
                                (best_patch - patch)).astype(numpy.float64)

    #Normalize
    if in_param['normalize']:
        new_center_patch = new_center_patch / linalg.norm(new_center_patch)
        center['evol_theta_original'][iteration] = \
        arccos(minimum(1, dot(new_center_patch, center['original_patch'])))
    else:
        center['evol_theta_original'][iteration] = \
        arccos(minimum(1, dot(new_center_patch / linalg.norm(new_center_patch),
               center['original_patch'].flatten())))

    if (center['n_firing'] <= in_param['learning_period']) and \
        mod(center['n_firing'], 10) == 0:
        center['learning_rate'] = float64(center['learning_rate'] * \
                                          (in_param['learning_rate_min'] / \
                                          in_param['learning_rate_max']) ** \
                                          (1.0 / (in_param['learning_period'] /\
                                          10.0)))

    center['activity_threshold']  = in_param['thr_coef'] * center['activity']
    center['evol_threshold'][iteration] = center['activity_threshold']

    return center, new_center_patch

#-Get Random Centers Function ----------------------------------------
def get_random_centers(n_size, rf_size, n_filter):
    '''return a random center matrix'''

    centers = zeros((rf_size[0], rf_size[1], n_filter, n_size), dtype=float64)

    for cent_ind in range(n_size):
        center = float64(random.random((rf_size[0], rf_size[1], n_filter)))
       
        center = maximum(0, center) 
        centers[:, :, :, cent_ind] = center / linalg.norm(center.flatten(1))

    return centers


#- Solve non selective inhib circuit function ----------
def solve_nonselect_inhib_circuit(input_data, alpha):
    '''Solve non selective inhibit circuits'''

    if alpha == 0:
        return input_data
    i_var = alpha * mean(maximum(0, input_data-0)) - 0
    y_var = maximum(0, input_data - i_var)

    return y_var



'''if __name__ == "__main__":
    common = {}
    common['on_off_map'] = load(paths['map'] + 'onoff.' + str(0) + '.' \
                                + str(0) + '.npy')
    center = learn_unlearn(len(common['on_off_map']), params['s']['n_iter'], [],
                           [], params['rec_filter'], params['s'], common)'''

