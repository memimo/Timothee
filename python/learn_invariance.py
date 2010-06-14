from numpy import *
from tools import *


def learn_invariance(weight, param, common):
    '''Timothee learning algorithm'''
    
    #Set random seed for test run
    if params['test_mode']:
        random.mtrand.seed(5489)

    #when true, reinforce when previous input and current output correlate
    #(eg Einhauser 2002)
    correl_prev_input_curr_output = False 
    geom_factor = (param['learning_rate_f'] \
                   / param['learning_rate_i']) ** \
        (1.0 / floor(param['n_iter'] / 1000.0))
    param['a_pos'] = param['learning_rate_i']
    param['a_neg'] = -param['a_pos'] / (param['n_feat'] * \
                                        prod(param['rf_size'])) * \
        param['correct_ratio']
   
    if weight == []:
        weight = zeros((param['rf_size'][0], param['rf_size'][1],
                       param['n_feat'], param['n']), float64)
        for n_ind in range(param['n']):
            weight[:, :, :, n_ind] = get_initial_weight(param['rf_size'],
                                                        param['n_feat'])
            
        #init middle layer
        if param['input_nu'] < inf:
            tr_input = 1.0 / param['input_nu'] * ones((param['rf_size'][0], 
                                                      param['rf_size'][1],
                                                      param['n_feat']), float64)
        else:
            tr_input = ones((param['rf_size'][0], param['rf_size'][1],
                            param['n_feat']), float64)

        if param['output_nu'] < inf:
            tr_output = 1.0 / param['output_nu'] * ones(param['n'], float64)
        else:
            tr_output = ones(param['n'], float64)

        input_threshold = zeros((param['rf_size'][0], param['rf_size'][1],
                                param['n_feat']), float64)
        output_threshold = ones(param['n'], float64)

        previous_winner_pair = -1.0 * ones((4), float64)

        # reporting
        n_firing_output = zeros((param['n']), float64)
        n_firing_input = zeros((param['rf_size'][0], param['rf_size'][1],
                               param['n_feat']), float64)
        n_above_thr = zeros((param['n_iter']), float64) 

        #init evol
        evol = zeros((param['n'], floor(param['n_iter'] / 100)), float64)

    print 'Learning invariances...'
    offset = floor(random.random() * param['n_iter'])

    
    for i_ind in range(param['n_iter']):
        
        effective = False
        #current input
        gross_input, common = get_s1_map(i_ind + offset, common)
    
        #update traces
        tr_input = gross_input / param['input_nu'] + \
        (1.0 - 1.0 / param['input_nu']) * tr_input
        #normalization by traces
        input_var = gross_input / tr_input
        n_above_thr[i_ind]  = (input_var > input_threshold).sum()
        
        #get max input
        input_max, idx = multi_dimensional_max(input_var)
        input_winner_i = idx[0]
        input_winner_j = idx[1]
        input_winner_f = idx[2]
        gross_input_max = gross_input[idx[0], idx[1], idx[2]]

        #current output
        gross_output = zeros(param['n'], float64)
        for n_ind in range(param['n']):
            if param['use_soft_max']:
                #softmax
                w = weight[:, :, :, n_ind]
                gross_output[n_ind] = sof_max(w, input_var, param['p'],
                                              param['q'], param['r'])
            else:
                gross_output[n_ind] = (weight[:, :, :, n_ind] * input_var).max()


        #update traces
        if i_ind == 0:
            tr_output = gross_output
        else:
            tr_output = gross_output / param['output_nu'] + \
            (1.0 - 1.0 / param['output_nu']) * tr_output

        #normalization by traces
        output = zeros(param['n'], float64)
        for n_ind in range(param['n']):
            output[n_ind] = gross_output[n_ind] / tr_output[n_ind]

        output_max = tr_output.max()
        output_winner = tr_output.argmax()

        #correlation between previous input and current output
        #(like Einhauser 2002)
        if correl_prev_input_curr_output: 
            #Learning
            if i_ind > 0 and previous_input_max > \
            input_threshold[previous_input_winner_i, previous_input_winner_j,
            previous_input_winner_f]:
                #update threshold
                input_threshold[previous_input_winner_i,
                previous_input_winner_j, previous_input_winner_f] = \
                param['thr_coef'] * previpus_input_max
                #update wiehgts: stpre winner-winner weight
                win_weight = weight[previous_input_winner_i,
                previous_input_winner_j, previous_intput_f, output_winner]
                #depress everyone
                weight[:, :, :, output_winner] = \
                weight[:, :, :, output_winner] + param['a_neg'] * \
                weight[:, :, :, output_winner] * (1-weight[:, :, :,
                                                  output_winner])

                #reproting
                n_firing_input[previous_input_winner_i, previous_winner_j,
                previous_input_winner_f] = \
                n_firing_input[previous_input_winner_i, previous_input_winner_k,
                previous_winner_f] + 1
                n_firing_output[output_winner] = \
                n_firing_output[output_winner] + 1
        else: #correlation between previous output and current input
            #learning
            if (i_ind > 0) and (input_max > input_threshold[input_winner_i,
                                input_winner_j, input_winner_f]):
                #update the
                input_threshold[input_winner_i, input_winner_j,
                input_winner_f] = param['thr_coef'] * input_max
                #store winner weight
                win_weight = weight[input_winner_i, input_winner_j,
                input_winner_f, previous_output_winner]
                #depress everyone
                weight[:, :, :, previous_output_winner] = \
                weight[:, :, :, previous_output_winner] + param['a_neg'] * \
                weight[:, :, :, previous_output_winner] * \
                (1.0-weight[:, :, :, previous_output_winner])
                #potentiate synapse between winners
                weight[input_winner_i, input_winner_j, input_winner_f,
                previous_output_winner] = win_weight + param['a_pos'] * \
                win_weight * (1.0 - win_weight)
                # reporting
                n_firing_input[input_winner_i, input_winner_j,
                input_winner_f] = n_firing_input[input_winner_i, input_winner_j,
                input_winner_f] + 1
                n_firing_output[previous_output_winner] = \
                n_firing_output[previous_output_winner] + 1


        weight = maximum(weight, 0)
        weight = minimum(weight, 1)

        
        if (mod(i_ind + 1, 100) == 0):
            evol[:, int((i_ind + 1) / 100) - 1] = weight.sum()

        #threshold decay
        input_threshold = input_threshold * (1.0 - param['thr_decay'])
        #previous := current
        previous_input_max = input_max
        previous_input_winner_i = input_winner_i
        previpus_input_winner_j = input_winner_j
        previous_input_winner_f = input_winner_f
        previous_output_winner = output_winner
        previous_output_max = output_max
        previous_gross_input_max = gross_input_max

        if mod(i_ind + 1, 1000) == 0:
            param['a_pos'] = param['a_pos'] * geom_factor
            param['a_neg'] = param['a_neg'] * geom_factor
        
        if mod(i_ind, round(param['n_iter'] / 100)) == 0:
            print '.'

    print '\ntr ' + str(tr_output)
    print 'thr ' + str(output_threshold)
    print 'out' + str(output)

    common['weight'] = weight 
    common['evol'] = evol 
    common['thr'] = input_threshold 
    common['n_firing_input'] = n_firing_input
    common['n_firing_output'] = n_firing_output 
    common['n_above_thr'] = n_above_thr
    return common


def get_initial_weight(patch_size, n_feat):

    weight = 0.75 * (1.0 + 0.0 * (random.random((patch_size[0], patch_size[1], n_feat))- 0.5))
    return weight 
