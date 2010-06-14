''' The config file'''

from numpy import *


#functions
def __filter_dog(filter_size, sigma1, sigma2):
    ''' Diffrence-of-Guassian (DoG) filter'''

    x_grid = tile(arange(1, filter_size + 1, dtype=float64), (filter_size, 1))
    y_grid = x_grid.T
    d2 = float64((x_grid - filter_size / 2.0 - 0.5) ** 2.0 + \
                 (y_grid - filter_size / 2.0 - 0.5) ** 2.0)

    filter_var = float64(1.0 / sqrt(2.0 * pi) * (1.0 / sigma1 * \
                         exp(-d2 / 2.0 / (sigma1 ** 2.0)) - 1.0 / sigma2 * \
                         exp(-d2 / 2.0 / (sigma2 ** 2.0))))

    filter_var = filter_var - mean(filter_var.flatten(1))
    filter_var = filter_var / linalg.norm(filter_var.flatten(1))

    return filter_var



#Network parameters
params = {}

#Stimuli
params['zoom'] = 1
params['crop_pos'] = []
for ii_ind in range(1, 10):
    for jj_ind in range(1, 12):
        params['crop_pos'].append([ii_ind * 25, jj_ind * 25])

#On-Off layer
params['on_off'] = {}
params['on_off']['rf_size'] = 7
params['on_off']['sub_sampling'] = 1

#Real DOG
sigma2 = float64(params['on_off']['rf_size'] / 5.0)
sigma1 = float64(sigma2 / 1.6)
params['dog_filter'] = __filter_dog(params['on_off']['rf_size'], sigma1, sigma2)

params['rec_filter'] = array([params['dog_filter'], -params['dog_filter']])

#Simple cell layer
params['s'] = {}
params['s']['n'] = 16
params['s']['rf_size'] = [7, 7]
params['s']['sub_sampling'] = 3
params['s']['n_iter'] = int(round(0.5 * 99 * 17009))
params['s']['kwta'] = 1
params['s']['learning_rate_max'] = float64(0.1)
params['s']['learning_rate_min'] = float64(0.01)
params['s']['thr_decay'] = float64(2.0) ** -15.0
params['s']['thr_coef'] = float64(1.0)
params['s']['push_away'] = False
params['s']['nu'] = float64(100.0)
params['s']['inhib'] = 0.0
params['s']['type'] = 'nS' + str(params['s']['n']) + '_inhibMean' + \
str(params['s']['inhib'])



#Complex cell layer
params['c'] = {}
params['c']['n'] = 4
params['c']['rf_size'] = [4, 4]
params['c']['n_iter'] = int(round(0.5 * 3 * 99 * 17009))
params['c']['learning_rate_i'] = 2.0 ** -3.0
params['c']['learning_rate_f'] = 2.0 ** -1.0
params['c']['block'] = round(50000.0 * 4 * 4 / prod(params['c']['rf_size']) * \
                             16.0 / params['s']['n'])
params['c']['correct_ratio'] = 1.5
params['c']['n_feat'] = params['s']['n']
#input
params['c']['thr_decay'] = 2.0 ** -15.0
params['c']['thr_coef'] = 1.0
params['c']['input_nu'] = inf
params['c']['tr_input_effective'] = False
#output
params['c']['output_nu'] = 1.5
params['c']['tr_output_effective'] = False
params['c']['type'] = 'nC' + str(params['c']['n'])
params['c']['use_soft_max'] = False

params['comp_name'] = 'ref_DoG' + str(params['on_off']['rf_size']) + '-' + \
str(params['on_off']['rf_size'] / sigma2) + '_S' + \
str(params['s']['rf_size'][0]) + '_C' + str(params['c']['rf_size'][0]) + \
'_shift' + str(params['s']['sub_sampling'])


#Test mode
params['test_mode'] = True

#File paths----------
paths = {}

if params['test_mode']:
    #paths['image'] = '../test/data/'
    paths['image'] = '/u2/mmirza/timothy/python/data/'
    paths['map'] = '../test/map/'
    paths['mat'] = '../test/mat/'
    paths['mat-conf'] = '../matlab/script/'
else:
    paths['home'] = '/u2/mmirza/timothy/python/'
    #paths['image'] = paths['home'] + 'data/'
    paths['image'] = '/u2/mmirza/timothy/LearnS1C1/movie01-06/'
    paths['map'] = paths['home'] + 'map/'

paths['matlab-sh'] = '/software/matlab2009a/bin/matlab'



