import matplotlib.pyplot as plt
import scipy.signal as signal
import pickle
from config import *

def analyze():
    cfile = open(paths['map'] + 'common.c1.' + params['comp_name'] + '.' + params['s']['type'] + '.pck', 'r')
    common = pickle.load(cfile)
    #print common
    #----Draw for each C


    n_c = 16
    n_r = 16
    margin = .2
    selected_bound = .9

    width = n_c*(1 + margin) + margin
    height = n_r * (1 + margin) + margin

    x_vec = common['weight'][:, :, :, 3]
    y_vec = (x_vec - x_vec.min()) / (x_vec.max() - x_vec.min())

    count = 1

    #plt.figure(1)
    for i_ind in range(4):
        for j_ind in range(4):
            for f_ind in range(4):
                if common['weight'][i_ind, j_ind, f_ind, 3] > selected_bound:
                    rec = rf_reconstruction(reshape(array(common['center'])[i_ind,j_ind,:,:,f_ind,:], (params['s']['rf_size'][0], params['s']['rf_size'][1], 2)), array(params['rec_filter']), False)
                    im = plt.imshow(rec, interpolation='bilinear')
                    plt.show()
                    aa





def rf_reconstruction(weight, convo, polarized):
    n_r, n_c, n_ori = shape(weight)

    reconstruction = 0
    if polarized:
        for i_ori in range(n_ori / 2):
            tmp = signal.convolve2d(weight[:, :, i_ori], convo[:,:,i_ori], 'same')
            reconstruction += tmp

        for i_ori in range(n_ori/2, n_ori):
            tmp = signal.convolve2d(weight[:, :, i_ori], -convo[:,:,i_ori - n_ori/2], 'same')
            reconstruction += tmp
    else:
        for i_ori in range(n_ori):
            tmp = signal.convolve2d(weight[:, :, i_ori], convo[:,:,i_ori], 'same')
            reconstruction += tmp


    return -reconstruction


if __name__ == "__main__":
    analyze()
