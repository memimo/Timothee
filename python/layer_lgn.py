from numpy import *
from config import *
import scipy as sp
import scipy.signal as signal
import glob

def compute_layer():
    '''Compute LGN ON-and OFF-center unit layer    '''

    #init
    internal_pos = [0, 0]
    
    for i_ind in range(params['c']['rf_size'][0]):
        internal_pos[0] =  (i_ind + 1 - floor((params['c']['rf_size'][0] +
                            1)/2.0)) * params['s']['sub_sampling'] * \
                            params['on_off']['sub_sampling']
        for j_ind in range(params['c']['rf_size'][1]):
            internal_pos[1] = (j_ind + 1 - floor((params['c']['rf_size'][1] +
                              1)/2.0)) * params['s']['sub_sampling'] * \
                              params['on_off']['sub_sampling']
            on_off_map = compute_map(params['dog_filter'], paths['image'],
                         params['crop_pos'], internal_pos,
                         params['s']['rf_size'],
                         params['on_off']['sub_sampling'], params['zoom'])
            save(paths['map'] + 'onoff.' + str(i_ind) + '.' + str(j_ind) + \
                 '.npy', on_off_map )
            del on_off_map



def compute_map(in_filter, img_path, crop_pos, internal_pos, map_size,
                sub_sampling, zoom):
    '''Compute the map for each unit'''
    
    #initialize variables
    filter_size = shape(in_filter)[0]
    map_var = zeros((map_size[0], map_size[1], 2), float64)
    crop_size = (array(sub_sampling) * (array(map_size) -1) + \
                (array(filter_size) * ones((1,2), float64)))[0]

    image_files = glob.glob(img_path + '*.tif')
    image_files.sort()
    num_img = len(image_files)

    map_list = [None]*(num_img*(shape(crop_pos)[0]))

    for i_ind in range(num_img):
        img = sp.misc.imread(image_files[i_ind])
        if zoom < 1:
            img = shrink(img, zoom)

        for j_ind in range(shape(crop_pos)[0]):
            im_cropped = img[crop_pos[j_ind][0] + internal_pos[0] - 1 - \
                         floor(0.5 * (crop_size[0] -1)):crop_pos[j_ind][0] + \
                         internal_pos[0] + ceil(0.5*(crop_size[0] - 1)), \
                         crop_pos[j_ind][1] + internal_pos[1] - 1 - \
                         floor(0.5*(crop_size[1] - 1)): crop_pos[j_ind][1] + \
                         internal_pos[1]  + ceil(0.5*(crop_size[1] - 1))]
            
            
            im_cropped = im_cropped.astype(float64) / 255.0
            map_var = signal.convolve2d(im_cropped, in_filter, 'valid')
            map_var = map_var[0::sub_sampling, 0::sub_sampling]
     
            map_list[j_ind*(num_img) + i_ind ] = map_var

        if mod(i_ind, 100000) == 0:
            print '.',
    
    return map_list




#if __name__ == "__main__":
    #compute_layer()
