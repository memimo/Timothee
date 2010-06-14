import unittest
import scipy.io as sio
from config import *
import layer_lgn
import layer_s1
import layer_c1

class timothee_model_test(unittest.TestCase):
    '''  All test are run with saved data from MAtlab implenetation on same data '''
    
    def setUp(self):
        
        #Parameter initializations
        self.test_data_path = 'test_mat_data2/'
        params['test_mode'] = True
        paths['image'] = paths['home'] + 'data/'
        


    #--- LGN Layer --------------------------
    def test_lgn(self):
        for i_ind in range(params['c']['rf_size'][0]):
            for j_ind in range(params['c']['rf_size'][1]):
                on_off_map  = load(paths['map'] + 'onoff.' + str(i_ind) + '.' + str(j_ind) + '.npy')

                m_file = sio.loadmat(self.test_data_path  + 'onOff.' + str(i_ind + 1) + '.' + str(j_ind + 1) + '.mat')
                map_matlab = m_file['onOffMap']

                for f_ind in range(len(on_off_map)):
                    self.assertTrue(allclose(map_matlab[0,f_ind], on_off_map[f_ind]))





    '''def test_layer_s1(self):
            

        s1_map = load(paths['map'] + 's1.0.nS16_inhibMean0.0.npy')
        m_file = sio.loadmat(self.test_data_path  + 's1.01.nS16_inhibMean0.00.mat')
        map_matlab = m_file['s1Map']

        self.assertTrue(allclose(s1_map, map_matlab))'''


if __name__ == "__main__":
    unittest.main()
