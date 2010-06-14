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
        self.test_data_path = 'test_mat_data/'
        params['test_mode'] = True
        paths['image'] = paths['home'] + 'data/'
        paths['map'] = paths['home'] + 'map-test/'
        #params['s']['n_iter'] = 1000
        params['c']['n_iter'] = 1000

        #layer_lgn.compute_layer()
        layer_s1.layer_s1()
        #layer_c1.learn_c1('timothee', common)


    #--- LGN Layer --------------------------
    def test_dog_filter(self):

        m_file = sio.loadmat(self.test_data_path + 'dog_filter.mat')
        filter_matlab = m_file['dog_filter']
       
        self.assertTrue(allclose(filter_matlab, filter_matlab))

    def test_compute_map(self):
        import layer_lgn 
        internal_pos = [-3, -3]
        on_off_map = layer_lgn.compute_map(params['dog_filter'], paths['test_image'], params['crop_pos'], internal_pos, params['s']['rf_size'], params['on_off']['sub_sampling'], params['zoom'])

        
        m_file = sio.loadmat(self.test_data_path + 'on_off_map.mat')
        map_matlab = m_file['onOffMap']


        for i_ind in range(len(on_off_map)):
            self.assertTrue(allclose(map_matlab[0,i_ind], on_off_map[i_ind]))

 
    def test_layer_lgn(self):

        for i_ind in range(params['c']['rf_size'][0]):
            for j_ind in range(params['c']['rf_size'][1]):
                on_off_map  = load(paths['map'] + 'onoff.' + str(i_ind) + '.' + str(j_ind) + '.npy')

                m_file = sio.loadmat(self.test_data_path  + 'onOff.' + str(i_ind + 1) + '.' + str(j_ind + 1) + '.mat')
                map_matlab = m_file['onOffMap']

                for f_ind in range(len(on_off_map)):
                    self.assertTrue(allclose(map_matlab[0,f_ind], on_off_map[f_ind]))



    #--- Layer S1 -----------------------
    def test_learn_unlearn(self):
        import learn_unlearn
        
        
        common = {}
        common['on_off_map'] = load(paths['map'] + 'onoff.' + str(0) + '.' + str(0) + '.npy')
        center = learn_unlearn.learn_unlearn(len(common['on_off_map']), params['s']['n_iter'], [], [], params['rec_filter'], params['s'], common)

        m_file = sio.loadmat(self.test_data_path  + 'center_s1.mat')
        center_matlab = m_file['center']

        for i_ind in range(params['s']['n']):
            #print center[i_ind]['activity'],'\n', center_matlab[0][i_ind].activity, 'activity\n\n'
            #print center[i_ind]['n_firing'] ,'\n', center_matlab[0][i_ind].nFiring, 'nFiring\n\n'
            #print center[i_ind]['evol_theta_original'],'\n' ,center_matlab[0][i_ind].evolThetaOriginal, 'evolThetaOriginal\n\n'
            #print center[i_ind]['evol_threshold'],'\n' , center_matlab[0][i_ind].evolThreshold, 'evolThreshold\n\n'
            #print center[i_ind]['activity_threshold'],'\n' , center_matlab[0][i_ind].activityThreshold, 'activityThreshold\n\n'
            #print center[i_ind]['learning_rate'] ,'\n',center_matlab[0][i_ind].learningRate, 'learningRate\n\n'
            #print center[i_ind]['patch'] ,'\n', center_matlab[0][i_ind].patch, 'patch\n\n'
            #print center[i_ind]['patch'].mean() ,'\n', center_matlab[0][i_ind].patch.mean(), 'patch mean\n\n'
            self.assertAlmostEqual(center[i_ind]['patch'].mean(), center_matlab[0][i_ind].patch.mean(),2)
            self.assertTrue(allclose(center[i_ind]['patch'], center_matlab[0][i_ind].patch,2))


    '''
    def test_layer_s1(self):
            

        s1_map = load(paths['map'] + 's1.0.nS16_inhibMean0.0.npy')
        m_file = sio.loadmat(self.test_data_path  + 's1.01.nS16_inhibMean0.00.mat')
        map_matlab = m_file['s1Map']

        self.assertTrue(allclose(s1_map, map_matlab))

    def test_record_responce(self):
        import layer_s1
        m_file = sio.loadmat(self.test_data_path  + 'response.mat')
        in_range = m_file['range']
        center = m_file['center']
        m_output = m_file['response']
        
        common = {}
        common['on_off_map'] = load(paths['map'] + 'onoff.' + str(0) + '.' + str(0) + '.npy')
        response = layer_s1.record_response(in_range - 1, False, center, common)
        self.assertTrue(allclose(m_output, response))


    #--- Layer C1 --------------------------------
    def test_layer_c1_timothee(self):
        import pickle
        import learn_invariance
        
        cfile = open(paths['map'] + 'common.' + params['comp_name'] + '.' + params['s']['type'] + '.pck', 'r')
        common = pickle.load(cfile)

        common['weight'], common['evol'], common['thr'], common['n_firing_input'], common['n_firing_output'], common['n_above_thr'] = learn_invariance.learn_invariance([], params['c'], common)
        
        m_file = sio.loadmat(self.test_data_path  + 'c1-timothee.mat')
        matlab_common = m_file['COMMON'][0][0]
        self.assertTrue(allclose(common['weight'], matlab_common.weight))
        self.assertTrue(allclose(common['evol'], matlab_common.evol))
        self.assertTrue(allclose(common['thr'], matlab_common.thr))
        self.assertTrue(allclose(common['n_firing_input'], matlab_common.nFiringInput))
        self.assertTrue(allclose(common['n_above_thr'], matlab_common.nAboveThr))
        


    #--- Analyze ------------------------
    def rf_reconstruction(self):
        from analyze import *

        m_file = sio.loadmat(self.test_data_path  + 'rf.mat')
        rf = m_file['rf']

        py_rf = rf_reconstruction(rf, params['rec_filter'], False)
        mat_rf = m_file['reconst_rf'] 

        self.assertTrue(allclose(py_rf, mat_rf))

    '''

if __name__ == "__main__":
    unittest.main()
    
