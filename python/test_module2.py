import unittest
import scipy.io as sio
import subprocess
import os
from config import *
import layer_lgn
import layer_s1
import layer_c1
import pickle
import learnS1C1


class timothee_model_test(unittest.TestCase):
    '''  All test are run with saved data from MAtlab implenetation on same data '''


    def test_init(self):
        ''' Run the matla to produce the templates as well as runnig the python codes to produce the python codes as welll for comparasion'''
  
         #Compute Python templates----------
        learnS1C1.run_all()

    
        #---Compute Matlab templates----------
        lgn_compute = True
        s1_compute = True
        c1_timothee_compute = True
        c1_foldiak_compute = True
        c1_einhauser_compute = True
        c1_spratling_compute = True
        
        #Check for computeense of mat files for lgn
        computes = zeros((params['c']['rf_size'][0], params['c']['rf_size'][1]), bool)
        for i_ind in range(params['c']['rf_size'][0]):
            for j_ind in range(params['c']['rf_size'][1]):
                if os.path.isfile(paths['mat'] + 'onOff.' + str(i_ind + 1) + '.' + str(j_ind + 1) + '.mat'):
                    computes[i_ind, j_ind] = True

 
        if computes.all():
            lgn_compute = False


        #Check computeense of mat files for layer S1---
        tmp_lgn = load(paths['map'] + 'onoff.' + str(0) + '.' + str(0) + '.npy')
        count = ceil(float(len(tmp_lgn)) / params['c']['block'])
        computes = zeros((count), bool)
        for block in range(count):
            if block < 9:
                name = paths['mat'] + 's1.0' + str(block + 1) + '.'  + params['s']['type'] + '0.mat'
            else:
                name = paths['mat'] + 's1.' + str(block + 1) + '.' + params['s']['type'] + '0.mat'
           
            if os.path.isfile(name):
                computes[block] = True        
        
        if computes.all():
            s1_compute = False
        

        #Check computeense of mat files for layer C1---
        #Timothee    
        if os.path.isfile(paths['mat'] + 'common_timothee.mat'):
            c1_timothee_compute = False
        #Foldiak
        if os.path.isfile(paths['mat'] + 'common_foldiak.mat'):
            c1_foldiak_compute = False
        #Einhasuer
        if os.path.isfile(paths['mat'] + 'common_einhauser.mat'):
            c1_einhauser_compute = False
        #Spratling
        if os.path.isfile(paths['mat'] + 'common_spratling.mat'):
            c1_spratling_compute = False

        #Run the matlab to compute necessary mat file templates
   
        fconf = open(paths['mat-conf'] + 'module_conf.m', 'w')
        fconf.write('PARAM.computeOnOff = ' + str(lgn_compute).lower() + ';\n')
        fconf.write('PARAM.computeS1 = ' + str(s1_compute).lower() + ';\n')
        fconf.write('PARAM.computeC1_timothee = ' + str(c1_timothee_compute).lower() + ';\n')
        fconf.write('PARAM.computeC1_foldiak = ' + str(c1_foldiak_compute).lower() + ';\n')
        fconf.write('PARAM.computeC1_einhauser = ' + str(c1_einhauser_compute).lower() + ';\n')
        fconf.write('PARAM.computeC1_spratling = ' + str(c1_spratling_compute).lower() + ';\n')

        fconf.close()

        
        subprocess.call([paths['matlab-sh'], '-nodesktop -nosplash -nojvm -r learnS1C1']) 
       

        self.assertTrue(True)

    
    '''def test_layer_lgn(self):
        Compare the lgn layer output

        for i_ind in range(params['c']['rf_size'][0]):
            for j_ind in range(params['c']['rf_size'][1]):
                on_off_map  = load(paths['map'] + 'onoff.' + str(i_ind) + '.' + str(j_ind) + '.npy')
                m_file = sio.loadmat(paths['mat']  + 'onOff.' + str(i_ind + 1) + '.' + str(j_ind + 1) + '.mat')
                map_matlab = m_file['onOffMap']
               
                for f_ind in range(len(on_off_map)):
                    self.assertTrue(allclose(map_matlab[0,f_ind], on_off_map[f_ind]))


    def test_layer_s1(self):
        Comapre the layer s1

        tmp_lgn = load(paths['map'] + 'onoff.0.0.npy')
        count = ceil(float(len(tmp_lgn)) / params['c']['block'])
        if count == 0: count = 1

        for i in range(count):
            s1map  = load(paths['map'] + 's1.' + str(i) + '.nS16_inhibMean0.0.npy')
            num = ''
            if i < 9:
                num = '0' + str(i + 1)
            else:
                num = str(i+ 1)
            m_file = sio.loadmat(paths['mat']  + 's1.' + num + '.nS16_inhibMean0.00.mat')
            map_matlab = m_file['s1Map']

            self.assertTrue(allclose(s1map, map_matlab))

       



    def test_layer_c1_timothee(self):
        Compare the layer c1 with timothee method

        cfile = open(paths['map'] + 'common.c_timothee.ref_DoG7-5.0_S7_C4_shift3.nS16_inhibMean0.0.pck', 'r')
        c1map = pickle.load(cfile)
        m_file = sio.loadmat(paths['mat']  + 'common_timothee.mat')
        matlab_common = m_file['COMMON'][0][0]

        self.assertTrue(allclose(c1map['weight'], matlab_common.weight))
        self.assertTrue(allclose(c1map['center'], matlab_common.center))
        

    '''
    def test_layer_c1_foldiak(self):
        '''Compare the layer c1 with foldiak method'''

        cfile = open(paths['map'] + 'common.c_foldiak.ref_DoG7-5.0_S7_C4_shift3.nS16_inhibMean0.0.pck', 'r')
        c1map = pickle.load(cfile)
        m_file = sio.loadmat(paths['mat']  + 'common_foldiak.mat')
        matlab_common = m_file['COMMON'][0][0]

        self.assertTrue(allclose(c1map['weight'], matlab_common.weight))
        self.assertTrue(allclose(c1map['center'], matlab_common.center))
       


    def test_layer_c1_einhauser(self):
        '''Compare the layer c1 with einhauser method'''

        cfile = open(paths['map'] + 'common.c_einhauser.ref_DoG7-5.0_S7_C4_shift3.nS16_inhibMean0.0.pck', 'r')
        c1map = pickle.load(cfile)
        m_file = sio.loadmat(paths['mat']  + 'common_einhauser.mat')
        matlab_common = m_file['COMMON'][0][0]

        self.assertTrue(allclose(c1map['weight'], matlab_common.weight))
        self.assertTrue(allclose(c1map['center'], matlab_common.center))
        


    def test_layer_c1_spratling(self):
        '''Compare the layer c1 with spartling method'''

        cfile = open(paths['map'] + 'common.c_spratling.ref_DoG7-5.0_S7_C4_shift3.nS16_inhibMean0.0.pck', 'r')
        c1map = pickle.load(cfile)
        m_file = sio.loadmat(paths['mat']  + 'common_spratling.mat')
        matlab_common = m_file['COMMON'][0][0]

        self.assertTrue(allclose(c1map['weight'], matlab_common.weight))
        self.assertTrue(allclose(c1map['center'], matlab_common.center))
        


if __name__ == "__main__":
    unittest.main()
    
