'''Compute the last layer or layer C1 based on four learning algorithems:
    Timothee
    Foldiak
    Spratling
    Einhauser'''

from analyze import *
from config import *
import pickle


def learn_c1(learning_type, common):
    '''Layer C1 computation'''

    if common == {}:
        cfile = open(paths['map'] + 'common.' + params['comp_name'] + '.' + \
                     params['s']['type'] + '.pck', 'r')
        common = pickle.load(cfile)
    
    #timothee
    if learning_type == 'timothee':
        import learn_invariance
        common = learn_invariance.learn_invariance([], params['c'],
                                                       common)

        cfile = open(paths['map'] + 'common.c_timothee.' + params['comp_name'] \
                     + '.' + params['s']['type'] + '.pck', 'w')
        pickle.dump(common, cfile)
        cfile.close()

    #Foldiak
    elif learning_type == 'foldiak':
        import foldiak
        common = foldiak.foldiak1991([], params['c'], common)
        
        cfile = open(paths['map'] + 'common.c_foldiak.' + params['comp_name'] \
                     + '.' + params['s']['type'] + '.pck', 'w')
        pickle.dump(common, cfile)
        cfile.close()

    #Spartling
    elif learning_type == 'spratling':
        import spratling
        common = spratling.spratling2005(common['n_map'],
                                             params['s']['n'], params['c']['n'],
                                             params['c']['rf_size'], [],
                                             params['c'], common)
        
        cfile = open(paths['map'] + 'common.c_spratling.' + \
                     params['comp_name'] + '.' + params['s']['type'] + \
                     '.pck', 'w')
        pickle.dump(common, cfile)
        cfile.close()

    #Einhauser
    elif learning_type == 'einhasuer':
        import einhauser
        common = einhauser.einhauser(common['n_map'], params['s']['n'],
                                         params['c']['n'],
                                         params['c']['rf_size'],
                                         [], params['c'], common)

        cfile = open(paths['map'] + 'common.c_einhauser.' + \
                     params['comp_name'] + '.' + params['s']['type'] + '.pck',
                     'w')
        pickle.dump(common, cfile)
        cfile.close()



if __name__ == "__main__":
    #learn_c1('timothee', common)
   learn_c1('foldiak', {})
    #learn_c1('einhasuer', common)
    #learn_c1('spratling', common)
