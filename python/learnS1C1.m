clear all

run ../matlab/script/param
run ../matlab/script/module_conf
run ../matlab/script/setPath

% addpath('../learnS1')
% addpath('..')

rand('state',PARAM.randState);

% global PARAM
global COMMON



if PARAM.computeOnOff
    run ../matlab/script/layerOnOff
end
if PARAM.computeS1
    run ../matlab/script/layerS1
end
if PARAM.computeC1_timothee
    run ../matlab/script/layerC1Timothee
end
run ../matlab/script/module_conf
if PARAM.computeC1_foldiak
    run ../matlab/script/layerC1Foldiak
end
run ../matlab/script/module_conf
if PARAM.computeC1_einhauser
    run ../matlab/script/layerC1Einhauser
end
run ../matlab/script/module_conf
if PARAM.computeC1_spratling
    run ../matlab/script/layerC1Spartling
end


quit()
