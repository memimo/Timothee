clear all

'sd'
run param
run module_conf
run setPath
'dsdsds'
% addpath('../learnS1')
% addpath('..')

rand('state',PARAM.randState);

% global PARAM
global COMMON

timedLogLn(['LEARN S1-C1 - ' PARAM.compName])

if PARAM.computeOnOff
    layerOnOff;
end
if PARAM.computeS1
    layerS1;
end
if PARAM.computeC1Timothee
    layerC1('timothee');
end
if PARAM.computeC1_foldiak
    layerC1('timothee');
end
if PARAM.computeC1_einhauser
    layerC1('einhauser');
end
if PARAM.computeC1_spraling
    layerC1('spratling');
end


quit()
