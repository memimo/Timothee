
if ~exist('COMMON','var') || isempty(COMMON)
    load([PATH.map 'COMMON.' PARAM.compName '.' PARAM.s.type '.mat'])
end

run param % over-ride the PARAM stored in COMMON

rand('state',PARAM.randState);

% for sigma2 = 0% 2.^[ -4 ]
% for dec = 2.^[-14 -15 -16]%PARAM.s.thrDecay;
% for r = [4/3 1.5 5/3]
%     PARAM.c.thrDecay = dec;
%     PARAM.c.correctRatio = r;
    
% originalCompName = PARAM.compName;
% compName = [compName  '.dec.2-' int2str(-log(dec)/log(2)) '.r' num2str(r) ]

% disp(PARAM.compName)

% COMMON.s1Map = exp((COMMON.s1Map-1)/2/sigma2);

%************************
% COMPLEX CELL LEARNING *
%************************


[COMMON.weight COMMON.evol COMMON.thr COMMON.nFiringInput COMMON.nFiringOutput COMMON.nAboveThr]= learnInvariance(@getS1Map,[],PARAM.c);
save([PATH.map 'common_timothee.mat'], 'COMMON')
%analyze

% save([PATH.map 'COMMON.' PARAM.compName '.' PARAM.s.type '.' PARAM.c.type '.mat'],'COMMON','PARAM')

% currentPath = pwd;
% cd(PATH.map)
% copyfile(['COMMON.' PARAM.compName '.' PARAM.s.type '.' PARAM.c.type '.mat'],['../COMMON.mat'])
% cd(currentPath); clear currentPath;


% clear evol
% compName = originalCompName; % retrieve original compName

% end
% end
% end

