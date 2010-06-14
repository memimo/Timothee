function LayerC1(method)

if ~exist('COMMON','var') || isempty(COMMON)
    load([PATH.map 'COMMON.' PARAM.compName '.' PARAM.s.type '.mat'])
end

run ../conf/param % over-ride the PARAM stored in COMMON

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

if method == 'timothee'
    [COMMON.weight COMMON.evol COMMON.thr COMMON.nFiringInput COMMON.nFiringOutput COMMON.nAboveThr]= learnInvariance(@getS1Map,[],PARAM.c);
elseif method == 'foldiak'
    [COMMON.weight COMMON.evol]= foldiak1991(@getS1Map,[],PARAM.c);
elseif method == 'einhauser'
    [COMMON.weight COMMON.evol COMMON.thr COMMON.nFiringInput COMMON.nFiringOutput]= einhauser2002(@getS1Map,COMMON.nMap,PARAM.s.n,PARAM.c.n,[PARAM.c.RFSize],[]);
elseif method = 'spartling'
    [COMMON.weight COMMON.evol COMMON.nFiringInput COMMON.nFiringOutput]= spratling2005(@getS1Map,COMMON.nMap,PARAM.s.n,PARAM.c.n,[PARAM.c.RFSize],[]);
end

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

