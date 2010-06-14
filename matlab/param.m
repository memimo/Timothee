global PARAM

PARAM.computeOnOff = false;
PARAM.computeS1 = false;
PARAM.computeC1_timothee = false;
PARAM.computeC1_foldiak = false;
PARAM.computeC1_einhauser = false;
PARAM.computeC1_spratling = false;

PARAM.randState = 0;

% Stimuli
PARAM.zoom = 1;
% PARAM.cropPos = [120 160]
% note: images are 320x240
% PARAM.cropPos = round( PARAM.zoom * ...
%     [ [50 50]  ; [100 50]  ; [150 50]  ; [200 50] ; ...
%       [50 100] ; [100 100] ; [150 100] ; [200 100] ; ...
%       [50 150] ; [100 150] ; [150 150] ; [200 150] ; ...
%       [50 200] ; [100 200] ; [150 200] ; [200 200] ; ...
%       [50 250] ; [100 250] ; [150 250] ; [200 250] ] );
PARAM.cropPos = [];
for i=1:9
    for j=1:11
        PARAM.cropPos = [PARAM.cropPos ; [i*25 j*25] ];
    end
end
  
% PARAM.cropPos = round( PARAM.zoom * ...
%                         [ [50 50] ; [150 50] ; [50 150] ; [150 150] ; [50 250] ; [150 250] ] );


% ON-OFF layer
PARAM.onOff.RFSize = 7;
PARAM.onOff.subSampling = 1; % no sub-sampling among ON-OFF maps for now

% Real DoG
% sigma2 = (PARAM.onOff.RFSize+1)/2/3; % so that DoG is approximately 0 just outside the filter
sigma2 = (PARAM.onOff.RFSize)/5; % /4 or 5 seems a good idea. 6 is too much
sigma1 = sigma2/1.6; % biological ratio
PARAM.DoGFilter = DoG(PARAM.onOff.RFSize,sigma1,sigma2);

% % Einhauser 2002
% LPF = [[1 2 1]; [2 4 2]; [1 2 1]];
% LAP = [[0 -1 0]; [-1 4 -1]; [0 -1 0]];
% PARAM.DoGFilter = conv2(LPF,LAP,'full');
% PARAM.DoGFilter = PARAM.DoGFilter / norm(PARAM.DoGFilter(:));
% PARAM.DoGFilter = LAP / norm(LAP(:));

PARAM.recFilter(:,:,:,1)=PARAM.DoGFilter;
PARAM.recFilter(:,:,:,2)=-PARAM.DoGFilter;

% Simple cell layer
PARAM.s.n = 16;
PARAM.s.RFSize = [7 7];
PARAM.s.subSampling = 3;
PARAM.s.nIter = round(0.5*99*17009);
PARAM.s.kWTA = 1; % only one cell can learn at a given location
% PARAM.s.kWTA = PARAM.s.n; % no kWTA
PARAM.s.learningRateMax =  0.1;
PARAM.s.learningRateMin =  0.01;
% PARAM.s.d2ThresholdMin = .50; % not used anymore
% PARAM.s.d2ThresholdMax = .50; % not used anymore
PARAM.s.thrDecay = 2^-15;%2^-11;
PARAM.s.thrCoef = 1;%.96;
PARAM.s.pushAway = false;
PARAM.s.nu = 100; % smooth temporal average. put Inf not to use
PARAM.s.inhib = 0.0;
PARAM.s.type = ['nS' int2str(PARAM.s.n) '_inhibMean' sprintf('%.2f',PARAM.s.inhib) ];
disp(['S type: ' PARAM.s.type])

% % anti hebbian learning
% PARAM.s.alpha = 2^-7;
% PARAM.s.beta = 2^-4;
% PARAM.s.gamma = 2^-5;
% PARAM.s.lambda = 10;% sigmoid parameter
% PARAM.s.p = 2^-7;%2^-9;


% Complex cell layer
PARAM.c.n = 4;
PARAM.c.RFSize = [4 4];
PARAM.c.nIter = round(.5*3*99*17009);
PARAM.c.learningRate_i = 2^-3; % avoid more than 2^-3
PARAM.c.learningRate_f = 2^-1; % avoid more than 2^-1
PARAM.c.block = round(50000 * 4*4/prod(PARAM.c.RFSize) * 16/PARAM.s.n);
PARAM.c.correctRatio = 1.5; % the greater the more unselected. should be around 4/3 or 1.5
PARAM.c.nFeat = PARAM.s.n;
% input
PARAM.c.thrDecay = 2^-15;% should be such that about 5% of stimuli are effective
PARAM.c.thrCoef = 1.0;%.95;
PARAM.c.inputNu = Inf;%100; % smooth temporal average. put Inf not to use
PARAM.c.trInputEffective = false; % only compute average with effective stim
% output
PARAM.c.outputNu = 1.5;%Inf; % smooth temporal average. put Inf not to use (Tim 12/2006: I recommend 1000)
PARAM.c.trOutputEffective = false; % only compute average with effective stim (Tim 12/2006: I recommend true)
PARAM.c.type = ['nC' int2str(PARAM.c.n)];
PARAM.c.useSoftMax = false;
disp(['C type: ' PARAM.c.type])

% computation name
PARAM.compName = [ 'ref_DoG' int2str(PARAM.onOff.RFSize) '-' num2str(PARAM.onOff.RFSize/sigma2) '_S' int2str(PARAM.s.RFSize(1)) '_C' int2str(PARAM.c.RFSize(1))  '_shift' int2str(PARAM.s.subSampling)];
