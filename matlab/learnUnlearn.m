function [center] = learnUnlearn(nMap,nIter,getMap,center,patch,recFilter,PARAM)
% Timothee Masquelier timothee.masquelier@alum.mit.edu Nov 2006
% nMap: number of input maps
% getMap: handle to a function i -> map
% center: eventually previous centers (to go on a computation), otherwise
% leave empty: []
% patch: eventually a patch to begin with (instead of random weights), otherwise leave empty: []
% recFilter: filter to use for reconstructions
% PARAM structure:
%   PARAM.n: number of cells
%   PARAM.RFSize: RF size of cells
%   PARAM.kWTA: limit the number of cells that can learn with a given image

rand('twister',5489);
%inhib
% PARAM.inhibWeight = 0.5;% horizontal inhibition
% PARAM.alpha = .01;
% PARAM.beta = .001;
% PARAM.gamma = .01;
% PARAM.lambda = 10;% sigmoid parameter
% PARAM.p = 1/PARAM.n;
PARAM.learningPeriod = 200;
PARAM.normalize = false; % not necessary to normalize weight (Tim 04/2007)
% PARAM.nu = 100;

tmp = feval(getMap,1);
nFilter = size(tmp{1},3);


if isempty(center)
    hasInitialPatch = ~isempty(patch);
    %instantiate main structure
    for iCenter=1:PARAM.n
        if ~hasInitialPatch
            patch{iCenter} = getRandomCenters(1,PARAM.RFSize,nFilter);
            center(iCenter).originalPatch = patch{iCenter};
            patch{iCenter} = patch{iCenter}(:); % for faster computations            
        else
            error('TO DO')
            center(iCenter).patch = reshape(patch(:,iCenter),[PARAM.RFSize nFilter]);
            center(iCenter).patch = center(iCenter).patch ./ norm(center(iCenter).patch(:));
        end
%         center(iCenter).minDistance = Inf;
        center(iCenter).activity = -1;
        %         center(iCenter).minBand = -1;
        %         center(iCenter).minI = -1;
        %         center(iCenter).minJ = -1;
        %         center(iCenter).distanceMap = {};
        center(iCenter).nFiring = 0;
        %         center(iCenter).evolDist = zeros(1,nIter);
        center(iCenter).evolThetaOriginal = zeros(1,nIter);
        center(iCenter).evolThreshold = zeros(1,nIter);
        %         center(iCenter).evolDeltaTheta = zeros(1,nIter);
        %         center(iCenter).evolTheta = zeros(1,nIter);
        %         center(iCenter).evolBand = zeros(1,nIter);
                center(iCenter).activityThreshold = 0;
        center(iCenter).learningRate = PARAM.learningRateMax;
    end
end

if PARAM.nu < Inf
    trace = 1/PARAM.nu*ones(1,PARAM.n);
else % mode with no trace
    trace = ones(1,PARAM.n);
end    

% for now
% threshold = zeros(1,PARAM.n);

        neutral = ones([PARAM.RFSize nFilter]);
        neutral = neutral(:)/ norm(neutral(:));


imOrder = randperm(nMap);
currentImIdx = 1;

idx = 1:PARAM.n;

disp(['learnUnlearn: size ' int2str(PARAM.RFSize) ' ...'])

tic

for iter=1:nIter
    
%     if mod(iter-1,100000)==0
%         %         fprintf(1,'.',iter);
%         for iCenter=1:PARAM.n
%             center(iCenter).patch = reshape(patch{iCenter}, [PARAM.RFSize nFilter]);
%         end
%         rec(center,recFilter,['iter ' int2str(iter-1)])
%     end

    map = getMap(imOrder(currentImIdx));
    %     map = feval(getMap,imOrder(currentImIdx));

    for iCenter=1:PARAM.n
        im = map{1};
        im = im(:);
        im = im / norm(im); % this normalization is very useful (Tim 04/2007)
        
        center(iCenter).rawActivity = dotProduct(im,patch{iCenter});
        center(iCenter).activity = center(iCenter).rawActivity/trace(iCenter);
%         if ~PARAM.normalize % note: if PARAM.normalize then this operation is unuseful, otherwise it helps balacing the competition between neurons
%                 center(iCenter).activity = center(iCenter).activity/norm(patch{iCenter});
%         end
    end
    
    % inhibition
%     apply inhibition to raw inputs
%     y = solveNonSelectiveMaxInhibCircuit([center.activity],PARAM.inhib);
    y = solveNonSelectiveInhibCircuit([center.activity],PARAM.inhib);
    for iCenter=1:PARAM.n % affect inhibited activites
        center(iCenter).activity = y(iCenter);
    end
    

    [maxActivity mostActive] = max([center.activity]);
    
    maxRawActivity = center(mostActive).rawActivity;

    if maxActivity >= center(mostActive).activityThreshold
%     if maxRawActivity >= center(mostActive).activityThreshold
        [center(mostActive) patch{mostActive}] = fire(center(mostActive), patch{mostActive}, im,iter,PARAM);
        
%         % reconstruction
%         center(mostActive).patch = reshape(patch{iCenter}, [PARAM.RFSize nFilter]);
%         rec(center(mostActive),recFilter,['iter ' int2str(iter-1)])
%     
%         if center.nFiring==50
%             error('stop')
%         end
        
        if PARAM.pushAway
            % push away other neurons
            for iCenter=1:PARAM.n
                if iCenter==mostActive
                    continue
                end
                if center(iCenter).activity >= center(iCenter).activityThreshold % would have fired later in the abscence of 1WTA
                    %                 newPatch = patch{iCenter} - center(iCenter).learningRate * 1 * ( im/norm(im(:)) - patch{iCenter} );
                    patch{iCenter} = patch{iCenter} + 1 * center(iCenter).learningRate * 1 * ( neutral - patch{iCenter} );
                    if PARAM.normalize
                        patch{iCenter} = patch{iCenter} / norm(patch{iCenter});
                    end
                end
            end
        end
    end
    
    % update traces
%     trace = [center.activity] / PARAM.nu + (1-1/PARAM.nu) * trace;
    trace = [center.rawActivity] / PARAM.nu + (1-1/PARAM.nu) * trace; % seems to be a better idea (Tim 04/2007)
    
%     for iCenter=1:PARAM.n % update thr for neuron that could have fired
%         if center(iCenter).activity>center(iCenter).activityThreshold
% %         if center(iCenter).rawActivity>center(iCenter).activityThreshold
%             center(iCenter).activityThreshold = PARAM.thrCoef*center(iCenter).activity;
% %             center(iCenter).activityThreshold = PARAM.thrCoef*center(iCenter).rawActivity;
%         end
%     end
        


    % update thr
        for iCenter=1:PARAM.n
            center(iCenter).activityThreshold = center(iCenter).activityThreshold * (1-PARAM.thrDecay);
        end


    %     fprintf(1,'.');
    if currentImIdx == nMap % last image
        currentImIdx = 1;
        imOrder = randperm(nMap); % reshuffle
    else
        currentImIdx = currentImIdx+1; % move on to next image
    end
    %             fprintf(1,'\n');



        if mod(iter,round(nIter/100))==0
            fprintf(1,'.');
        end
    if iter==600 && sum([center.nFiring])==0
        disp('Neurons do not fire.')
        break;
    end

end % iteration

fprintf(1,'\n');
toc

for iCenter=1:PARAM.n
    center(iCenter).patch = reshape(patch{iCenter}, [PARAM.RFSize nFilter]);
end

function [center newCenterPatch] = fire(center,patch,bestPatch,iter,PARAM)
% RFSize = size(bestPatch,1);
% alpha = (PARAM.size/4)^2; % see Mutch CVPR 2006
alpha = 1;

% inc firing counts
center.nFiring = center.nFiring+1;

% Hebbian modification: dW = a.y.(X-W)
y = 1;
% y = (1-center.minDistance/PARAM.firingThreshold/PARAM.size)^1;
% y = exp( -center.minDistance / 2 / center.sigma^2 / alpha );
% y = max(0,1/center.d2Threshold*(center.d2Threshold-center.minDistance));
% y = y^2 * 3;
% newCenterPatch = patch + PARAM.learningRate * y * bestPatch;
newCenterPatch = patch + center.learningRate * y * ( bestPatch - patch );
% newCenterPatch = patch + center.learningRate * y * ( bestPatch - patch );

% normalize
if PARAM.normalize
    newCenterPatch = newCenterPatch / norm(newCenterPatch);
end

% Reporting
% % compute delta theta
% center.evolDeltaTheta(iter) = acos(min(1,dotProduct(patch(:),newCenterPatch(:))));
% % compute Theta
% center.evolTheta(iter) = acos(min(1,dotProduct(patch(:),normedBestPatch(:))));
% % dist
% center.evolDist(iter) = center.maxActivity;
% % band
% center.evolBand(iter) = center.minBand;
% % dist
if PARAM.normalize
    center.evolThetaOriginal(iter) = acos(min(1,dotProduct(newCenterPatch,center.originalPatch)));
else
    center.evolThetaOriginal(iter) = acos(min(1,dotProduct(newCenterPatch/norm(newCenterPatch),center.originalPatch)));
end

% % store new center
% patch = newCenterPatch;
% 
% compute new threshold if necessary
if center.nFiring <= PARAM.learningPeriod && mod(center.nFiring,10)==0
    % geometric
    center.learningRate = center.learningRate * (PARAM.learningRateMin/PARAM.learningRateMax)^(1/(PARAM.learningPeriod/10));
    %     center.d2Threshold = center.d2Threshold * (PARAM.d2ThresholdMin/PARAM.d2ThresholdMax)^(1/(PARAM.learningPeriod/10));
    % arithmetic
    %     center.learningRate = center.learningRate - (PARAM.learningRateMax-PARAM.learningRateMin)/(PARAM.learningPeriod/10);
    %     center.d2Threshold = center.d2Threshold - (PARAM.d2ThresholdMax-PARAM.d2ThresholdMin)/(PARAM.learningPeriod/10);

    %     center.activityThreshold = exp(-center.d2Threshold / (2*PARAM.sigma^2*alpha));
end
% center.activityThreshold = center.maxActivity * center.activityThreshold ;

center.activityThreshold = PARAM.thrCoef*center.activity;
% center.activityThreshold = PARAM.thrCoef*center.rawActivity;

% center.activityThreshold = center.activity;
center.evolThreshold(iter) = center.activityThreshold;

% function center = computeMaxActivity(center,sigma,alpha)
%
% center.minDistance = Inf;
% for iBand=1:length(center.distanceMap)
%     [distance idx] = min(center.distanceMap{iBand}(:));
%     if distance < center.minDistance % found a better match
%         center.minDistance = max(0,distance);% just to avoid tiny negative values
%         center.maxActivity = exp( -distance / 2 / (sigma^2) / alpha );
%         center.minBand = iBand;
%         center.minJ = floor((idx-1)/size(center.distanceMap{iBand},1))+1;
%         center.minI = idx - (center.minJ-1) * size(center.distanceMap{iBand},1);
%     end
% end

function centers = getRandomCenters(n,RFSize,nFilter)
MEAN = 20;
STD = 1;

centers = [];

for iCenter=1:n
    %     center = MEAN + STD * randn([PARAM.size PARAM.size nFilter]);
    %         center = rand([PARAM.size PARAM.size nFilter]);
    %     center = ones([PARAM.size PARAM.size nFilter]);

    center = 0*1+1*rand([RFSize nFilter]);

    %     sigma = RFSize(1)/4;
    %     x = abs( [1:RFSize(1)] - (RFSize(1)+1)/2 )' * ones(1,RFSize(2));
    %     y = abs( [1:RFSize(2)] - (RFSize(2)+1)/2 )' * ones(1,RFSize(1));
    %     y = y';
    %     g = exp(-( x.^2 + y.^2 )/2/(sigma^2));
    %
    %     for n=1:nFilter
    %         center(:,:,n) = g;
    %     end
    %
    center = max(0,center); % rectify negative values
    centers(:,:,:,iCenter) = center ./ norm(center(:));
end

