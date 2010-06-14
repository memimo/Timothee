function [center, nCouldHaveFired, nFired, w ] = antiHebbianLearning(nMap,nIter,getMap,center,patch,w,recFilter,PARAM)
% Timothee Masquelier timothee.masquelier@alum.mit.edu Nov 2006
% see Foldiak 1990
% nMap: number of input maps
% getMap: handle to a function i -> map
% center: eventually previous centers (to go on a computation), otherwise leave empty: []
% patch: eventually a patch to begin with (instead of random weights), otherwise leave empty: []
% recFilter: filter to use for reconstructions
% PARAM structure:
%   PARAM.n: number of cells
%   PARAM.RFSize: RF size of cells
%   PARAM.kWTA: limit the number of cells that can learn with a given image


%inhib
% PARAM.inhibWeight = 0.5;% horizontal inhibition
% PARAM.alpha = .01;
PARAM.beta = .001;
% PARAM.gamma = .01;
PARAM.lambda = 10;% sigmoid parameter
% PARAM.p = 1/PARAM.n;

tmp = feval(getMap,1);
nFilter = size(tmp{1},3);


if isempty(center)
    %instantiate main structure
    for iCenter=1:PARAM.n
        if isempty(patch)
            center(iCenter).patch = getRandomCenters(1,PARAM.RFSize,nFilter);
        else
            center(iCenter).patch = reshape(patch(:,iCenter),[PARAM.RFSize nFilter]);
            center(iCenter).patch = center(iCenter).patch ./ norm(center(iCenter).patch(:));
        end
        center(iCenter).originalPatch = center(iCenter).patch;
        center(iCenter).minDistance = Inf;
        center(iCenter).maxActivity = -1;
        center(iCenter).minBand = -1;
        center(iCenter).minI = -1;
        center(iCenter).minJ = -1;
        center(iCenter).distanceMap = {};
        center(iCenter).nFiring = 0;
        %         center(iCenter).evolDist = zeros(1,nIter);
        center(iCenter).evolThetaOriginal = zeros(1,nIter);
        center(iCenter).evolThreshold = zeros(1,nIter);
        %         center(iCenter).evolDeltaTheta = zeros(1,nIter);
        %         center(iCenter).evolTheta = zeros(1,nIter);
        %         center(iCenter).evolBand = zeros(1,nIter);
%         center(iCenter).activityThreshold = .5;
        center(iCenter).learningRate = PARAM.learningRateMax;
    end
end

% for now
% threshold = zeros(1,PARAM.n);
threshold = 0*ones(1,PARAM.n);

if PARAM.nu < Inf
    trace = 1/PARAM.nu*ones(1,PARAM.n);
else % mode with no trace
    trace = ones(1,PARAM.n);
end    

if isempty(w)
%     w = zeros(PARAM.n,PARAM.n);
    w = -.0 * 1/PARAM.n * (ones(PARAM.n,PARAM.n)-eye(PARAM.n));
%     disp(w)
end

nCouldHaveFired(nIter)=0;
nFired(nIter)=0;

imOrder = randperm(nMap);
currentImIdx = 1;

idx = 1:PARAM.n;

disp(['antiHebbianLearning: size ' int2str(PARAM.RFSize) ' ...'])

tic

for iter=1:nIter
    map = getMap(imOrder(currentImIdx));
    %     map = feval(getMap,imOrder(currentImIdx));
    
    for iCenter=1:PARAM.n
        im = map{1};
        im = im / norm(im(:));
        %             center(iCenter).activity = dotProduct(im(:),center(iCenter).patch(:))/norm(im(:));
%         input(iCenter) = dotProduct(im(:),center(iCenter).patch(:));
        input(iCenter) = dotProduct(im(:),center(iCenter).patch(:)) / trace(iCenter);
    end

%     y = solveInhibitoryNetwork([input-threshold]',w,zeros(PARAM.n,1),PARAM.lambda);
%     y = solveInhibitoryNetwork([input]',w,zeros(PARAM.n,1),PARAM.lambda);
    y = input;
    
%     yBin = y>.5;
    yBin = y>=threshold;
    
    nCouldHaveFired(iter) = sum(yBin);
    nFired = min(PARAM.kWTA, nCouldHaveFired); % not very useful...
    
    if PARAM.kWTA<PARAM.n % kWTA competition. we need to sort activities
        [y idx] = sort(-[y]); y = -y;
        yBin(PARAM.kWTA+1:PARAM.n) = 0; % do not let fire the last ones in the list
    end

    
    for i=1:PARAM.n
        if yBin(i)% firing criterion
            % hebbian learning
            if iter>0%500
                center(idx(i)).patch = center(idx(i)).patch + PARAM.beta*(im-center(idx(i)).patch);
                center(idx(i)).patch = center(idx(i)).patch / norm(center(idx(i)).patch(:)); % normalize
                center(idx(i)).evolThetaOriginal(iter) = acos(min(1,dotProduct(center(idx(i)).patch(:)/norm(center(idx(i)).patch(:)),center(idx(i)).originalPatch(:))));
                center(idx(i)).nFiring = center(idx(i)).nFiring+1;
                threshold(idx(i)) = PARAM.thrCoef * y(idx(i));
            end
        else % firing criterion
            if PARAM.kWTA<PARAM.n % activities are sorted, meaning we don't need to examine farther
                break
            end
        end
        center(idx(i)).evolThreshold(iter) = threshold(idx(i)); % th storage
    end

%     % anti hebbian learning
%     if iter > 0%500
%         w = w - PARAM.alpha*(double(yBin')*yBin-PARAM.p^2);
%         for i=1:PARAM.n
%             w(i,i)=0;
%         end
% %         w = w .* (1-eye(PARAM.n));
%         w = min(0,w);
%     end
    
%     % threshold modif & storage
%     threshold = threshold + PARAM.gamma*(yBin-PARAM.p);
    threshold = threshold * (1-PARAM.thrDecay);

    % update traces
    trace = y / PARAM.nu + (1-1/PARAM.nu) * trace;


    %     fprintf(1,'.');
    if currentImIdx == nMap % last image
        currentImIdx = 1;
        imOrder = randperm(nMap); % reshuffle
    else
        currentImIdx = currentImIdx+1; % move on to next image
    end
    %             fprintf(1,'\n');



    if mod(iter,10000)==0
        fprintf(1,'.',iter);
        %         reconstruction(center,filters,fSiz,nFilter,['iter ' int2str(iter)])
    end
    if iter==600 && sum([center.nFiring])==0
        disp('Neurons do not fire. Increase distance threshold and/or increase sigma')
        break;
    end

%     for iCenter=1:PARAM.n
%         center(iCenter).activityThreshold = center(iCenter).activityThreshold * (1-PARAM.thrDecay);
%     end
end % iteration

fprintf(1,'\n');
toc

% function center = fire(center,bestPatch,iter,PARAM)
% % RFSize = size(bestPatch,1);
% % alpha = (PARAM.size/4)^2; % see Mutch CVPR 2006
% alpha = 1;
% 
% % inc firing counts
% center.nFiring = center.nFiring+1;
% 
% normedBestPatch =  bestPatch/norm(bestPatch(:));
% 
% % Hebbian modification: dW = a.y.(X-W)
% y = 1;
% % y = (1-center.minDistance/PARAM.firingThreshold/PARAM.size)^1;
% % y = exp( -center.minDistance / 2 / center.sigma^2 / alpha );
% % y = max(0,1/center.d2Threshold*(center.d2Threshold-center.minDistance));
% % y = y^2 * 3;
% % newCenterPatch = center.patch + PARAM.learningRate * y * bestPatch;
% newCenterPatch = center.patch + center.learningRate * y * ( normedBestPatch - center.patch );
% % newCenterPatch = center.patch + center.learningRate * y * ( bestPatch - center.patch );
% 
% % normalize
% newCenterPatch = newCenterPatch / norm(newCenterPatch(:));
% 
% % Reporting
% % % compute delta theta
% % center.evolDeltaTheta(iter) = acos(min(1,dotProduct(center.patch(:),newCenterPatch(:))));
% % % compute Theta
% % center.evolTheta(iter) = acos(min(1,dotProduct(center.patch(:),normedBestPatch(:))));
% % % dist
% % center.evolDist(iter) = center.maxActivity;
% % % band
% % center.evolBand(iter) = center.minBand;
% % % dist
% center.evolThetaOriginal(iter) = acos(min(1,dotProduct(newCenterPatch(:),center.originalPatch(:))));
% 
% % store new center
% center.patch = newCenterPatch;
% 
% % compute new threshold if necessary
% if center.nFiring <= PARAM.learningPeriod && mod(center.nFiring,10)==0
%     % geometric
%     center.learningRate = center.learningRate * (PARAM.learningRateMin/PARAM.learningRateMax)^(1/(PARAM.learningPeriod/10));
%     %     center.d2Threshold = center.d2Threshold * (PARAM.d2ThresholdMin/PARAM.d2ThresholdMax)^(1/(PARAM.learningPeriod/10));
%     % arithmetic
%     %     center.learningRate = center.learningRate - (PARAM.learningRateMax-PARAM.learningRateMin)/(PARAM.learningPeriod/10);
%     %     center.d2Threshold = center.d2Threshold - (PARAM.d2ThresholdMax-PARAM.d2ThresholdMin)/(PARAM.learningPeriod/10);
% 
%     %     center.activityThreshold = exp(-center.d2Threshold / (2*PARAM.sigma^2*alpha));
% end
% % center.activityThreshold = center.maxActivity * center.activityThreshold ;
% center.activityThreshold = center.maxActivity;
% 
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

    center = 1+.05*rand([RFSize nFilter]);

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

