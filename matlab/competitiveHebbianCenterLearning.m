function [center, nCouldHaveFired, nFired ] = competitiveHebbianCenterLearning(nMap,nIter,getMap,weightSharingMode,center,patch,recFilter,PARAM,cropPos)
% Timothee Masquelier timothee.masquelier@alum.mit.edu August 2006
% nMap: number of input maps
% getMap: handle to a function i -> map
% center: eventually previous centers (to go on a computation), otherwise
% leave empty: []
% patch: eventually a patch to begin with (instead of random weights),
% otherwise leave empty: []
% recFilter: filter to use for reconstructions
% PARAM structure:
%   PARAM.n: number of cells
%   PARAM.RFSize: RF size of cells
%   PARAM.kWTA: limit the number of cells that can learn with a given image
%   PARAM.d2ThresholdMax: (should be around .45) threshold at the begining of computation in terms of the square of the max distance between center and current input for the cell to fire
%   PARAM.d2ThresholdMin: (should be around .35) idem at the end of computation
%   PARAM.learningRateMax: (should be around .2) learningRate at the begining of computation
%   PARAM.learningRateMin: (should be around .02) learningRate at the end of computation

% for now
PARAM.oneByOne = false;
PARAM.learningPeriod = 200;
% PARAM.sigma = 2^-2; % not used anymore
PARAM.inhibCoef = 0;

%inhib
PARAM.inhibWeight =0;% horizontal inhibition
PARAM.lambda = 10;% sigmoid parameter


tmp = feval(getMap,1);
nFilter = size(tmp{1},3);

% alpha = (PARAM.RFSize/4)^2; % see Mutch CVPR 2006
alpha = 1;

% offset between image and convolution map
% offset = floor((PARAM.RFSize-1)/2); % use this offset with 'same' convo type
offset = 0; % use this offset with 'valid' convo type

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
        %         center(iCenter).sigma = PARAM.sigma;
        center(iCenter).maxActivity = -1;
        center(iCenter).minBand = -1;
        center(iCenter).minI = -1;
        center(iCenter).minJ = -1;
        center(iCenter).distanceMap = {};
        center(iCenter).nFiring = 0;
%         center(iCenter).evolDist = zeros(1,nIter);
        center(iCenter).evolThetaOriginal = zeros(1,nIter);
%         center(iCenter).evolDeltaTheta = zeros(1,nIter);
%         center(iCenter).evolTheta = zeros(1,nIter);
%         center(iCenter).evolBand = zeros(1,nIter);
        center(iCenter).d2Threshold = PARAM.d2ThresholdMax; %max tolerance at the beginning
%         center(iCenter).activityThreshold = exp(-center(iCenter).d2Threshold / (2*PARAM.sigma^2*alpha));
        center(iCenter).activityThreshold = 0.7;
        center(iCenter).learningRate = PARAM.learningRateMax;
    end
end

nCouldHaveFired(nIter)=0;
nFired(nIter)=0;

imOrder = randperm(nMap);
currentImIdx = 1;

if PARAM.oneByOne
    nActiveCenter = 1;
else
    nActiveCenter = PARAM.n;
end

% % debug
% dirList = dir(['..\movie01-06\*.tif']);
% nIm = length(dirList);
% f=figure;


disp(['competitiveHebbianCenterLearning: size ' int2str(PARAM.RFSize) ' ...'])
tic
for iter=1:nIter
%     map = getMap(imOrder(currentImIdx));
    map = feval(getMap,imOrder(currentImIdx));

    %     % tmp : useful for patch learning
    %     mapSize = round(PARAM.mapSize*PARAM.RFSize);
    %     [nR nC tmp] = size(map{1});
    %     randOffset = [ ceil(rand*(nR-mapSize)) ceil(rand*(nC-mapSize)) ];
    %     map{1} = map{1}(randOffset(1):randOffset(1)+mapSize-1,randOffset(2):randOffset(2)+mapSize-1,:);

    for iCenter=1:nActiveCenter
        if weightSharingMode
            flippedCenter = center(iCenter).patch(end:-1:1,end:-1:1,:); % patch is flipped because WindowedPatchDistance uses conv2 (and not filer2)
            for iBand=1:length(map)
                center(iCenter).distanceMap{iBand} = WindowedPatchDistance(map{iBand},flippedCenter,'valid');
            end
            center(iCenter) = computeMaxActivity(center(iCenter),PARAM.sigma,alpha);
        else
            im = map{1};
            normIm = norm(im(:));
            center(iCenter).minDistance = 1-dotProduct(im(:),center(iCenter).patch(:))/normIm;
            center(iCenter).maxActivity = 1-center(iCenter).minDistance;
%             center(iCenter).maxActivity = exp( -center(iCenter).minDistance / 2 / (PARAM.sigma^2) / alpha ) / center(iCenter).activityThreshold ;
        end
    end
    
    if PARAM.inhibWeight>0 % horizontal inhibition (see Foldiak 1990) before competion
        y = solveInhibitoryNetwork([[center.maxActivity]-[center.activityThreshold]]',PARAM.inhibWeight/nActiveCenter*ones(nActiveCenter,nActiveCenter),zeros(nActiveCenter,1),PARAM.lambda);
        [activity idx] = sort(-[y]); activity = -activity;
        firingThreshold = .5*ones(1,nActiveCenter);
    else % competition only
        [activity idx] = sort(-[center.maxActivity]); activity = -activity;
        firingThreshold = [center.activityThreshold];
    end

%         i = 1;
    %     disp(activity)
        nCouldHaveFired(iter) = sum(activity >= firingThreshold);
    nFiredThisTime = 0;
        for i=1:PARAM.n
%     while activity(i) >= 0.5% firing criterion
                if activity(i) >= firingThreshold(idx(i))
        %         disp([dist ;idx])
        %         disp(['Firing cell: ' int2str(idx(i))])

        % retrieve firing patch (indices were painful to get right, but they are now!)
        if weightSharingMode
            bestPatch = zeros([PARAM.RFSize nFilter]);
            [nR nC tmp] = size(map{center(idx(i)).minBand});
            bestPatch(  max(1,offset-center(idx(i)).minI+2):min(PARAM.RFSize(1),nR-center(idx(i)).minI+offset+1), ...
                max(1,offset-center(idx(i)).minJ+2):min(PARAM.RFSize(2),nC-center(idx(i)).minJ+offset+1), : ) ...
                = map{center(idx(i)).minBand}(  max(1,center(idx(i)).minI-offset):min(nR,center(idx(i)).minI-offset+PARAM.RFSize(1)-1), ...
                max(1,center(idx(i)).minJ-offset):min(nC,center(idx(i)).minJ-offset+PARAM.RFSize(2)-1), : );
        else
            bestPatch = im;
        end

        %         % debug
        %         d1 = 1 - dotProduct(bestPatch(:),center(idx(i)).patch(:)) / norm(bestPatch(:));
        %         d2 = center(idx(i)).minDistance;
        %         disp([d1-d2])
        %         if      offset-center(idx(i)).minI+2 > 1 ...
        %             ||  nR-center(idx(i)).minI+offset+1 < PARAM.RFSize ...
        %             ||  offset-center(idx(i)).minJ+2 > 1 ...
        %             ||  nC-center(idx(i)).minJ+offset+1 < PARAM.RFSize
        %             fprintf(1,'o',iter); % Best patch is partly out
        %         else
        %             fprintf(1,'i',iter); % Best patch is completely in
        %         end

        % fire the unit and trigger learning rule
        %         disp([int2str(idx(i)) ' fires:'])
        %         fprintf(1,'%d.',idx(i));
        center(idx(i)) = fire(center(idx(i)),bestPatch,iter,PARAM);
        
%         % debug
%         cropNumber = floor((imOrder(currentImIdx)-1)/nIm)+1;
%         imNumber = imOrder(currentImIdx) - (cropNumber-1)*nIm;
%         image = imread(['..\movie01-06\' dirList(imNumber).name]);
%         
% %         subplot(1,3,1);imagesc(image);colormap(gray);
% 
%         internalPos(1)=1;
%         internalPos(2)=1;        
%         cropSize = 1*(80-1) + (7*ones(1,2));
%         imCropped = image(  cropPos(cropNumber,1)+internalPos(1)-1-floor(.5*(cropSize(1)-1)):cropPos(cropNumber,1)+internalPos(1)-1+ceil(.5*(cropSize(1)-1)), ...
%                             cropPos(cropNumber,2)+internalPos(2)-1-floor(.5*(cropSize(2)-1)):cropPos(cropNumber,2)+internalPos(2)-1+ceil(.5*(cropSize(2)-1)));
%                         
%         subplot(1,3,1);imagesc(imCropped);colormap(gray);
%         title('image');
%         
%         rec = RFReconstruction(im, recFilter, false);
%         subplot(1,3,2);imagesc(rec);colormap(gray)
%         title('convo')
%         
%         rec = RFReconstruction(center(idx(i)).patch, recFilter, false);
%         subplot(1,3,3);imagesc(rec);colormap(gray)
%         title('RF')
% 
%         set(f,'Name',['nFiring=' int2str(center(idx(i)).nFiring) ' - image ' dirList(imNumber).name ]);
%         pause
%                         

        nFired(iter) = nFired(iter)+1;
        nFiredThisTime = nFiredThisTime+1;
        %         if nFired > nCouldHaveFired % debug
        %             error('nFired > nCouldHaveFired')
        %         end

        %         % inhib other (less active) units
        %         minDistHaveChanged = false;
        %         iRange = max(1,center(idx(i)).minI-(maxOverlap-1)):min(nR,center(idx(i)).minI+(maxOverlap-1));
        %         jRange = max(1,center(idx(i)).minJ-(maxOverlap-1)):min(nC,center(idx(i)).minJ+(maxOverlap-1));
        %         for j=i+1:PARAM.n
        %             % inhib other distance maps
        %             center(idx(j)).distanceMap{center(idx(i)).minBand}(iRange,jRange) = Inf;
        %             % recompute min distance if necessary
        %             if      center(idx(j)).minBand == center(idx(i)).minBand ... %only inhib same band maps
        %                 &&  abs(center(idx(j)).minI-center(idx(i)).minI) < maxOverlap ...
        %                 &&  abs(center(idx(j)).minJ-center(idx(i)).minJ) < maxOverlap
        %                     center(idx(j)) = computeMinDistance(center(idx(j)));
        %                     dist(j) = center(idx(j)).minDistance;
        %                     minDistHaveChanged = true;
        % %                     disp(['Impacted min for cell ' int2str(idx(j))])
        %             end
        %         end

        if nFiredThisTime==PARAM.kWTA % no need to inhibit in that case
            break;
        end

%         % new mode of inhibition
%         maxActivityMightHaveChanged(1:PARAM.n) = false;
%         for iBand = 1:length(map)
%             %             inhibitionMap = ( center(idx(i)).distanceMap{iBand} < 5*dist(i) );
%             %             inhibitionMap = ( center(idx(i)).distanceMap{iBand} < 1.25*PARAM.firingThreshold * PARAM.RFSize );
%             %             inhibitionMap = ( center(idx(i)).distanceMap{iBand} < 1.0 * (-2*center(idx(i)).sigma^2*alpha*log(PARAM.firingThreshold)) );
%             inhibitionMap = center(idx(i)).distanceMap{iBand} <= 1 - (1-PARAM.inhibCoef)*(1-center(idx(i)).d2Threshold);
%             %             inhibitionMap = 1/10 * .1 * ( 1 - center(idx(i)).distanceMap{iBand} );
%             % debug
%             %             if idx(i)==27
%             %                 disp(['Band ' int2str(center(idx(i)).minBand) ' fired: Band ' int2str(iBand) ' - inhib ' int2str(sum(inhibitionMap(:))) ' pix (' num2str(100*mean(inhibitionMap(:)>0)) '%)'])
%             %             end
%             %             if iBand == center(idx(i)).minBand
%             %                   disp(['Cell ' int2str(idx(i)) ' fired: Band ' int2str(iBand) ' - inhib ' int2str(sum(inhibitionMap(:))) ' pix (' num2str(100*mean(inhibitionMap(:)>0)) '%)'])
%             %             end
% 
%             for j=i+1:PARAM.n
%                 if idx(j) > nActiveCenter % don't need to inhibit unactive units
%                     continue
%                 end
%                 % inhib by setting distance to Inf
%                 center(idx(j)).distanceMap{iBand} = center(idx(j)).distanceMap{iBand} + 10 * inhibitionMap;
%                 % put a flag to recompute min distance if necessary
%                 if      center(idx(j)).minBand == iBand ... % min was reached with current band
%                         &&  inhibitionMap( center(idx(j)).minI, center(idx(j)).minJ )  % min distance might have changed
%                     maxActivityMightHaveChanged(idx(j)) = true;
%                     %                     disp(['Impacted min for cell ' int2str(idx(j))])
%                 end
%             end
%         end
%         %         sum(maxActivityMightHaveChanged(:))
%         for j=i+1:PARAM.n % recompute min distance where necessary
%             if maxActivityMightHaveChanged(idx(j))
%                 center(idx(j)) = computeMaxActivity(center(idx(j)),PARAM.sigma,alpha);
%                 %                 if center(idx(j)).maxActivity > activity(j) %debug
%                 %                     error('activity greater after inhibition. this shouldn''t happen.')
%                 %                 end
%                 activity(j) = center(idx(j)).maxActivity;
%             end
%         end
% 
%         % reorder distances (note: only terms i+1:end may have changed)
%         if sum(maxActivityMightHaveChanged)>0
%             [partialActivity partialIdx] = sort(-activity(i+1:end)); partialActivity = -partialActivity;
%             activity = [ activity(1:i) partialActivity];
%             tmp = idx;
%             for j=1:PARAM.n-i
%                 idx(i+j) = tmp(i+partialIdx(j));
%             end
% 
%             %             % debug
%             %             [activity2 idx2] = sort(-[center.maxActivity]); activity2 = -activity2;
%             %             if sum(activity~=activity2)
%             %                 error('sort pb')
%             %             end
%         end

        i = i+1;
        %         if i>kWTA
        %             break
        %         end
                end
    end % firing criterion

    %     fprintf(1,'.');
    if currentImIdx == nMap % last image
        currentImIdx = 1;
        imOrder = randperm(nMap); % reshuffle
    else
        currentImIdx = currentImIdx+1; % move on to next image
    end
    %             fprintf(1,'\n');

    if PARAM.oneByOne
        if center(nActiveCenter).nFiring>=PARAM.learningPeriod
            nActiveCenter = nActiveCenter + 1;
            if nActiveCenter > PARAM.n
                disp(['All neurons have learned. Exiting.'])
                break
            else
                disp(['Now ' int2str(nActiveCenter) ' neurons in the pool'])
            end
        end
    end


    if mod(iter,10000)==0
        fprintf(1,'.',iter);
        %         reconstruction(center,filters,fSiz,nFilter,['iter ' int2str(iter)])
    end
    if iter==100 && sum([center.nFiring])==0
        disp('Neurons do not fire. Increase distance threshold and/or increase sigma')
        break;
    end

%     for iCenter=1:nActiveCenter
%         center(iCenter).activityThreshold = center(iCenter).activityThreshold * (1-PARAM.thrDecay);
%     end
end % iteration

fprintf(1,'\n');
toc

function center = fire(center,bestPatch,iter,PARAM)
% RFSize = size(bestPatch,1);
% alpha = (PARAM.size/4)^2; % see Mutch CVPR 2006
alpha = 1;

% inc firing counts
center.nFiring = center.nFiring+1;

normedBestPatch =  bestPatch/norm(bestPatch(:));

% Hebbian modification: dW = a.y.(X-W)
y = 1;
% y = (1-center.minDistance/PARAM.firingThreshold/PARAM.size)^1;
% y = exp( -center.minDistance / 2 / center.sigma^2 / alpha );
% y = max(0,1/center.d2Threshold*(center.d2Threshold-center.minDistance));
% y = y^2 * 3;
% newCenterPatch = center.patch + PARAM.learningRate * y * bestPatch;
newCenterPatch = center.patch + center.learningRate * y * ( normedBestPatch - center.patch );
% newCenterPatch = center.patch + center.learningRate * y * ( bestPatch - center.patch );

% normalize
newCenterPatch = newCenterPatch / norm(newCenterPatch(:));

% Reporting
% % compute delta theta
% center.evolDeltaTheta(iter) = acos(min(1,dotProduct(center.patch(:),newCenterPatch(:))));
% % compute Theta
% center.evolTheta(iter) = acos(min(1,dotProduct(center.patch(:),normedBestPatch(:))));
% % dist
% center.evolDist(iter) = center.maxActivity;
% % band
% center.evolBand(iter) = center.minBand;
% % dist
center.evolThetaOriginal(iter) = acos(min(1,dotProduct(newCenterPatch(:),center.originalPatch(:))));

% store new center
center.patch = newCenterPatch;

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
% center.activityThreshold = .95*center.maxActivity;

function center = computeMaxActivity(center,sigma,alpha)

center.minDistance = Inf;
for iBand=1:length(center.distanceMap)
    [distance idx] = min(center.distanceMap{iBand}(:));
    if distance < center.minDistance % found a better match
        center.minDistance = max(0,distance);% just to avoid tiny negative values
        center.maxActivity = exp( -distance / 2 / (sigma^2) / alpha );
        center.minBand = iBand;
        center.minJ = floor((idx-1)/size(center.distanceMap{iBand},1))+1;
        center.minI = idx - (center.minJ-1) * size(center.distanceMap{iBand},1);
    end
end

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

