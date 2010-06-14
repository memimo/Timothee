%***********************
% SIMPLE CELL LEARNING *
%***********************

rand('state',PARAM.randState);

for i=1:PARAM.c.RFSize(1)
    for j=1:PARAM.c.RFSize(2)
        
        timedLogLn(['S1 MAP (' int2str(i) ',' int2str(j) ')']);

        load([PATH.map 'onOff.' int2str(i) '.' int2str(j) '.mat']);
        COMMON.onOffMap = onOffMap;
        clear onOffMap;

        % LEARNING
        if true % learn S1 selectivity
            %             [center, nCouldHaveFired, nFired] = competitiveHebbianCenterLearning(length(COMMON.onOffMap),PARAM.s.nIter,@getOnOffMap,false,[],[],PARAM.recFilter,PARAM.s,PARAM.cropPos);
            %             [center, nCouldHaveFired, nFired] = competitiveHebbianCenterLearning(length(COMMON.onOffMap),PARAM.s.nIter,@getOnOffMap,false,[],[],PARAM.recFilter,PARAM.s);
%                         [center, nCouldHaveFired, nFired, w ] = antiHebbianLearning(length(COMMON.onOffMap),PARAM.s.nIter,@getOnOffMap,[],[],[],PARAM.recFilter,PARAM.s);
            [center] = learnUnlearn(length(COMMON.onOffMap),PARAM.s.nIter,@getOnOffMap,[],[],PARAM.recFilter,PARAM.s);
%                                 rec(center,PARAM.recFilter,'')
%                                 dp = analyzeCenter(center,'');
%                                 error('stop')

            % free memory
            if isfield(center,'evolDist')
                center = rmfield(center,'evolDist');
            end
            if isfield(center,'evolThetaOriginal')
                center = rmfield(center,'evolThetaOriginal');
            end
            if isfield(center,'evolDeltaTheta')
                center = rmfield(center,'evolDeltaTheta');
            end
            if isfield(center,'evolTheta')
                center = rmfield(center,'evolTheta');
            end
            if isfield(center,'evolBand')
                center = rmfield(center,'evolBand');
            end
            clear nFired nCouldHaveFired

        else % use hard wired gabor filter

            [fSiz,filters,c1OL,numSimpleFilters] = init_gabor(180/(PARAM.s.n/2)*[0:PARAM.s.n/2-1], PARAM.s.RFSize(1), .5, .5);
            for r=1:(PARAM.s.n/2)

                filter(:,:,r) = reshape(filters(:,r),[PARAM.s.RFSize]);
                center(r).patch(:,:,1) = max(0,filter(:,:,r));
                center(r).patch(:,:,2) = -min(0,filter(:,:,r));
                center(r+PARAM.s.n/2).patch(:,:,1) = -min(0,filter(:,:,r));
                center(r+PARAM.s.n/2).patch(:,:,2) = max(0,filter(:,:,r));


                %                 figure
                %                 imagesc(filter(:,:,r))
                %                 colormap(gray)
                %                 disp(filter(:,:,r))
                %                 pause
                %         filter(:,:,r) = filter(end:-1:1,end:-1:1,r); % flip in order to use conv2
            end
        end
        % format i_s x j_s x i_onOff x j_onOfff x iCenter x 2
        COMMON.center(i,j,:,:,:,:) = reshape([center.patch],[PARAM.s.RFSize length(center) 2]);

%         recAll(COMMON.center,PARAM.recFilter,[ 'All S1' PARAM.compName])
%         error('stop')

%        error('stop')

        % RECORD RESPONSE
        for block=1:ceil( length(COMMON.onOffMap) / PARAM.c.block )
            name = ['s1.' sprintf('%02d',block) '.' PARAM.s.type '.mat'];
            if i==1 && j==1
                s1Map = zeros([PARAM.c.RFSize PARAM.s.n PARAM.c.block]);
            else
                load([PATH.map name])
            end
            range = 1+(block-1)*PARAM.c.block:min(length(COMMON.onOffMap),block*PARAM.c.block);
%             s1Map(i,j,:,1:length(range)) = s1(:,range);% format i_s x j_s x iCenter x iMap
            s1Map(i,j,:,1:length(range)) = recordResponse(range,@getOnOffMap,false,COMMON.center(i,j,:,:,:,:));% format iCenter x iMap
            save([PATH.map name],'s1Map');
            clear s1Map
        end
        COMMON.nMap = length(COMMON.onOffMap);
        COMMON = rmfield(COMMON,'onOffMap');
%         clear s1
    end
end
save([PATH.map 'COMMON.' PARAM.compName '.' PARAM.s.type '.mat'],'COMMON','PARAM')
currentPath = pwd;
cd(PATH.map)
copyfile(['COMMON.' PARAM.compName '.' PARAM.s.type '.mat'],['../COMMON.mat'])
cd(currentPath); clear currentPath;

recAll(COMMON.center,PARAM.recFilter,[ 'All S1' PARAM.compName])

% error('stop')
