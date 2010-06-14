function [weight evol thresholdAM nFiringInput nFiringOutput]=einhauser2002(getMap,nIter,nFeat,nNode,patchSize,weight)
% Timothee Masquelier timothee.masquelier@alum.mit.edu August 2006
% see Einhauser 2002
% global PARAM
rand('twister',5489);
PARAM.nu = Inf;
% PARAM.thresoldAM = .05;

% % Einhauser : 
% PARAM.learningRate = 2^-5;
% PARAM.thrDecay = 10^-4;
% Masquelier :
PARAM.learningRate = 5*10^-4;
PARAM.thrDecay = 2^-15;

if length(weight)==0
    for n=1:nNode
        % format: i x j x feat x node
        weight(:,:,:,n) = getInitialWeight(patchSize,nFeat);
    end
    % init middle layer
    if PARAM.nu < Inf
        trAM = 1/PARAM.nu*ones([patchSize nFeat]);% format: i x j x feat
    else
        trAM = ones([patchSize nFeat]);% format: i x j x feat
    end
    thresholdAM = .001*ones([patchSize nFeat]);% format: i x j x feat
    previousMaxAM = 0;
    previousWinnerMI = 1; % has no impact
    previousWinnerMJ = 1; % has no impact
    previousWinnerMF = 1; % has no impact
    
    % init top layer
    if PARAM.nu < Inf
        trAT = 1/PARAM.nu*ones(1,nNode);
    else
        trAT = ones(1,nNode);
    end
    
    previousMaxAT = -1;
    previousWinnerT = -1;
    
    %reporting
    nFiringOutput = zeros(1,nNode);
    nFiringInput = zeros([patchSize nFeat]);
end

disp('Einhauser 2002...')
tic
offset = floor(rand*length(nIter));
for i=1:nIter
%     x = getMap(i); % format: i x j x feat
    x = feval(getMap,i+offset); % format: i x j x feat
        
    AM = x ./ trAM; % format: i x j x feat
%     AM = x; % format: i x j x feat

    [maxAM winnerM] = max(AM(:));
    winnerMF = floor((winnerM-1)/patchSize(1)/patchSize(2))+1;
    winnerMJ = floor((winnerM-(winnerMF-1)*patchSize(1)*patchSize(2)-1)/patchSize(1))+1;
    winnerMI = winnerM - (winnerMF-1)*patchSize(1)*patchSize(2) - (winnerMJ-1)*patchSize(1);
    
%     % Spratling normalization
%     normFactor = sum(weight,4);
%     for n=1:nNode
%         %                 weight(:,:,:,n) = weight(:,:,:,n) ./ sum(sum(sum(weight(:,:,:,n))));
%         weight(:,:,:,n) = weight(:,:,:,n) ./ normFactor;
%     end
%     % Spratling : compute Z
%     nodeMax = [ max(max(max(weight,[],1),[],2),[],3) ];
%     nodeMax = nodeMax(:);
%     inputMax = max(weight,[],4);
%     inputMax = inputMax + ~inputMax; % avoid dividing by 0. anyway inputMax(i,j,f)=0 => for all nodes n weight(i,j,f,n)=0 => for all nodes Z(i,j,f,n)=0
%     for n=1:nNode
%         Z = x .* ( weight(:,:,:,n).^2 ) ./ inputMax / nodeMax(n); % format: i x j x feat
%         AT(n) = max(max(max( Z ))) / trAT(n);
%     end    
    for n=1:nNode
        AT(n) = max(max(max( weight(:,:,:,n).*AM ))) / trAT(n);
%         AT(n) = max(max(max( weight(:,:,:,n).*AM )));
    end
    % find winner in top layer
    [maxAT winnerT] = max(AT);

    % learning
%     if previousMaxAM > thresholdAM(previousWinnerMI,previousWinnerMJ,previousWinnerMF)
%         % update thr
%         thresholdAM(previousWinnerMI,previousWinnerMJ,previousWinnerMF) = previousMaxAM;
%         % update weights:
%         % depress everyone
%         weight(:,:,:,winnerT) = (1-PARAM.learningRate) * weight(:,:,:,winnerT);
%         % potentiate synapse between winners
%         weight(previousWinnerMI,previousWinnerMJ,previousWinnerMF,winnerT) =  weight(previousWinnerMI,previousWinnerMJ,previousWinnerMF,winnerT) + PARAM.learningRate;
% %         % normalize
% %         weight(:,:,:,winnerT) = weight(:,:,:,winnerT) / sum(sum(sum(weight(:,:,:,winnerT))));
%         
%         % reporting
%         nFiringInput(previousWinnerMI,previousWinnerMJ,previousWinnerMF) = nFiringInput(previousWinnerMI,previousWinnerMJ,previousWinnerMF)+ 1;
%         nFiringOutput(winnerT) = nFiringOutput(winnerT) + 1;
% %         if sum(nFiringInput(:)>0) == length(nFiringInput(:))
% %             error('all input have fired at least once')
% %         end
%     end
    % learning
    if i>1 && maxAM > thresholdAM(winnerMI,winnerMJ,winnerMF)
        % update thr
        thresholdAM(winnerMI,winnerMJ,winnerMF) = maxAM;
        % update weights:
        % depress everyone
        weight(:,:,:,previousWinnerT) = (1-PARAM.learningRate) * weight(:,:,:,previousWinnerT);
        % potentiate synapse between winners
        weight(winnerMI,winnerMJ,winnerMF,previousWinnerT) =  weight(winnerMI,winnerMJ,winnerMF,previousWinnerT) + PARAM.learningRate;
%         % normalize
%         weight(:,:,:,winnerT) = weight(:,:,:,winnerT) / sum(sum(sum(weight(:,:,:,winnerT))));
        
        % reporting
        nFiringInput(previousWinnerMI,previousWinnerMJ,previousWinnerMF) = nFiringInput(previousWinnerMI,previousWinnerMJ,previousWinnerMF)+ 1;
        nFiringOutput(winnerT) = nFiringOutput(winnerT) + 1;
%         if sum(nFiringInput(:)>0) == length(nFiringInput(:))
%             error('all input have fired at least once')
%         end
    end
    
    
    if mod(i,100)==0
%         evol(:,:,i/100) = reshape(sum(sum(weight,1),2),[PARAM.nFeat PARAM.n]); % format: feat x node x iter
        evol(1,:,i/100) = sum(sum(sum(weight,1),2),3); % format: 1 x node x iter
    end

    % update traces
    trAM = AM / PARAM.nu + (1-1/PARAM.nu) * trAM;
    trAT = AT / PARAM.nu + (1-1/PARAM.nu) * trAT;
    % thr decay
    thresholdAM = thresholdAM*(1-PARAM.thrDecay);
    % previous := current
    previousMaxAM = maxAM;
    previousWinnerMI = winnerMI;
    previousWinnerMJ = winnerMJ;
    previousWinnerMF = winnerMF;
    previousWinnerT = winnerT;
    previousMaxAT = maxAT;

    if mod(i,10000)==0
        fprintf(1,'.');
    end
end
toc
fprintf(1,'\n');

function weight = getInitialWeight(patchSize,nFeat)
% weight = 1 + .01/100*rand([patchSize nFeat]);
% weight = max(0,weight);
% weight = weight / sum(weight(:));
weight = .75*( 1 + .0 * (rand([patchSize nFeat])-.5) );

% function input = tmp(patchSize,nFeat)
% input = rand(patchSize,patchSize,nFeat);
