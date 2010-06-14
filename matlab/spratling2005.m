function [weight evol nFiringInput nFiringOutput]=spratling2005(getMap,nIter,nFeat,nNode,patchSize,weight)
% Timothee Masquelier timothee.masquelier@alum.mit.edu August 2006
% see Spratling 2005
rand('twister',5489);
% global PARAM

PARAM.learningRate = .001;

if length(weight)==0
    for n=1:nNode
        % format: i x j x feat x node
        weight(:,:,:,n) = getInitialWeight(patchSize,nFeat);
    end
    formerY = ones(nNode);
    meanFormerY = mean(formerY);
    
    %reporting
    nFiringOutput = zeros(1,nNode);
    nFiringInput = zeros([patchSize nFeat]);
end

disp('Spratling 2005...')
tic
offset = floor(rand*length(nIter));

for i=1:nIter
%     x = getMap(i); % format: i x j x feat
    x = feval(getMap,i+offset); % format: i x j x feat

    %     x = tmp(patchSize,nFeat);
    nodeMax = [ max(max(max(weight,[],1),[],2),[],3) ];
    nodeMax = nodeMax(:);
    inputMax = max(weight,[],4);
    inputMax = inputMax + ~inputMax; % avoid dividing by 0. anyway inputMax(i,j,f)=0 => for all nodes n weight(i,j,f,n)=0 => for all nodes Z(i,j,f,n)=0
    for n=1:nNode
        if nodeMax(n)==0
            error('nodeMax(n)==0');
        end
        if sum(sum(inputMax==0))>0
            error('sum(sum(inputMax==0))>0');
        end
        Z = x .* ( weight(:,:,:,n).^2 ) ./ inputMax / nodeMax(n); % format: i x j x feat
        %         Z = x .* ( weight(:,:,:,n).^2 ) ;
        [y(n) winner] = max(Z(:));

        if formerY(n) > meanFormerY % only active node learn
            % decrease everybody (including the winner)
            weight(:,:,:,n) = weight(:,:,:,n) - PARAM.learningRate * (formerY(n)-meanFormerY) / (nNode*meanFormerY) * x ;
            % increase the winner twice
            winnerF = floor((winner-1)/patchSize(1)/patchSize(2))+1;
            winnerJ = floor((winner-(winnerF-1)*patchSize(1)*patchSize(2)-1)/patchSize(1))+1;
            winnerI = winner - (winnerF-1)*patchSize(1)*patchSize(2) - (winnerJ-1)*patchSize(1);
            weight(winnerI,winnerJ,winnerF,n) = weight(winnerI,winnerJ,winnerF,n) + 2 * PARAM.learningRate * x(winnerI,winnerJ,winnerF) * (formerY(n)-meanFormerY) / (nNode*meanFormerY);
           
            % inc counters
            nFiringInput(winnerI,winnerJ,winnerF) = nFiringInput(winnerI,winnerJ,winnerF)+ 1;
            nFiringOutput(n) = nFiringOutput(n) + 1;
        end
    end
    
    % clip
    weight = max(weight,0);
    % normalize
    normFactor = sum(weight,4);
    for n=1:nNode
        %                 weight(:,:,:,n) = weight(:,:,:,n) ./ sum(sum(sum(weight(:,:,:,n))));
        weight(:,:,:,n) = weight(:,:,:,n) ./ normFactor;
    end

    if mod(i,100)==0
%         evol(:,:,i/100) = reshape(sum(sum(weight,1),2),[PARAM.nFeat PARAM.n]); % format: feat x node x iter
        evol(1,:,i/100) = sum(sum(sum(weight,1),2),3); % format: 1 x node x iter
    end

    formerY = y;
    meanFormerY = mean(formerY);
    if mod(i,10000)==0
        fprintf(1,'.');
    end
end
toc
% fprintf(1,'\n');

function weight = getInitialWeight(patchSize,nFeat)
weight = 1 + .01*rand([patchSize nFeat]);
weight = max(0,weight);
weight = weight / sum(weight(:));

% function input = tmp(patchSize,nFeat)
% input = rand(patchSize,patchSize,nFeat);
