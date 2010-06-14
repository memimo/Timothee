function [weight evol]=foldiak1991(getMap,weight,PARAM)
% Timothee Masquelier timothee.masquelier@alum.mit.edu Sept 2006
% see Foldiak 1991
rand('twister',5489);
PARAM.nu = 5;
PARAM.learningRate = .02;

if length(weight)==0
    for n=1:PARAM.n
        % format: i x j x feat x node
        weight(:,:,:,n) = getInitialWeight(PARAM.RFSize,PARAM.nFeat);
    end
    
    % init top layer
    trY = 1/PARAM.nu*ones(1,PARAM.n);
    
end

disp('Foldiak 1991...')
tic
offset = floor(rand*length(PARAM.nIter));
for i=1:PARAM.nIter
%     x = getMap(i); % format: i x j x feat
    x = feval(getMap,i+offset); % format: i x j x feat
        
    for n=1:PARAM.n
%         y(n) = max(max(max( weight(:,:,:,n).*x )));
        y(n) = sum(sum(sum( weight(:,:,:,n).*x )));
    end
 
    % find winner in top layer
    [maxY winner] = max(y);

    % update winner's weights
    weight(:,:,:,winner) = weight(:,:,:,winner) + PARAM.learningRate * trY(winner)*(x-weight(:,:,:,winner));
    

    
    if mod(i,100)==0
%         evol(:,:,i/100) = reshape(sum(sum(weight,1),2),[PARAM.nFeat PARAM.n]); % format: feat x node x iter
        evol(1,:,i/100) = sum(sum(sum(weight,1),2),3); % format: 1 x node x iter
    end

    % update traces
    trY = (1:PARAM.n==winner) / PARAM.nu + (1-1/PARAM.nu) * trY;

    if mod(i,10000)==0
        fprintf(1,'.');
    end
end
toc
fprintf(1,'\n');

function weight = getInitialWeight(patchSize,nFeat)
weight = 0 + .1*rand([patchSize nFeat]);
% weight = max(0,weight);
% weight = 103 * weight / sum(weight(:));

% function input = tmp(patchSize,nFeat)
% input = rand(patchSize,patchSize,nFeat);
