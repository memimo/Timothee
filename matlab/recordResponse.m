% function [response ] = recordResponse(nMap,getMap,weightSharingMode,center)
function [response ] = recordResponse(range,getMap,weightSharingMode,center)
% Timothee Masquelier timothee.masquelier@alum.mit.edu August 2006
% nMap: number of input maps
% getMap: handle to a function i -> map
% weightSharingMode: if true convolution mode + max response, if false simple local dot product
% center: structure that must contains at least center.patch

% tmp = getMap(1);
% nFilter = size(tmp{1},3);


% offset between image and convolution map
% offset = floor((PARAM.RFSize-1)/2); % use this offset with 'same' convo type
% offset = 0; % use this offset with 'valid' convo type

sz = size(center); % format: i_s x j_s x i_onOff x j_onOfff x iCenter x 2
nCenter = sz(1)*sz(2)*sz(5);

count=1;
for i=1:sz(1)
    for j=1:sz(2)
        for c=1:sz(5)
            tmp = center(i,j,:,:,c,:);
            patch{count} = tmp(:);
            count=count+1;
        end
    end
end

    

response = zeros(nCenter,length(range));
disp(['Recording responses...'])
tic
if weightSharingMode
    error('this does not work yet')
%     flippedCenter = center(iCenter).patch(end:-1:1,end:-1:1,:); % patch is flipped because WindowedPatchDistance uses conv2 (and not filer2)
%     for iBand=1:length(map)
%         center(iCenter).distanceMap{iBand} = WindowedPatchDistance(map{iBand},flippedCenter,'valid');
%     end
%     center(iCenter) = computeMaxActivity(center(iCenter),PARAM.sigma,alpha);
else
    for m=range
        %     map = getMap(m);
        map = feval(getMap,m);
        im = map{1};
        im = im / norm(im(:)); % this normalization is very useful (Tim 04/2007)
        for iCenter=1:nCenter
            response(iCenter,m-range(1)+1) = dotProduct(im(:),patch{iCenter});
        end
    end
    if mod(m-range(1)+1,10000)==0
        fprintf(1,'.');
    end
end
toc
