function mapList = computeOnOffMap(filter,path,cropPos,internalPos,mapSize,subSampling,zoom)

filterSize = size(filter,1);
mask = ones(filterSize,1);


% initialize variables
map = zeros([mapSize 2]);
cropSize = subSampling*(mapSize-1) + (filterSize*ones(1,2));

dirList = dir([path '*.tif']);
nIm = length(dirList);
mapList{nIm*size(cropPos,1)} = {};


disp(['Computing ' int2str(nIm*size(cropPos,1)) ' ON-OFF maps (from ' int2str(nIm) ' images)...'])
tic
for i=1:nIm
    name = dirList(i).name(1:end-4);
    im8 = imread([path dirList(i).name]);
    if zoom<1
        im8 = shrink(im8,zoom);
    end

    for j=1:size(cropPos,1)

        % tmp: crop
%        disp(int2str(   [ cropPos(j,1)+internalPos(1)-1-floor(.5*(cropSize(1)-1)),cropPos(j,1)+internalPos(1)-1+ceil(.5*(cropSize(1)-1));
%                             cropPos(j,2)+internalPos(2)-1-floor(.5*(cropSize(2)-1)),cropPos(j,2)+internalPos(2)-1+ceil(.5*(cropSize(2)-1))]));
%                         pause
        imCropped = im8(    cropPos(j,1)+internalPos(1)-1-floor(.5*(cropSize(1)-1)):cropPos(j,1)+internalPos(1)-1+ceil(.5*(cropSize(1)-1)), ...
                            cropPos(j,2)+internalPos(2)-1-floor(.5*(cropSize(2)-1)):cropPos(j,2)+internalPos(2)-1+ceil(.5*(cropSize(2)-1)));
        im = double(imCropped) / 255;

        map = conv2(im,filter,'valid');
        map = map(1:subSampling:end,1:subSampling:end); % sub-sampling
    
%         % normalization
%         normMap = conv2(mask,mask,im.^2,'valid') .^ 0.5;
%         normMap = normMap + ~normMap; % avoid dividing by zero. anyway normMap=0 => im=0 => map=0
%         normMap = normMap(1:subSampling:end,1:subSampling:end); % sub-sampling
% 
%         map = map ./ normMap;

%         onMap = max(0,map);
%         offMap = - min(0,map);
%     
%         map(:,:,1) = onMap;
%         map(:,:,2) = offMap;
        
        mapList{(j-1)*nIm+i} = map;
        
    end %j
    
    if mod(i,10000)==0
        fprintf(1,'.');
    end
end
toc
