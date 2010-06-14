function mapList = computeOnOffMap(filter,path,pos,size)

mask = ones(size,1);


% initialize variables
map = zeros([size 2]);
sf = size(filter);
sf = sf(1:2);
cropSize = size + (sf-1);

dirList = dir(path);
nIm = length(dirList);
mapList{nIm*size(PARAM.cropPos,1)} = {};

disp(['Computing ' int2str(nIm*size(PARAM.cropPos,1)) ' S1 maps (from ' int2str(nIm) ' images)...'])
tic
for i=1:nIm
    name = dirList(i).name(1:end-4);
    im8 = imread([PATH.image dirList(i).name]);

    for j=1:size(PARAM.cropPos,1)

        % tmp: crop
        imCropped = im8(PARAM.cropPos(j,1):PARAM.cropPos(j,1)+cropSize(1)-1,PARAM.cropPos(j,2):PARAM.cropPos(j,2)+cropSize(2)-1);
        im = double(imCropped) / 255;


        %     % USE CONV2 AND THEN SUBSAMPLE
        %     % normalization
        %     normMap = conv2(mask,mask,im.^2,'valid') .^ 0.5;
        %     normMap = normMap + ~normMap; % avoid dividing by zero. anyway normMap=0 => im=0 => map=0
        %     % subsampling
        %     normMap = normMap(1:PARAM.s1SubSampling:end,1:PARAM.s1SubSampling:end);
        %
        %     for r=1:length(PARAM.rot);
        %         rawMap = conv2(im,filter(:,:,r),'valid');
        %         % subsampling
        %         rawMap = rawMap(1:PARAM.s1SubSampling:end,1:PARAM.s1SubSampling:end);
        %         map(:,:,r) = abs( rawMap ./ normMap );
        %     end

        % ONLY COMPUTE NEEDED VALUES
        for I=1:PARAM.s1MapSize(1)
            beginI = (I-1)*PARAM.s1SubSampling+1;
            endI = beginI+PARAM.RF_siz-1;
            for J=1:PARAM.s1MapSize(2)
                beginJ = (J-1)*PARAM.s1SubSampling+1;
                endJ = beginJ+PARAM.RF_siz-1;
                patch = im(beginI:endI,beginJ:endJ);
                normPatch = ( sum(patch(:).^2) )^.5;
                for r=1:length(PARAM.rot);
                    map(I,J,r) = abs(sum(sum(filter(:,:,r).*patch))) / normPatch;
                end
            end
        end
%         save([PATH.map 's1.' sprintf('%06d',(j-1)*nIm+i) '.mat'],'map');
        mapList{(j-1)*nIm+i} = map;
    end
    if mod(i,1000)==0
        fprintf(1,'.');
    end
end
toc
