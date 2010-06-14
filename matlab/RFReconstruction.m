function reconstruction = RFReconstruction(weight,convo,polarized)

[nR nC nOri] = size(weight);

reconstruction = zeros(nR, nC);


if polarized
    for iOri = 1:nOri/2
        tmp = conv2(weight(:,:,iOri), convo(:,:,iOri), 'same');
        reconstruction	= reconstruction + tmp;
    end
    
    for iOri = nOri/2+1:nOri
        tmp = conv2(weight(:,:,iOri), -convo(:,:,iOri-nOri/2), 'same');
        reconstruction	= reconstruction + tmp;
    end
else
    for iOri = 1:nOri
        tmp = conv2(weight(:,:,iOri), convo(:,:,iOri), 'same');
        reconstruction	= reconstruction + tmp;
    end
    % tim 04/2006 : in the convo kernels it is easier to define positive value on edges
    reconstruction = -reconstruction;
end
