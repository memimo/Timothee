function dp=analyzeCenter(allCenters, name, patchSizes, numPatchesPerSize)

if nargin<=2
    patchSizes=size(allCenters(1).patch,1);
    numPatchesPerSize = length(allCenters);
end
if nargin==1
    name = ''
end

for s=1:length(patchSizes)
    center = allCenters((s-1)*numPatchesPerSize+1:s*numPatchesPerSize);


    nIter = length(center(1).evolThetaOriginal);
    nCenter = length(center);

    set(0,'DefaultAxesLineStyleOrder','-|--|:')
    for c=1:nCenter
        leg{c} = int2str(c);
    end


    figure('Name',[name ' - center analyze (size ' int2str(patchSizes(s)) ')' ],'MenuBar','none')
    nR=1;
    nC=3;
    i=1;

    % subplot(nR,nC,i)
    % dist = reshape([center.evolDist],nIter,nCenter);
    % plot(dist,'.')
    % title('min d(w,x)')
    % i=i+1;

    if isfield(center,'evolTheta')
        subplot(nR,nC,i)
        theta = reshape(180/pi*[center.evolTheta],nIter,nCenter);
        plot(theta,'.')
        title('Theta(w,x)')
        i=i+1;

        subplot(nR,nC,i)
        avgPeriod = round(.8*nIter):nIter;
        avgTheta = sum(theta(avgPeriod,:))./sum(theta(avgPeriod,:)>0);
        hist(avgTheta)
        title(['Mean theta (last ' int2str(length(avgPeriod)) ' iter)'])
        % disp(['Mean theta over last ' int2str(length(avgPeriod)) ' iterations = ' num2str() ])
        i=i+1;

    end

    subplot(nR,nC,i)
    evolThetaOriginal = reshape(180/pi*[center.evolThetaOriginal],nIter,nCenter);
    plot(evolThetaOriginal,'.')
    title('Theta(w,w0)')
    legend(char(leg),'Location','EastOutside');
    i=i+1;

    if isfield(center,'evolDeltaTheta')
        subplot(nR,nC,i)
        deltaTheta = reshape(180/pi*[center.evolDeltaTheta],nIter,nCenter);
        plot(deltaTheta,'.')
        title('d Theta')
        i=i+1;
    end
    if isfield(center,'evolBand')
        subplot(nR,nC,i)
        band = reshape([center.evolBand],nIter,nCenter);
        plot(band,'.')
        title('Scale band')
        i=i+1;
    end
    subplot(nR,nC,i)
    hist([center.nFiring])
    title('nFiring')
    i=i+1;

    
    if isfield(center,'evolThreshold')
        subplot(nR,nC,i)
        thr = reshape([center.evolThreshold],nIter,nCenter);
        plot(thr,'.')
        title('Threshold')
        i=i+1;
    end

    meanDotProd=0;
    dp = zeros(nCenter,nCenter);
    for m=1:nCenter
        if center(m).nFiring == 0 % did not learn
            continue;
        end
        for n=m+1:nCenter
            if center(n).nFiring == 0 % did not learn
                continue;
            end

            %             dotProduct(center(n).patch(:),center(m).patch(:))
            dp(m,n) = dotProduct(center(n).patch(:)/norm(center(n).patch(:)),center(m).patch(:)/norm(center(m).patch(:)));
            dp(n,m) = dp(m,n);
            if dp(m,n) > .75
                %                 warning([int2str(m) ' and ' int2str(n) ' similar : dot prod=' num2str(dp(m,n))])
            end
            meanDotProd = meanDotProd + dp(m,n);
        end
    end
    meanDotProd = meanDotProd / (nCenter*(nCenter-1)/2);
    disp(['Mean dot prod = ' num2str(meanDotProd)])

    [maxDp maxIdx] = max(dp(:));
    n = floor((maxIdx-1)/nCenter)+1;
    m = maxIdx - (n-1)*nCenter;
    disp(['Max dot prod = ' num2str(maxDp) ' (between centers ' int2str(n) ' and ' int2str(m) ')'])
    disp(['Dot prod > .5 = ' num2str( 100 * sum(dp(:)>.5)/2/(nCenter*(nCenter-1)/2) ) '%'])

end
