
if isfield(COMMON,'evol')
    figure('Name',['Weight sum - ' PARAM.compName],'MenuBar','none')
    for n=1:PARAM.c.n
        subplot(1,PARAM.c.n,n)
        %     plot(reshape(evol(:,n,:),[size(evol,1) size(evol,3)])') % note: format of evol: feat x node x iter
        plot(reshape(sum(COMMON.evol(:,n,:),1),[1 size(COMMON.evol,3)])); % note: format of COMMON.evol: feat x node x iter
        %     print('-djpeg',[PATH.result 'rec.' PARAM.compName '.jpg'])
    end
end

figure('Name',['Weight distrib - ' PARAM.compName],'MenuBar','none')
hist(COMMON.weight(:))

%reporting
if isfield(COMMON,'nFiringOutput')
    disp([num2str(100*sum(COMMON.nFiringOutput)/PARAM.c.nIter) '% effective stim' ])
    disp(['Outputs fired ' num2str(mean(COMMON.nFiringOutput(:))) ' times on avg'])
    f=COMMON.nFiringOutput / mean(COMMON.nFiringOutput(:));
    disp(['Std dev among output firing rates = ' num2str(100*std(f(:))) '%' ])
end
if isfield(COMMON,'nFiringInput')
    disp(['Inputs fired ' num2str(mean(COMMON.nFiringInput(:))) ' times on avg'])
    f=COMMON.nFiringInput / mean(COMMON.nFiringInput(:));
    disp(['Std dev among input firing rates = ' num2str(100*std(f(:))) '%' ])    
end


% COMMON.weight format : i x j x feat x node
% criterion = .95;
% criterion = min(COMMON.weight(:)) + (max(COMMON.weight(:))-min(COMMON.weight(:)))/2;
selectedBound = .9;
% unselectedBound = .1;

margin = .2;
% disp([int2str(sum(COMMON.weight(:)>criterion)) ' selected S1 (out of ' int2str(size(COMMON.weight,1)*size(COMMON.weight,1)*size(COMMON.weight,)) ')'])
for n=1:size(COMMON.weight,4)
    nActive = sum(sum(sum(COMMON.weight(:,:,:,n)>selectedBound)));
    disp(['Complex cell ' int2str(n) ': ' int2str(nActive) ' afferents'])
    if nActive>0
        nC = round(sqrt(4/3*nActive));
        nR = ceil(nActive/nC);
%         figure('Name',['Cell ' int2str(n) ' - ' PARAM.compName],'MenuBar','none','PaperPositionMode','auto','Position',4*[1 1 (nC*12+2) (4+nR*12+2)])
        figure('Name',['Cell ' int2str(n) ' - ' PARAM.compName],'MenuBar','none')
        width = nC*(1+margin)+margin;
        height = nR*(1+margin)+margin;
        k = 1;
        for i=1:size(COMMON.weight,1)
            for j=1:size(COMMON.weight,2)
                for f=1:size(COMMON.weight,3)
                    if COMMON.weight(i,j,f,n)>selectedBound
                        r = floor((k-1)/nC)+1;
                        c = k - (r-1)*nC;
                        subplot('Position',[ ((c-1)*(1+margin)+margin)/width ((nR-r)*(1+margin)+margin)/height 1/width 1/height ])

                        %                         subplot(nR,nC,k)
                        rec = RFReconstruction(reshape(COMMON.center(i,j,:,:,f,:),[PARAM.s.RFSize 2]), PARAM.recFilter, false); % format i_s x j_s x i_onOff x j_onOfff x iCenter x 2
                        imagesc(rec)
                        colormap(gray)
                        axis off
                        

                        k=k+1;
                    end
                end
            end
        end
%         print('-dtiff',['C1_' sprintf('%02d',n) ]);
    end
end

% COMMON.weight format : i x j x feat x node
nS1 = size(COMMON.weight,1)*size(COMMON.weight,2)*size(COMMON.weight,3);
selected = sum(COMMON.weight>selectedBound,4);
overSelected = sum(selected(:)>1);
unselected = sum(selected(:)==0);
disp([int2str(overSelected) ' S1 selected by more than 1 C1'])
disp([int2str(unselected) ' S1 not selected by any C1'])
disp([int2str(sum(selected(:)==1)) ' S1 selected by exactly 1 C1'])

if unselected > 0
    figure('Name',['Unselected S1'],'MenuBar','none')
    nC = round(sqrt(4/3*unselected));
    nR = ceil(unselected/nC);
    width = nC*(1+margin)+margin;
    height = nR*(1+margin)+margin;
    k = 1;
    for i=1:size(COMMON.weight,1)
        for j=1:size(COMMON.weight,2)
            for f=1:size(COMMON.weight,3)
                if sum(COMMON.weight(i,j,f,:)>selectedBound)==0
                    r = floor((k-1)/nC)+1;
                    c = k - (r-1)*nC;
                    subplot('Position',[ ((c-1)*(1+margin)+margin)/width ((nR-r)*(1+margin)+margin)/height 1/width 1/height ])
                    rec = RFReconstruction(reshape(COMMON.center(i,j,:,:,f,:),[PARAM.s.RFSize 2]), PARAM.recFilter, false); % format i_s x j_s x i_onOff x j_onOfff x iCenter x 2
                    imagesc(rec)
                    colormap(gray)
                    axis off
                    k=k+1;
                end
            end
        end
    end
end

if overSelected > 0
    figure('Name',['Overselected S1'],'MenuBar','none')
    nC = round(sqrt(4/3*overSelected));
    nR = ceil(overSelected/nC);
    width = nC*(1+margin)+margin;
    height = nR*(1+margin)+margin;
    k = 1;
    for i=1:size(COMMON.weight,1)
        for j=1:size(COMMON.weight,2)
            for f=1:size(COMMON.weight,3)
                if sum(COMMON.weight(i,j,f,:)>selectedBound)>1
                    r = floor((k-1)/nC)+1;
                    c = k - (r-1)*nC;
                    subplot('Position',[ ((c-1)*(1+margin)+margin)/width ((nR-r)*(1+margin)+margin)/height 1/width 1/height ])
                    rec = RFReconstruction(reshape(COMMON.center(i,j,:,:,f,:),[PARAM.s.RFSize 2]), PARAM.recFilter, false); % format i_s x j_s x i_onOff x j_onOfff x iCenter x 2
                    imagesc(rec)
                    colormap(gray)
                    axis off
                    k=k+1;
                end
            end
        end
    end
end

if isfield(COMMON,'nFiringInput')
    figure('Name','nFiringInput - Distribution');hist(COMMON.nFiringInput(:));
end
if isfield(COMMON,'nFiringOutput')
    figure('Name','nFiringOutput - Distribution');hist(COMMON.nFiringOutput(:));
end