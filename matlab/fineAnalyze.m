
margin =.1;
margin_loc =.1;
selectedBound = .9;


for n=1:size(COMMON.weight,4)
    nActive = sum(sum(sum(COMMON.weight(:,:,:,n)>selectedBound)));
    if nActive>0
%         figure('Name',['Cell ' int2str(n) ' - ' PARAM.compName],'MenuBar','none')
        figure('Name',['Cell ' int2str(n) ' - ' PARAM.compName],'PaperPositionMode','auto','Position',4*[1 1 (4*12+2) (4+4*12+2)])
        nR = size(COMMON.weight,1);
        nC = size(COMMON.weight,2);
        width = nC*(1+margin)+margin;
        height = nR*(1+margin)+margin;

        nActive_loc_max = max(max(sum(COMMON.weight(:,:,:,n)>selectedBound,3)));
        nC_loc = 4;%round(sqrt(4/3*nActive_loc_max));
        nR_loc = 4;%ceil(nActive_loc_max/nC_loc);
        width_loc = nC_loc*(1+margin_loc)+margin_loc;
        height_loc = nR_loc*(1+margin_loc)+margin_loc;

        for i=1:size(COMMON.weight,1)
%             if i<size(COMMON.weight,1)
%                 subplot('Position',[ 0 ((nR-i)*(1+margin)+.5*margin)/height 1 margin/10/height ])
%             end

            for j=1:size(COMMON.weight,2)
%                 if j<size(COMMON.weight,2)
%                     subplot('Position',[((j-1)*(1+margin)+.5*margin)/width 0 margin/10/height  1 ])
%                 end

                framePos = [((j-1)*(1+margin)+margin)/width ((nR-i)*(1+margin)+margin)/height];% left, bottom
                frameSize = [1/width 1/height];
                
%                 rectangle('Position',[ framePos frameSize])
%                 pause
                
                k = 1;
                for f=1:size(COMMON.weight,3)
                    if COMMON.weight(i,j,f,n)>selectedBound
                        r = floor((k-1)/nC_loc)+1;
                        c = k - (r-1)*nC_loc;
                        pos = framePos + [ ((c-1)*(1+margin_loc)+margin_loc)/width_loc/width ((nR_loc-r)*(1+margin_loc)+margin_loc)/height_loc/height ];
                        sze = frameSize .* [1/width_loc 1/height_loc];
                        subplot('Position',[ pos sze ])

                        %                         subplot(nR,nC,k)
                        rec = RFReconstruction(reshape(COMMON.center(i,j,:,:,f,:),[PARAM.s.RFSize 2]), PARAM.recFilter, false); % format i_s x j_s x i_onOff x j_onOfff x iCenter x 2
                        imagesc(rec)
                        colormap(gray)
                        axis off

%                         k=k+1;
                    end
                    k=k+1;
                end

            end
        end
        print('-deps',['C1_' sprintf('%02d',n) ]);
    end
end

% rec all

% figure('Name',['All cells - ' PARAM.compName],'MenuBar','none')
figure('Name',['All cells - ' PARAM.compName],'PaperPositionMode','auto','Position',4*[1 1 (4*12+2) (4+4*12+2)])
nR = size(COMMON.weight,1);
nC = size(COMMON.weight,2);
width = nC*(1+margin)+margin;
height = nR*(1+margin)+margin;

nActive_loc_max = size(COMMON.weight,3);
nC_loc = round(sqrt(1*nActive_loc_max));
nR_loc = ceil(nActive_loc_max/nC_loc);
width_loc = nC_loc*(1+margin_loc)+margin_loc;
height_loc = nR_loc*(1+margin_loc)+margin_loc;

for i=1:size(COMMON.weight,1)
    for j=1:size(COMMON.weight,2)

        framePos = [((j-1)*(1+margin)+margin)/width ((nR-i)*(1+margin)+margin)/height];% left, bottom
        frameSize = [1/width 1/height];

        k = 1;
        for f=1:size(COMMON.weight,3)
            r = floor((k-1)/nC_loc)+1;
            c = k - (r-1)*nC_loc;
            pos = framePos + [ ((c-1)*(1+margin_loc)+margin_loc)/width_loc/width ((nR_loc-r)*(1+margin_loc)+margin_loc)/height_loc/height ];
            sze = frameSize .* [1/width_loc 1/height_loc];
            subplot('Position',[ pos sze ])

            %                         subplot(nR,nC,k)
            rec = RFReconstruction(reshape(COMMON.center(i,j,:,:,f,:),[PARAM.s.RFSize 2]), PARAM.recFilter, false); % format i_s x j_s x i_onOff x j_onOfff x iCenter x 2
            imagesc(rec)
            colormap(gray)
            axis off

            k=k+1;
        end
    end
end

print('-deps','All_S1');



% unselected
figure('Name',['Unselected - ' PARAM.compName],'PaperPositionMode','auto','Position',4*[1 1 (4*12+2) (4+4*12+2)])
nR = size(COMMON.weight,1);
nC = size(COMMON.weight,2);
width = nC*(1+margin)+margin;
height = nR*(1+margin)+margin;

nActive_loc_max = size(COMMON.weight,3);
nC_loc = round(sqrt(1*nActive_loc_max));
nR_loc = ceil(nActive_loc_max/nC_loc);
width_loc = nC_loc*(1+margin_loc)+margin_loc;
height_loc = nR_loc*(1+margin_loc)+margin_loc;

for i=1:size(COMMON.weight,1)
    for j=1:size(COMMON.weight,2)

        framePos = [((j-1)*(1+margin)+margin)/width ((nR-i)*(1+margin)+margin)/height];% left, bottom
        frameSize = [1/width 1/height];

        k = 1;
        for f=1:size(COMMON.weight,3)

            if max(COMMON.weight(i,j,f,:))<=selectedBound % unselected criterion
                r = floor((k-1)/nC_loc)+1;
                c = k - (r-1)*nC_loc;
                pos = framePos + [ ((c-1)*(1+margin_loc)+margin_loc)/width_loc/width ((nR_loc-r)*(1+margin_loc)+margin_loc)/height_loc/height ];
                sze = frameSize .* [1/width_loc 1/height_loc];
                subplot('Position',[ pos sze ])

                %                         subplot(nR,nC,k)
                rec = RFReconstruction(reshape(COMMON.center(i,j,:,:,f,:),[PARAM.s.RFSize 2]), PARAM.recFilter, false); % format i_s x j_s x i_onOff x j_onOfff x iCenter x 2
                imagesc(rec)
                colormap(gray)
                axis off
            end
            k = k+1;

        end
    end
end

print('-deps','Unselected');





