function rec(centers,filter,title)

if nargin==2
    title='Reconstruction';
end

s = size(centers); % format i_s x j_s x i_onOff x j_onOfff x iCenter x 2

N = s(1)*s(2)*s(5);

figure('Name',title,'MenuBar','none')
nC = round(sqrt(4/3*N));
nR = ceil(N/nC);
margin = .2;
width = nC*(1+margin)+margin;
height = nR*(1+margin)+margin;

n = 0;
for i =1:s(1)
    for j =1:s(2)
        for k = 1:s(5);
            n=n+1;
            r = floor((n-1)/nC)+1;
            c = n - (r-1)*nC;
            subplot('Position',[ ((c-1)*(1+margin)+margin)/width ((nR-r)*(1+margin)+margin)/height 1/width 1/height ])
            patch = reshape(centers(i,j,:,:,k,:),[s(3) s(4) s(6)]);

            % filter options
            %     patch = keepOnlyBestOri(patch.*(patch>mean(patch(:))));
            %     patch = patch.*(patch>mean(patch(:)));

            %     patch = ( patch - min(patch(:)) ) / (max(patch(:)) -min(patch(:)));
            rec = RFReconstruction(patch, filter, false);
            imagesc(rec)
            colormap(gray)
            %     text(0,0,0,int2str(centers(i).nFiring))
            %     pause
            %     title(['nF = ' int2str(centers(i).nFiring)])
            %     xlabel('test')
            axis off
        end
    end
end
