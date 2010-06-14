function rec(centers,filter,title)

if nargin==2
    title='Reconstruction';
end

figure('Name',title,'MenuBar','none')
nC = round(sqrt(4/3*length(centers)));
nR = ceil(length(centers)/nC);
margin = .2;
width = nC*(1+margin)+margin;
height = nR*(1+margin)+margin;

for i =1:length(centers)
    r = floor((i-1)/nC)+1;
    c = i - (r-1)*nC;
    subplot('Position',[ ((c-1)*(1+margin)+margin)/width ((nR-r)*(1+margin)+margin)/height 1/width 1/height ])
    patch = centers(i).patch;
    
    % filter options
%     patch = keepOnlyBestOri(patch.*(patch>mean(patch(:))));
%     patch = patch.*(patch>mean(patch(:)));
    
%     patch = ( patch - min(patch(:)) ) / (max(patch(:)) -min(patch(:))); 
    rec = RFReconstruction(patch, filter, false);
    imagesc(rec)
    colormap(gray)
    text(0,0,0,int2str(centers(i).nFiring))
%     pause
%     title(['nF = ' int2str(centers(i).nFiring)])
%     xlabel('test')
    axis off
end
