function [ val idx ] =  multiDimensionalMax(A)
% Extracts maximum from a multi-dimensional array and returns corresponding coordinates in idx
% Timothee Masquelier timothee.masquelier@alum.mit.edu , Oct 2006

sz = size(A);
nDim = length(sz);
idx = zeros(1,nDim);

modulo = zeros(1,nDim);
modulo(1) = 1;
for d=1:nDim-1
    modulo(d+1) = modulo(d)*sz(d);
end

[val i] = max(A(:));

rest = i;
for d=nDim:-1:1
    idx(d) = floor((rest-1)/modulo(d))+1;    
    rest = rest - (idx(d)-1)*modulo(d);
end

%                 mA{1} = A;
% 
% for d=1:nDim
%     [mA{d+1} mI{d+1}] = max(mA{d});
% end
% 
% val = mA{nDim+1};
% idx(nDim) = mI{nDim+1};
% 
% for d=nDim:-1:2
%     idx(d-1) = mI{d}(idx(d)); 
% end
% 
%             [mA mI] = max(A);%max(RESULT.(CLASSIFIER_PARAM.classifier{c}).rocArea);
%             [mmA mmI] = max(mA);
%             [mmmA mmmI] = max(mmA);
%             
%             optimalArea = mmmA
%             optimalParam = mmmI(1,1,1)
%             optimalSph =    mmI(1,1,mmmI)
%             optimalMode =    mI(1,mmI(mmmI),mmmI)
% 
