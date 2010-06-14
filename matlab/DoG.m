function filter = DoG(size,sigma1,sigma2)
x = repmat([1:size],size,1);
y = x';
d2 = (x-size/2-.5).^2 + (y-size/2-.5).^2;

filter = 1/sqrt(2*pi) * ( 1/sigma1 * exp(-d2/2/(sigma1^2)) - 1/sigma2 * exp(-d2/2/(sigma2^2)) );

% sum of weight must be 0
filter = filter - mean(filter(:));
filter = filter / norm(filter(:));

% % tmp
% figure
% subplot(1,2,1)
% mesh(filter);
% 
% subplot(1,2,2)
% plot(filter(round(size/2+.5),:));
% pause
