function y = solveNonSelectiveMaxInhibCircuit(input,alpha)
% Assumes inhibition I is non-selective and proportional to the max
% activity: 
% I = alpha * max(y)
% and:
% y(j) = max(0, input - alpha*max_i!=j(y(i))  );

global Input
global Alpha

if alpha==0 % special case: no inhib
    y = input;
    return
end

Input = input;
Alpha = alpha;

% f(input)

% solve for y
y = fsolve(@f,input,optimset('Display','off'));


function z = f(y)
global Input
global Alpha
n = length(y);
maxData = repmat(y,n,1) .* (ones(n,n)-eye(n)); % this is to compute max over i!=j
% max(maxData,[],2)'
z = max(0, Input - Alpha * max(maxData,[],2)') - y;
