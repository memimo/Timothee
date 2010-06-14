function y = solveNonSelectiveInhibCircuit(input,alpha)
% Assumes inhibition I is non-selective and proportional to the mean
% activity: 
% I = alpha * mean(y)
% and:
% y = max(0,input-I);

global Input
global Alpha

if alpha==0 % special case: no inhib
    y = input;
    return
end

Input = input;
Alpha = alpha;

% solve for inhibition
I = fsolve(@f,0,optimset('Display','off'));

y = max(0,input-I);


function y = f(x)
global Input
global Alpha
y = Alpha*mean(max(0,Input-x)) - x;
