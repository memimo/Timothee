function y = solveInhibitoryNetwork(input,w,y0,lambda)
% Timothee Masquelier timothee.masquelier@alum.mit.edu July 2006
% See Foldiak 1990
% Computes the stable state (after an initial transient) for a neural
% network with n units sigmoidal tranfert function and inhibitory feedback connections.
% Param:
% input = n x 1 inputs at each node (after substracting threshold)
% w = n x n symmetric matrix of feedback inhibitory weights
% y0 = n x 1 initial output values
% lambda = param for the sigmoid transfer function


if sum(sum(abs(w)))==0 % particular case: no-inhib, no ODE to solve
    y = sigmoid(input,lambda)';
    return
end

n = length(input);
Y0 = [ y0; w(:); input; n; lambda];

% main function handle
f = @f;


tspan = [0 20];
for i=1:5
    [T,Y] = ode45(f,tspan,Y0);
    shift = round(1/4*size(Y,1));
    stability = max( abs( (Y(end,1:n)-Y(end-shift,1:n)) / (max(Y(end,:))-min(Y(end,:))) ));
    if stability < 10^-3
%         % tmp
%         figure('MenuBar','none')
%         subplot(1,2,1)
%         plot(input,'+')
%         subplot(1,2,2)
%         plot(T,Y(:,1:n));
%         pause

        % return final values
        y = Y(end,1:n);
        return
    else
        tspan(2) = 2*tspan(2);
    end
end

warning('Network did not stabilize')
        figure
        plot(T,Y(:,1:n));
% return final values
y = Y(end,1:n);



function result = f(t,Y)
% Timothee Masquelier timothee.masquelier@alum.mit.edu July 2006
% See Foldiak 1990
% main equa diff function Y'=f(t,Y)
% Y in reality contains [y w input n  lambda], but only y is variable

% fixed parameters
lambda=Y(end);
n=Y(end-1);
w = reshape(Y(n+1:n+n^2),[n n]);
input = Y(n+n^2+1:2*n+n^2);

% retrieve real y
y = Y(1:n);

% compute y' for non-fixed param (see Foldiak 1990)
result = [ sigmoid(w*y+input,lambda)-y ; zeros(n^2+n+2,1)];


function y = sigmoid(x,lambda)
y = 1 ./ (1+exp(-lambda*x));
