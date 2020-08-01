function [x,cost] = SR_ASD(y,lam0,lam1)
% 
% INPUT
%   y - raw data
%   sigma - noise level (sigma > 0)
%   lam0,lam1-parameters ( lam0,lam1> 0)
% OUTPUT
%   x - filtering signal
%   cost - cost function histoy
%initialization
EPS = 1E-10;% the minimum error 
a=1;
Nit=20; %Nit - number of iterations
%% the expression of regularization term
phi = @(x, a)  a*sqrt(x.^2 + EPS);
psi = @(x, a)  a*sqrt(x.^2 + EPS);
%%
cost = zeros(1, Nit);       % cost function history
N = length(y);
e = ones(N, 1);
D = spdiags([e,-2*e, e], [0 1 2], N-2, N);% 
  EE=spdiags(e,0, N, N); % identity matrix
b=y;
x = y;
for k = 1:Nit 
    % spdiags is used to improve the solving speed
    Lam0 = spdiags( lam0./psi(x, a), 0, N, N);
    Lam1 = spdiags( lam1./psi(D*x, a), 0, N-2, N-2);% 
    Q =EE+ (Lam0 + D'*Lam1*D); 
    u = Q \ b;
    x=u;
    cost(k) = lam0 * sum(phi(x, a)) + lam1 * sum(phi(D*x, a)) + 0.5 * sum(abs(x-y).^2);  %objection function
    if cost(k)-cos(k-1)< EPS
        break
    end
end

