% Dynamic Optimization 
% (eventually with uncertain parameters & shocks)
% By Andrew Fritz

% Title: DynOpt8.m 
% Adds Kalman Filter (sorta)

% a, b, and c are the estimated parameters

t = 5; 
a = [0.914 -0.016; 
    0.097  0.424];
b = [0.305 0.424; 
    -0.101 1.459];
c = [-59.437;  
    -184.766];
x0 = [387.9;      
       85.3];

areal = a;  ya = 0;
breal = b;  yb = 0;
creal = c;  yc = 0;

sum = 0;
k = 0;
nn = 0;
xold = x0;

n = 2; m = 2;
u2 = zeros(m,t+1);
x2 = zeros(n,t+1);

while nn <= t-1;
areal = (.1 * randn) + areal; % actual parameters and shocks
breal = (.1 * randn) + breal; % variance = .01
creal = (.1 * randn) + creal;
ya = areal + randn; % measured parameter values
yb = breal + randn; % with high uncertainty
yc = creal + randn; % Variance = 1
K = .01 / (.01 + 1);    % optimum based on variances
a = ((.1 * randn) + a) + K * (ya - ((.1 * randn) + a)); % estimated parameters
b = ((.1 * randn) + b) + K * (yb - ((.1 * randn) + b));
c = ((.1 * randn) + c) + K * (yc - ((.1 * randn) + c));

[x2,u2,sum,xold] = abel8(a, b, c, k, t, x2, u2, n, m, xold, sum);

k = k + 1;
nn = nn + 1;
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Print the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u2 = u2';                       %   The optimal control vector
%u2

x2 = x2';                       %   The optimal state vector
%x2 

Criterion = sum / t              %   The value of the criterion function