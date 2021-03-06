% Dynamic Optimization with uncertain parameters % shocks
% By Andrew Fritz

% Abel 3: Functions exactly like qlpabel.m
% 

t = 7; 
a = [0.914 -0.016; 
    0.097  0.424];
b = [0.305 0.424; 
    -0.101 1.459];
c = [-59.437;  
    -184.766];

sum = 0;
k = 0;
nn = 0;
xold = x0;

n = 2; m = 2;
u2 = zeros(m,t+1);
x2 = zeros(n,t+1);

while nn <= t-1;
[x2,u2,sum,xold] = abel3(a, b, c, k, t, x2, u2, n, m, xold, sum);
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

Criterion = sum               %   The value of the criterion function