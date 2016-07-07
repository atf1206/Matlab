% Dynamic Optimization with uncertain parameters % shocks
% By Andrew Fritz

% Abel 2: shocks, no stochastic control
% watch for high sums!

a = [0.914 -0.016; 
    0.097  0.424];
b = [0.305 0.424; 
    -0.101 1.459];
c = [-59.437;  
    -184.766];

[x,u,sum] = abel2(a, b, c);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Print the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = u';                       %   The optimal control vector
u

x = x';                       %   The optimal state vector
x 

Criterion = sum               %   The value of the criterion function