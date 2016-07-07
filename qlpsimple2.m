% Title: Simple Quadratic-Linear Problem
% Program name: qlpsimple.m
% Based on qlpabel.m program by Amman, Hattori, Kendrick and Salas

%   Preliminaries
t = 5;   
a = 0.7;   b = -0.3;   
sum = 0;

u = zeros(1,t);   x = zeros(1,t);
%y = zeros(1,t);

x0 = -1;
xold = x0;
x2old = x0;
%y(1,0) = x0;

kg = .5;
%   The Forward loop
k = 0;
while k <= t;    
  noise1 = 5 * randn;
  noise2 = 5 * randn;
  glarge = - a / b;
  uopt = glarge * x0;
  
  xnew = a * x2old + b * uopt + noise1;
  ynew = xnew + noise2;
  xnew2 = xnew + kg * (ynew - xnew);
  
  sum = sum + 0.5 * xold^2;
  x(1,k+1) = xold;
  u(1,k+1) = uopt;
  y(1,k+1) = x2old;
  xold = xnew;
  x2old = xnew2;
  k = k+1;
end;                        

%   Print the solution
u = u'                        
x = x'
y = y'
Criterion = sum        
 

%graph
%t = 0 : dt : duration;
%t = t';
%plot(t,pos,'r',t,poshat,'g',t,posmeas,'b');
plot(t,x','r',t,y','g');
grid;
xlabel('Time (sec)');
ylabel('something');
title('Kalman Filter Performance');