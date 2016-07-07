% Modified Dynamic Optimization program
% Run using: DynOpt4.m

% Taken From:
% Title: Quadratic-Linear Tracking problem for Abel
% Program name: qlpabel.m

% Based on the Chapter 4 of Stochastic Control for Economic Models
% example by David Kendrick.
% GAUSS version by Hans Amman, modified to MATLAB by Huber Salas
% with subsequent changes by Miwa Hattori and David Kendrick
% to implement a deterministic,
% two-control version of the Abel (1975) model (Jan 2005)
% Computes the optimal cost-to-go, control and state vectors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x2,u2,sum, xold] = abel4(a, b, c, k, t, x2, u2, n, m, xold, sum);

k2 = k;

xtar = [387.9;  
       85.3];
utar = [110.4; 
       147.17];
w = [0.0625 0;  
         0  1];
wn = [6.25 0;  
       0 100];
f = [0 0;   
     0 0];
lambda = [1    0; 
          0 0.444];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The Riccati Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kold = wn;                    %   Boundary condition
pold = -wn*xtar*(1.0075)^t;   %   Boundary condition

kstore = zeros(n,n,t);      % storage for dynamic Riccati matrices
pstore = zeros(n,t);        % storage for dynamic Riccati vectors

u = zeros(m,t+1);
x = zeros(n,t+1);

kstore(:,:,t) = kold(:,:);
pstore(1:n,t) = pold;

k = t-1;
while k >= 1;
  utark = (1.0075^k).*utar;   %   Time dependent targets
  xtark = (1.0075^k).*xtar;
  wsmall = -w*xtark;
  lambdas = -lambda*utark;

  knew = a'*kold*a+w-(a'*kold*b+f)*inv(b'*kold*b+lambda')*(f'+b'*kold*a);
                            %   Computing the Riccati matrices
  pnew = -(a'*kold*b+f)*inv(b'*kold*b+lambda')*(b'*(kold*c+pold)+lambdas)+...
             a'*(kold*c+pold)+wsmall;
                            %   Computing the tracking equation

  kold = knew;                %   Setup next period
  pold = pnew;
  kstore(:,:,k) = knew(:,:);
  pstore(1:n,k) = pnew;
  k = k-1;
end;                        %   End of the Riccati loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The Forward loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = k2;

  %a = a * randn;  
    
  utark = (1.0075^k).*utar;
  xtark = (1.0075^k).*xtar;

  wsmall = -w*xtark;
  lambdas = -lambda*utark;

  kold(:,:) = kstore(:,:,k+1);
  pold = pstore(1:n,k+1);
  glarge = -inv(b'*kold*b+lambda')*(f'+b'*kold*a);
  gsmall = -inv(b'*kold*b+lambda')*(b'*(kold*c+pold)+lambdas);
  
  uopt = glarge*xold+gsmall;
  xnew = a*xold+b*uopt+c;
  sum = sum+0.5*(xold-xtark)'*w*(xold-xtark)+0.5*(uopt-utark)'*lambda*(uopt-utark);
  %sum
  x2(1:n,k+1) = xold;
  u2(1:m,k+1) = uopt;
  xold = xnew;
  k = k+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   The Last Period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if k == t;
x2(1:n,t+1) = xold;
utark = (1.0075^k).*utar;
xtark = (1.0075^k).*xtar;
sum = sum+0.5*(xold-xtark)'*wn*(xold-xtark);
end;
