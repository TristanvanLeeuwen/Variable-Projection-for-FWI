function [x,info] = GQN(fh,x0,options)
% Simple Quasi/Gauss-Newton method for penalty method
%
% use:
%   [xn,info] = GQN(fh,x0,lambda,options)
%
% input:
%   fh - function handle to misfit of the form [f1,f2,g] = fh(x)
%        where f = f1 + lambda*f2 is the function value, g is the gradient of the same size
%        as the input vector x. 
%   x0 - initial guess
%
%   options.maxiter - max iterations [default 10]
%   options.tol     - tolerance on 2-norm of gradient [1e-6]
%   options.M       - history size [5]
%   options.fid     - file id for output [1]
%   options.method  - {'LBFGS','GN'}
%   options.cgtol   - CG tolerance for GN subproblems
%   options.cgprec  - precondition CG with LBFGS Hessian
%
% output:
%   xn   - final estimate
%   info - iteration history (# iter, # eval, stepsize, f , ||g(x)||_2, cgiter)
%
% -------------------------------------------------------------------------
%      Copyright (C) 2013 Tristan van Leeuwen
%                         Centrum Wiskunde & Informatica
%                         Tristan.van.Leeuwen@cwi.nl
%  
%      This program is free software: you can redistribute it and/or modify
%      it under the terms of the GNU General Public License as published by
%      the Free Software Foundation, either version 3 of the License, or
%      (at your option) any later version.
%  
%      This program is distributed in the hope that it will be useful,
%      but WITHOUT ANY WARRANTY; without even the implied warranty of
%      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%      GNU General Public License for more details.
%  
%      You should have received a copy of the GNU General Public License
%      along with this program.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

if nargin<3
    options = [];
end

maxiter = getoption(options,'maxiter',10);
M       = getoption(options,'M',5);
tol     = getoption(options,'tol',1e-6);
fid     = getoption(options,'fid',1);
method  = getoption(options,'method','LBFGS');
cgtol   = getoption(options,'cgtol',1e-2);
cgmax   = getoption(options,'cgmax',100);
cgprec  = getoption(options,'cgprec',0);
write   = getoption(options,'write',0);

% initialization
fprintf(fid,'# iter, # eval, stepsize, f          , ||g(x)||_2, cgiter\n');
outputstr = '%6d, %6d, %1.2e, %1.5e, %1.5e, %6d\n';
n         = length(x0);
converged = 0;
iter      = 0;
x         = x0;
alpha0    = 1;
cgiter    = 0;
S         = zeros(n,M);
Y         = zeros(n,M);

% initial evaluation

[f,g,H] = fh(x);
tol     = tol*norm(g);
nfeval  = 1;

% output
info(1,:) = [iter,nfeval,alpha0,f,norm(g),cgiter];
fprintf(fid,outputstr,info(1,:));
if write
    dlmwrite(['x_' num2str(iter) '.dat'],x);
end

% main loop
while ~converged
    alpha0 = 1;
    % compute search direction
    switch method
        case 'LBFGS'% L-BFGS
            if iter == 0;
                alpha0 = 1/norm(g);
            end
            
            s = B(-g,S,Y,[]);
            
            p = -(s'*g)/(g'*g);
            if (p <= 0)
                fprintf(fid,'Loss of descent, reset history\n');
                S = zeros(n,M);
                Y = zeros(n,M);
                alpha0 = 1/norm(g);
            end
            
        case 'GN'
            
            if cgprec
                [s,flag,~,cgiter] = pcg(H,-g,cgtol,cgmax,@(x)B(x,S,Y,[]));
            else
                [s,flag,~,cgiter] = pcg(H,-g,cgtol,cgmax);
            end

            if flag
                fprintf(fid,'CG did not converge... (flag = %d)\n', flag);
            end
            
            p = -(s'*g)/(g'*g);
            if (p <= 0)
                fprintf(fid,'Loss of descent, resort to gradient-step\n');
                s      = -g;
                alpha0 = 1/norm(g);
            end
            
    end
    
    % linesearch
    [ft,gt,H,alpha,lsiter] = wWolfeLS(fh,x,f,g,s,alpha0);
    nfeval = nfeval + lsiter;
    
    % update
    xt = x + alpha*s;
    
    S = circshift(S,[0 -1]);
    Y = circshift(Y,[0 -1]);
    S = [S(:,1:end-1) xt - x];
    Y = [Y(:,1:end-1) gt - g];  

    f = ft;
    g = gt;
    x = xt;
    
    iter = iter + 1;
    
    % output
    info(iter + 1,:) = [iter,nfeval,alpha,f,norm(g),cgiter];
    fprintf(fid,outputstr,info(iter + 1,:));
    if write
        dlmwrite(['x_' num2str(iter) '.dat'],x);
    end
    
    % check convergence
    converged = (iter>maxiter)||(norm(g)<tol)||(alpha<1e-10)||(norm(s)<1e-10);
    
end

end

function [ft,gt,Ht,alpha,lsiter] = wWolfeLS(fh,x0,f0,g0,s0,alpha0)
% Simple Wolfe linesearch, adapted from
% (http://cs.nyu.edu/overton/mstheses/skajaa/msthesis.pdf, algorihtm 3).
%
%

lsiter = 0;
c1     = 1e-2;
c2     = 0.9;
done   = 0;
mu     = 0;
nu     = inf;
alpha  = alpha0/2;
while ~done
    if nu < inf
        alpha = (nu + mu)/2;
    else
        alpha = 2*alpha;
    end
    
    if lsiter < 20
        xt         = x0 + alpha*s0;
        [ft,gt,Ht] = fh(xt);
        lsiter     = lsiter + 1;
    else
        alpha = 0;
        break;
    end
      
    if ft > f0 + c1*alpha*g0'*s0
        nu = alpha;
    elseif gt'*s0 < c2*g0'*s0
        mu = alpha;
    else
        done = 1;
    end
end
end

function z = B(x,S,Y,H0)
% apply lbfgs inverse Hessian to vector
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
% use:
%   z = B(x,S,Y,b0,Binit)
%
% input:
%   x - vector of length n
%   S - history of steps in n x M matrix
%   Y - history of gradient differences in n x M matrix
%
% output
%   z - vector of length n
%

J = find(sum(abs(S),1));
S = S(:,J);
Y = Y(:,J);
M = size(S,2);
n = length(x);

if isempty(H0)&&(M>0)
    H0 = norm(Y(:,end))^2/(S(:,end)'*Y(:,end))*ones(n,1);
else
    H0 = ones(n,1);
end

alpha = zeros(M,1);
rho   = zeros(M,1);
for k = 1:M
    rho(k) = 1/(Y(:,k)'*S(:,k));
end
q = x;
% first recursion
for k = M:-1:1
    alpha(k) = rho(k)*S(:,k)'*q;
    q        = q - alpha(k)*Y(:,k);
end

% apply `initial' Hessian
z = q./H0;
% second recursion
for k = 1:M
    beta = rho(k)*Y(:,k)'*z;
    z    = z + (alpha(k) - beta)*S(:,k);
end
end




