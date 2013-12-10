function [f,g,H] = misfit(m,Q,D,model,srcest,corr)
% least-squares misfit for FWI with source estimation
%
% use:
%   [f,g] = misfit(m,Q,D,model,srcest)
%
% input:
%   m        - model [s^2/km^2]
%   Q        - sources as 2D array
%   D        - data as 3D array (receivers x sources x frequencies)
%   model.h  - [dz,dx] gridspacing in z and x direction [m]
%   model.n  - [nz,nx] number of gridpoints in z and x direction
%   model.f  - frequencies
%   model.Is - source location indices
%   model.Ir - receiver location indices
%   srcest   - source estimation on/off
%   corr     - include correction term in reduced Hessian
%
% output:
%   f - misfit
%   g - gradient
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

% get sampling operators
Ps = getP(model.n,model.Is{:});
Pr = getP(model.n,model.Ir{:});

% initialize misfit and gradient
f = 0;
g = 0*m;

% loop over frequencies 
for k = 1:length(model.f)
	% get Helmholtz operator
	Ak = getA(model.f(k),m,model.h,model.n);
	Gk = getG(model.f(k),m,model.h,model.n);
    
	% solve forward wave-equation and sample data
	Uk = Ak\(Ps'*Q);
	Sk = Pr*Uk;
    
	% compute source-weight
	if srcest
	  ck = real(sum(conj(Sk).*D(:,:,k),1)./sum(conj(Sk).*Sk,1))';     
	end
	Ck = spdiags(ck,0,size(Sk,2),size(Sk,2));
    
	% solve adjoint wave-equation
	Vk = Ak'\(Pr'*(Sk*Ck - D(:,:,k)));

	% compute misfit and gradient
	f = f + .5*norm(Sk*Ck - D(:,:,k),'fro')^2;
	g = g - sum(conj(Uk*Ck).*(Gk'*Vk),2);
end
g = real(g);
H = @(x)Hfun(x,m,Q,D,model,corr,srcest);

function y = Hfun(x,m,Q,D,model,corr,srcest)

% get sampling operators
Ps = getP(model.n,model.Is{:});
Pr = getP(model.n,model.Ir{:});

y = 0*x;

% loop over frequencies 
for k = 1:length(model.f)
	% get Helmholtz operator
	Ak = getA(model.f(k),m,model.h,model.n);
    Gk = getG(model.f(k),m,model.h,model.n);
    
	% solve forward wave-equation and sample data
	Uk = Ak\(Ps'*Q);
    Sk = Pr*Uk;
    
    % source-weight
    if srcest
        ck = real(sum(conj(Sk).*D(:,:,k),1)./sum(conj(Sk).*Sk,1))';     
    end
    Ck = spdiags(ck,0,size(Sk,2),size(Sk,2));
    
    % dc/dm
    dck = real( bsxfun(@times,(conj(Uk).*(Gk'*(Ak'\(Pr'*(D(:,:,k) - 2*Sk*Ck))))),1./sum(abs(Sk).^2,1)) );
    
    % dS = w*J*x + (dc'*x)*S
    dSk = (Pr*(Ak\(Gk*spdiags(x,0,length(x),length(x))*Uk)))*Ck;
    if corr
        dSk = dSk + Sk*spdiags((dck'*x),0,size(Sk,2),size(Sk,2));
    end
    
    % g = w*J'*p + (s'*p)*dc
    g = sum(conj(Uk*Ck).*(Gk'*(Ak'\(Pr'*dSk))),2);
    if corr
        g = g + dck*transpose(sum(conj(Sk).*dSk,1));
    end
    
    y = y + real(g);
end