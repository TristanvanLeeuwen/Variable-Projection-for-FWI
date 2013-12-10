function A = getA(f,m,h,n)
% 5-point discretization of the 2D Helmholtz operator with absorbing
% boundary conditions.
%
% use:
%   A = getA(f,m,h,n);
%
% input:
%   f - frequency [Hz]
%   m - model [s^2/km^2]
%   h - [dz,dx] gridspacing in z and x direction [m]
%   n - [nz,nx] number of gridpoints in z and x direction
%
% ouput:
%   A - sparse matrix
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

% angular frequency (note factor 1e-3 because m is in [s^2/km^2])
omega = 1e-3*2*pi*f;

% number of gridpoints
N     = prod(n);

% Stiffness matrices
D1 = spdiags(ones(n(1),1)*[1 -2 1]/h(1).^2,[-1:1],n(1),n(1)); D1(1,1:2) = [1 -1]/h(1); D1(end,end-1:end) = [-1 1]/h(1);
D2 = spdiags(ones(n(2),1)*[1 -2 1]/h(2).^2,[-1:1],n(2),n(2)); D2(1,1:2) = [1 -1]/h(2); D2(end,end-1:end) = [-1 1]/h(2);

S = kron(speye(n(2)),D1) + kron(D2,speye(n(1)));

% Mass matrix including ABC's
w = [0 ones(1,n(1)-2) 0];
w = w(:)*[0 ones(1,n(2)-2) 0];
w = w(:);

M = omega^2*spdiags(w.*m,0,N,N) + 1i*omega*spdiags((1-w).*sqrt(m),0,N,N);

% 
A = M + S;



