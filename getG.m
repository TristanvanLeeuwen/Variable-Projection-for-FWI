function G = getG(f,m,h,n)
% Jacobian of A*u, see getA.m.
% The Jacobian of A*u can be written as d(A*u)/dm = G*diag(m)*u
%
% use:
%   G = getG(f,m,h,n);
%
% input:
%   f - frequency [Hz]
%   m - model [s^2/km^2]
%   h - [dz,dx] gridspacing in z and x direction [m]
%   n - [nz,nx] number of gridpoints in z and x direction
%
% ouput:
%   G - sparse matrix
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
N     = prod(n);

% 
w = [0 ones(1,n(1)-2) 0];
if n(2)>1
    w = w(:)*[0 ones(1,n(2)-2) 0];
end
w = w(:);

% construct sparse matrix; we ignore the contributions from the boundary for the time being to ensure
% that all quantities are real.
G = omega^2*spdiags(w,0,N,N);   % + 1i*omega*spdiags((1-w)./(2*sqrt(m)),0,N,N);

