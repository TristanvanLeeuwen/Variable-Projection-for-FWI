function S = F(m,Q,model)
% Modeling operator based on 5-point discretization of Helmholtz operator
% 
%   S = PA(m)^{-1}Q
%
% use:
%   S = F(m,Q,model)
%
% input:
%   m        - model [s^2/km^2]
%   Q        - sources as 2D array (each column is a vector of source-weights, typically Q = I)
%   S        - data as 3D array (receivers x sources x frequencies)
%   model.h  - [dz,dx] gridspacing in z and x direction [m]
%   model.n  - [nz,nx] number of gridpoints in z and x direction
%   model.f  - frequencies
%   model.Is - source location indices Is{1} - z locations, Is{2} - x locations
%   model.Ir - receiver location indices
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

S = zeros(size(Pr,1),size(Ps,1),length(model.f));

% loop over frequencies 
for k = 1:length(model.f)
	% get Helmholtz operator
	Ak = getA(model.f(k),m,model.h,model.n);

	% solve forward wave-equation and sample data
	Uk       = Ak\Ps'*Q;
	S(:,:,k) = Pr*Uk;
end