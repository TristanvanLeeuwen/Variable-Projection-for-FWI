function P = getP(n,Iz,Ix)
% sampling operator, selects gridpoints from 2D grid:
%
% use:
%   P = getP(n,Iz,Ix);
%
% input:
%   n     - [nz,nx] number of gridpoints in z and x direction
%   Iz,Ix - indices
%
% ouput:
%   P - sparse matrix
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

I1 = speye(n(1));
I2 = speye(n(2));
P  = kron(I2(Ix,:),I1(Iz,:));

