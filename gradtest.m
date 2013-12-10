function [ep,em] = gradtest(fh,m,dm)
% finite-difference gradient test based on Taylor expansion:
%   
%   f(m + alpha*dm) = f(m) + alpha*df(m)'*dm + O(alpha^2)
%   
% The error e(alpha) = f(m + alpha*dm) - f(m) - alpha*df(m)'*dm
% should be O(alpha^2) which can be verified by computing it for various
% alpha.
%
% use:
%   [ep,em] = gradtest(fh,m,dm)
% input:
%   fh - function handle [f,g] =fh(m)
%   m  - reference model
%   dm - perturbation
%
% output
%   ep - error based on forward difference
%   em - error based on backward difference
%
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

alpha = [1 .1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];
[f,g] = fh(m);
for k = 1:length(alpha)
    fp = fh(m + alpha(k)*dm);
    fm = fh(m - alpha(k)*dm);
    
    ep(k) = f+alpha(k)*dm'*g - fp;
    em(k) = f-alpha(k)*dm'*g - fm; 
end
