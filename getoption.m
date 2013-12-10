function s = getoption(options,field,default)
% parse option from struct
%
% use
%   s = getoption(options,field,default)
%
% input
%   options - struct with option fields
%   field   - specify which field to parse
%   default - default value in case field does not exist
%
% output
%   s - value
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

s = default;
if isfield(options,field)
    s = getfield(options,field);
end