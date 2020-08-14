%{
Copyright © 2020 Alexey A. Shcherbakov. All rights reserved.

This file is part of GratingFMM.

GratingFMM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

GratingFMM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GratingFMM. If not, see <https://www.gnu.org/licenses/>.
%}
%% description:
% return a 2D Toeplitz matrix (2D analog of the matlab function 'toeplitz')
%% input:
% M: vector of size (2*ny-1)*(2*nx-1)
% nx, ny: block dimensions of Toeplitz matrix
%% output:
% T: 2D block-Toeplitz matrix of size nx*ny
%% implementation:
function [T] = toeplitz2(M, nx, ny)
	CT = cell(1,2*nx-1);
	ind = toeplitz(linspace(nx,2*nx-1,nx),flip(linspace(1,nx,nx)));		
	for i = 1:2*nx-1
		CT{1,i} = toeplitz(M(ny:2*ny-1,i),M(ny:-1:1,i));
	end
	T = cell2mat(CT(ind));
end
%
% end of function toeplitz2
%