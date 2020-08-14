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
% calculate a permittivity Fourier matrix of a 2D lamellar grating
% being periodic in x and y dimensions of the 3D Cartesian coordinates
%% input:
% xno, yno: numbers of Fourier harmonics
% ax, ay: aspect ratios, should be between 0 and 1
% eps1: pitch permittivity (occupies ax*ay volume fraction of the period)
% eps2: permittivity of a medium which surrounds the pitch
%   (occupies (1-ax)*(1-ay) volume fraction of the period)
%% output:
% FE: cell array containing two Fourier matrices of the permittivity and
% inverse permittivity
%% implementation:
function [FE] = calc_emntd_lam(xno, yno, ax, ay, eps1, eps2)
	FB = zeros(2*yno-1,2*xno-1);
	FE = cell(1,2);
	FE{1,1} = eps2 + FB;
	FE{1,2} = 1/eps2 + FB;

	ix = linspace(1,xno-1,xno-1);
	iy = linspace(1,yno-1,yno-1);
	[IX,IY] = meshgrid(ix,iy);
		
	fx = sin(ix*pi*ax)./(pi*ix);
	fy = sin(iy*pi*ay)./(pi*iy);
	FX = sin(IX*pi*ax)./(pi*IX);
	FY = sin(IY*pi*ay)./(pi*IY);

	FB(yno+1:2*yno-1,xno) = ax*fy;
	FB(yno-1:-1:1,xno) = FB(yno+1:2*yno-1,xno);

	FB(yno,xno+1:2*xno-1) = ay*fx;
	FB(yno,xno-1:-1:1) = FB(yno,xno+1:2*xno-1);

	FB(yno+1:2*yno-1,xno+1:2*xno-1) = FX.*FY;
	FB(yno+1:2*yno-1,xno-1:-1:1) = FB(yno+1:2*yno-1,xno+1:2*xno-1);
	FB(yno-1:-1:1,xno+1:2*xno-1) = FB(yno+1:2*yno-1,xno+1:2*xno-1);
	FB(yno-1:-1:1,xno-1:-1:1) = FB(yno+1:2*yno-1,xno+1:2*xno-1);

	FE{1,1} = (eps1 - eps2)*FB;
	FE{1,2} = (1/eps1 - 1/eps2)*FB;
	FE{1,1}(yno,xno) = eps2 + (eps1 - eps2)*ax*ay;
	FE{1,2}(yno,xno) = 1/eps2 + (1/eps1 - 1/eps2)*ax*ay;
end
%
% end of calc_emntd_lam
%