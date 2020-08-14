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
% calculate a permittivity Fourier matrix of a 2D binary grating
% being periodic in x and y dimensions of the 3D Cartesian coordinates
%% input:
% xno, yno: numbers of Fourier harmonics
% cx, cy: rows of centers of a 2D rectangular mesh filling the grating period along
%   x and y dimensions normalized by the period (each value should be between -0.5 and 0.5)
% dx, dy: rows of widths of a 2D rectangular mesh elements along
%   x and y dimensions normalized by the period (each value should be between 0 and 1)
% eps: row of permittivities for each mesh element (length(eps) should be
%   equal to length(cx)*length(cy))
%% output:
% FE: cell array containing two Fourier matrices of the permittivity and
% inverse permittivity
%% implementation:
function [FE] = calc_emntd_bin(xno, yno, cx, cy, dx, dy, eps)
	nx = length(cx);
	ny = length(cy);
	if (length(cx)~=length(dx)) || (length(cy)~=length(dy)) || (length(eps)~=(nx*ny))
		error("incorrect binary grating definition");
	end

	FE = cellmat(1,2,2*yno-1,2*xno-1);

	[CX,CY] = meshgrid(cx,cy);
	[DX,DY] = meshgrid(dx,dy);

	ix = linspace(1,xno-1,xno-1);
	iy = linspace(1,yno-1,yno-1);
	[IX,IY] = meshgrid(ix,iy);
	
	for ip = 1:nx*ny
		fx = (sin(ix*pi*DX(ip))./(pi*ix)).*exp((-2*pi*1i*CX(ip))*ix);
		fy = (sin(iy*pi*DY(ip))./(pi*iy)).*exp((-2*pi*1i*CY(ip))*iy);
		FX = (sin(IX*pi*DX(ip))./(pi*IX)).*exp((-2*pi*1i*CX(ip))*IX);
		FY = (sin(IY*pi*DY(ip))./(pi*IY)).*exp((-2*pi*1i*CY(ip))*IY);

		M = zeros(2*yno-1,2*xno-1);

		M(yno+1:2*yno-1,xno) = DX(ip)*fy;
		M(yno-1:-1:1,xno) = conj(M(yno+1:2*yno-1,xno));
		M(yno,xno+1:2*xno-1) = DY(ip)*fx;
		M(yno,xno-1:-1:1) = conj(M(yno,xno+1:2*xno-1));

		M(yno+1:2*yno-1,xno+1:2*xno-1) = FX.*FY;
		M(yno+1:2*yno-1,xno-1:-1:1) = conj(FX).*FY;
		M(yno-1:-1:1,xno+1:2*xno-1) = FX.*conj(FY);
		M(yno-1:-1:1,xno-1:-1:1) = conj(FX.*FY);

		FE{1,1} = FE{1,1} + eps(ip)*M;
		FE{1,2} = FE{1,2} + (1/eps(ip))*M;
		FE{1,1}(yno,xno) = FE{1,1}(yno,xno) + DX(ip)*DY(ip)*eps(ip);
		FE{1,2}(yno,xno) = FE{1,2}(yno,xno) + DX(ip)*DY(ip)/eps(ip);
	end
end
%
% end of calc_emntd_cyl
%