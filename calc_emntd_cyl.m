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
% calculate a permittivity Fourier matrix of a 2D grating with cylindrical
% pithches or holes being periodic in x and y dimensions of the 3D
% Cartesian coordinates
%% input:
% xno, yno: numbers of Fourier harmonics
% rpx, rpy: radius-to-period ratios
% eps_c: cylinder permittivity
% eps_m: permittivity of a medium which surrounds the cylinder
%% output:
% FE: cell array containing two Fourier matrices of the permittivity and
% inverse permittivity
%% implementation:
function FE = calc_emntd_cyl(xno, yno, rpx, rpy, eps_c, eps_m)
	FE = cellmat(1,2,2*yno-1,2*xno-1);

	ix = linspace(1,xno-1,xno-1);
	iy = linspace(1,yno-1,yno-1);
	[IX,IY] = meshgrid(ix,iy);
	
	fx = rpy*besselj(1,(2*pi*rpx)*ix) ./ ix;
	fy = rpx*besselj(1,(2*pi*rpy)*iy) ./ iy;
	FXY = sqrt((rpx*IX).^2 + (rpy*IY).^2);
	FXY = (rpx*rpy)*besselj(1,(2*pi)*FXY) ./ FXY;

	M = zeros(2*yno-1,2*xno-1);

	M(yno,xno) = pi*rpx*rpy;

	M(yno+1:2*yno-1,xno) = fy;
	M(yno-1:-1:1,xno) = M(yno+1:2*yno-1,xno);
	M(yno,xno+1:2*xno-1) = fx;
	M(yno,xno-1:-1:1) = M(yno,xno+1:2*xno-1);

	M(yno+1:2*yno-1,xno+1:2*xno-1) = FXY;
	M(yno+1:2*yno-1,xno-1:-1:1) = FXY;
	M(yno-1:-1:1,xno+1:2*xno-1) = FXY;
	M(yno-1:-1:1,xno-1:-1:1) = FXY;

	FE{1,1} = FE{1,1} + (eps_c - eps_m)*M;
	FE{1,2} = FE{1,2} + (1/eps_c - 1/eps_m)*M;
	FE{1,1}(yno,xno) = FE{1,1}(yno,xno) + eps_m;
	FE{1,2}(yno,xno) = FE{1,2}(yno,xno) + 1/eps_m;
end
%
% end of calc_emntd_cyl
%