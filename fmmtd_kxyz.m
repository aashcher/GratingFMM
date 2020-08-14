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
% calculate grating vectors for 2D grating (periodicity in x and y directions)
% and plane wave propagation constants in homogeneous media below and 
% above the grating
%% input:
% xno, yno: numbers of Fourier harmonics in X and Y dimensions
% kx0, ky0: wavevector projections of an incident plane wave
% kgx, kgy: wavelength-to-period ratios
% eps1: permittivity of a medium below the grating (substrate)
% eps2: permittivity of a medium above the grating (superstrate)
%% output:
% kz1: row of propagation constants in the substrate
% kz2: row of propagation constants in the superstrate
% kx: row of grating vectors in x-direction
% ky: row of grating vectors in y-direction
% kxy: row of grating vectors in xy plane
%% implementation:
function [kz1, kz2, kx, ky, kxy] = fmmtd_kxyz(xno, yno, kx0, ky0, kgx, kgy, eps1, eps2)
	[kx,ky] = meshgrid(kx0 + kgx*(linspace(1,xno,xno) - ceil(xno/2)), ...
											 ky0 + kgy*(linspace(1,yno,yno) - ceil(yno/2)));
	kx = reshape(kx,1,[]);
	ky = reshape(ky,1,[]);
	kxy = kx.^2 + ky.^2;
	
	kz1 = sqrt(eps1 - kxy);
	kz2 = sqrt(eps2 - kxy);
	ind = angle(kz1) < -1e-12;
	kz1(ind) = -kz1(ind);
	ind = angle(kz2) < -1e-12;
	kz2(ind) = -kz2(ind);

	kxy = sqrt(kxy);
end


