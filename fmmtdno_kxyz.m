%%
%{
Copyright © 2021 Alexey A. Shcherbakov. All rights reserved.

This file is part of GratingFMM

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
% calculate grating vectors for 2D grating with nonorthogonal periods
% and plane wave propagation constants in homogeneous media below and 
% above the grating
%% input:
% no1, no1: numbers of Fourier harmonics
% kx0, ky0: X and Y wavevector projections of an incident plane wave
% kg1, kg2: wavelength-to-period ratios
% psi: angle between periodicity directions (first direction is supposed 
%  to be parallel to axis X)
% eps1: permittivity of a medium below the grating (substrate)
% eps2: permittivity of a medium above the grating (superstrate)
%% output:
% kz1: row of propagation constants in the substrate
% kz2: row of propagation constants in the superstrate
% kx: row of grating vectors in x-direction
% ky: row of grating vectors in y-direction
% kxy: row of grating vectors in xy plane
%% implementation:
function [kz1, kz2, kx, ky, kxy] = fmmtdno_kxyz(no1, no2, kx0, ky0, kg1, kg2, psi, eps1, eps2)
	[N1, N2] = meshgrid(linspace(1,no1,no1) - ceil(no1/2), ...
											linspace(1,no2,no2) - ceil(no2/2));
	td = 1/abs(sin(psi));
	kx = kx0 + td * ( kg1 * sin(psi) * N1 );
	ky = ky0 + td * ( -kg1 * cos(psi) * N1 + kg2 * N2 );
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


%%% END OF FILE