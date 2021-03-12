%{
Copyright © 2021 Alexey A. Shcherbakov. All rights reserved.

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
% pithches or holes being periodic in two dimensions of the 3D space
%% input:
% no1, no2: numbers of Fourier harmonics
% rp1, rp2: radius-to-period ratios
% psi: angle between periodicity directions (first direction is supposed 
% to be parallel to axis X)
% eps_c: cylinder permittivity
% eps_m: permittivity of a medium which surrounds the cylinder
%% output:
% FE: cell array containing two Fourier matrices of the permittivity and
% inverse permittivity
%% implementation
function FE = calc_emntdno_cyl(no1, no2, rp1, rp2, psi, eps_c, eps_m)
	FE = cellmat(1,2,2*no2-1,2*no1-1);

	cD = 2*pi / abs(sin(psi));

	m1 = linspace(1,2*no1-1,2*no1-1) - no1;
	m2 = linspace(1,2*no2-1,2*no2-1) - no2;
	[M1,M2] = meshgrid(m1,m2);

	M = cD*sqrt((rp1^2)*(M1.^2) + (rp2^2)*(M2.^2) - (2*rp1*rp2*cos(psi))*(M1.*M2));
	ind = abs(M) > 1e-12;
	M(ind) = (cD*rp1*rp2) * besselj(1,M(ind)) ./ M(ind);
	M(~ind) = 0.5*cD*rp1*rp2;

	FE{1} = (eps_c - eps_m) * M;
	FE{2} = (1/eps_c - 1/eps_m) * M;
	FE{1}(no2,no1) = FE{1}(no2,no1) + eps_m;
	FE{2}(no2,no1) = FE{2}(no2,no1) + 1/eps_m;
end

%%% END OF FILE