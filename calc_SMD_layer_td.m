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
% calculate a diagonal S-matrix of a homogeneous layer for a 2D set of
% wavevector projections
%% input:
% xno, yno: numbers of Fourier harmonics in x and y dimensions
% kx0, ky0: zero order wavevector projections
% kgx, kgy: wavevector steps in x and y dimensions
% kh: layer thickness multiplied by the vacuum wavenumber
% eps: layer permittivity
%% output:
% SMD: diagonal interface S-matrix of size (2*no,2,2), where no = xno*yno
% block SMD(:,1,1) corresponds to refelection from below to below
% block SMD(:,2,2) corresponds to refelection from above to above
% block SMD(:,2,1) corresponds to transmission from below to above
% block SMD(:,1,2) corresponds to transmission from above to below
% central harmonic index is ind_0 = (ceil(xno/2)-1)*yno+ceil(yno/2)
% first (no) components of the S-matrix correspond to the TE polarization,
%  and indeces from (no+1) to (2*no) correspond to the TM polarization
%% implementation:
function [SMD] = calc_SMD_layer_td(xno, yno, kx0, ky0, kgx, kgy, kh, eps)
	no = xno*yno;
	SMD = zeros(2*no,2,2);
		% wavevector projections
	kx = kx0 + kgx*(linspace(1,xno,xno) - ceil(xno/2));
	ky = ky0 + kgy*(linspace(1,yno,yno) - ceil(yno/2));
	[kkx,kky] = meshgrid(kx,ky);
	kkx = reshape(kkx,1,[]);
	kky = reshape(kky,1,[]);
	kkxy = kkx.^2 + kky.^2;
	
	kz = sqrt(eps - kkxy);
	ind = angle(kz) < -1e-12;
	kz(ind) = -kz(ind);
	
	SMD(:,2,1) = exp((1i*kh)*kz);
	SMD(:,1,2) = SMD(:,2,1);
end
%
% end of calc_SMD_layer_td
%