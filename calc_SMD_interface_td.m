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
% calculate a diagonal S-matrix of an interface between two homogeneous
% isotropic media for two polarizations and a 2D set of wavevector
% projections
%% input:
% xno, yno: numbers of Fourier harmonics in x and y dimensions
% kx0, ky0: zero order wavevector projections
% kgx, kgy: wavevector steps in x and y dimensions
% eps1, eps2: permittivities of media below and above the interface
%% output:
% SMD: diagonal interface S-matrix of size (2*no,2,2), where no = xno*yno
% block SMD(:,1,1) corresponds to refelection from substrate to substrate
% block SMD(:,2,2) corresponds to refelection from superstrate to superstrate
% block SMD(:,2,1) corresponds to transmission from substrate to superstrate
% block SMD(:,1,2) corresponds to transmission from superstrate to substrate
% central harmonic index is ind_0 = (ceil(xno/2)-1)*yno+ceil(yno/2)
% first no components in each of the two first dimensions if the S-matrix
%  correspond to the TE polarization, and indeces from no+1 to 2*no 
%  correspond to the TM polarization
%% implementation:
function [SM] = calc_SMD_interface_td(xno, yno, kx0, ky0, kgx, kgy, eps1, eps2)
	no = xno*yno;
	ind_e = 1:no;
	ind_h = no+1:2*no;
		% propagation constants:
	[kz1, kz2] = fmmtd_kxyz(xno, yno, kx0, ky0, kgx, kgy, eps1, eps2);

	SM = zeros(2*no,2,2);
		% TE:
	SM(ind_e,1,1) = (kz1-kz2)./(kz1+kz2);
	SM(ind_e,2,1) = 1 + SM(ind_e,1,1);
	SM(ind_e,2,2) = -SM(ind_e,1,1);
	SM(ind_e,1,2) = 1 + SM(ind_e,2,2);
		% TM:
	SM(ind_h,1,1) = (eps2*kz1-eps1*kz2)./(eps2*kz1+eps1*kz2);
	SM(ind_h,2,1) = 1 + SM(ind_h,1,1);
	SM(ind_h,2,2) = -SM(ind_h,1,1);
	SM(ind_h,1,2) = 1 + SM(ind_h,2,2);
end
%
% end of calc_SMD_interface_td
%