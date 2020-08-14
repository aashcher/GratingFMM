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
% isotropic media for a fixed polarization and a 1D set of wavevector
% projections
%% input:
% no: number of Fourier harmonics
% kx0: zero harmonic wavevector projection
% kg: wavevector step
% eps1, eps2: permittivities of media below and above the interface
% pol: polarization, either 'TE' or 'TM'
%% output:
% SMD: diagonal interface S-matrix of size (no,2,2)
% block SMD(:,1,1) corresponds to refelection from below to below
% block SMD(:,2,2) corresponds to refelection from above to above
% block SMD(:,2,1) corresponds to transmission from below to above
% block SMD(:,1,2) corresponds to transmission from above to below
% central harmonic index is ind_0 = ceil(no/2)
%% implementation:
function [SMD] = calc_SMD_interface(no, kx0, kg, eps1, eps2, pol)
	[kz1, kz2] = fmm_kxz(no, kx0, 0, kg, eps1, eps2);

	SMD = zeros(no,2,2);
	if strcmp(pol,'TE')
		SMD(:,1,1) = (kz1-kz2)./(kz1+kz2);
		SMD(:,2,1) = 1 + SMD(:,1,1);
		SMD(:,2,2) = -SMD(:,1,1);
		SMD(:,1,2) = 1 + SMD(:,2,2);
	elseif strcmp(pol,'TM')
		SMD(:,1,1) = (eps2*kz1-eps1*kz2)./(eps2*kz1+eps1*kz2);
		SMD(:,2,1) = 1 + SMD(:,1,1);
		SMD(:,2,2) = -SMD(:,1,1);
		SMD(:,1,2) = 1 + SMD(:,2,2);
	else
		error('function calc_SM_interface: unknown polarization');
	end
end
%
% end of calc_SMD_interface
%