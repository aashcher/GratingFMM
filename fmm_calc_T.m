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
% T-matrix of an interface between a photonic crystal and a homogeneous medium
% the case of the collinear diffraction on 1D gratings
%% input:
% no: number of Fourier harmonics
% EV: matrix of the electric field eigenvectors for the photonic crystal
% HV: matrix of the magnetic field eigenvectors for the photonic crystal
% kz: row of plane wave propagation constants in the homogeneous medium
% eps: permittivity the homogeneous medium
% pol: polarization
%% output:
% T: T-matrix of size 4*no by 4*no
%% implementation:
function [T] = fmm_calc_T(no, EV, HV, kz, eps, pol)
	T = zeros(2*no,2*no);
	ib1 = [1:no];
	ib2 = [no+1:2*no];
	if strcmp(pol,'TE')
		ikz = transpose(0.5./kz);
		T(ib1,ib1) = -HV.*ikz;
		T(ib2,ib1) = 0.5*EV - T(ib1,ib1);
		T(ib1,ib1) = T(ib1,ib1) + 0.5*EV;
		T(ib1,ib2) = T(ib2,ib1);
		T(ib2,ib2) = T(ib1,ib1);
	elseif strcmp(pol,'TM')
		eikz = transpose((0.5*eps)./kz);
		T(ib1,ib1) = EV.*eikz;
		T(ib2,ib1) = 0.5*HV - T(ib1,ib1);
		T(ib1,ib1) = T(ib1,ib1) + 0.5*HV;
		T(ib1,ib2) = -T(ib2,ib1);
		T(ib2,ib2) = -T(ib1,ib1);
	end
end

