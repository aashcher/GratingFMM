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
% the case of the non-collinear diffraction on 1D gratings
%% input:
% no: number of Fourier harmonics
% EV: matrix of the electric field eigenvectors for the photonic crystal
% HV: matrix of the magnetic field eigenvectors for the photonic crystal
% kx: row of grating vector projections in x direction
% ky0: Bloch wavevector y-projection
% kxy: row of grating vector projections in xy plane
% kz: row of plane wave propagation constants in the homogeneous medium
% eps: permittivity the homogeneous medium
%% output:
% T: T-matrix of size 4*no by 4*no
%% implementation:
function [T] = fmmnc_calc_T(no, EV, HV, kx, ky0, kxy, kz, eps)
	T = zeros(4*no,4*no);

		% pre-calculate combinations of wavevector projections
	ikxy = 0.5./kxy;
	ky_ikxy = ky0*ikxy;
	kx_ikxy = kx.*ikxy;
	ikxyz = ikxy./kz;
	ky_ikxyz = transpose(ky0*ikxyz);%
	kx_ikxyz = transpose(kx.*ikxyz);%
	eky_ikxyz = eps*ky_ikxyz;
	ekx_ikxyz = eps*kx_ikxyz;
	kx_ikxy = transpose(kx_ikxy);%
	ky_ikxy = transpose(ky_ikxy);%

		% block indices
	ib1 = 1:no;
	ib2 = no+1:2*no;
	ib3 = 2*no+1:3*no;
	ib4 = 3*no+1:4*no;

		% fill the T-matrix blocks:
	T(ib1,ib1) = EV(ib1,ib1).*ky_ikxy - EV(ib2,ib1).*kx_ikxy;
	T(ib1,ib3) = T(ib1,ib1);
	TM = HV(ib1,ib1).*kx_ikxyz + HV(ib2,ib1).*ky_ikxyz;
	T(ib1,ib1) = T(ib1,ib1) + TM;
	T(ib1,ib3) = T(ib1,ib3) - TM;
	T(ib3,ib3) = T(ib1,ib1);
	T(ib3,ib1) = T(ib1,ib3);

	T(ib1,ib2) = EV(ib1,ib2).*ky_ikxy - EV(ib2,ib2).*kx_ikxy;
	T(ib1,ib4) = T(ib1,ib2);
	TM = HV(ib1,ib2).*kx_ikxyz + HV(ib2,ib2).*ky_ikxyz;
	T(ib1,ib2) = T(ib1,ib2) + TM;
	T(ib1,ib4) = T(ib1,ib4) - TM;
	T(ib3,ib2) = T(ib1,ib4);
	T(ib3,ib4) = T(ib1,ib2);

	T(ib2,ib1) = -EV(ib1,ib1).*ekx_ikxyz - EV(ib2,ib1).*eky_ikxyz;
	T(ib2,ib3) = T(ib2,ib1);
	TM = HV(ib1,ib1).*ky_ikxy - HV(ib2,ib1).*kx_ikxy;
	T(ib2,ib1) = T(ib2,ib1) + TM;
	T(ib2,ib3) = T(ib2,ib3) - TM;
	T(ib4,ib1) = -T(ib2,ib3);
	T(ib4,ib3) = -T(ib2,ib1);

	T(ib2,ib2) = -EV(ib1,ib2).*ekx_ikxyz - EV(ib2,ib2).*eky_ikxyz;
	T(ib2,ib4) = T(ib2,ib2);
	TM = HV(ib1,ib2).*ky_ikxy - HV(ib2,ib2).*kx_ikxy;
	T(ib2,ib2) = T(ib2,ib2) + TM;
	T(ib2,ib4) = T(ib2,ib4) - TM;
	T(ib4,ib2) = -T(ib2,ib4);
	T(ib4,ib4) = -T(ib2,ib2);
end
