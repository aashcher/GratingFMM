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
% calculation of a grating S-matrix by the Fourier Modal Method
%  in the case of the collinear diffraction by 1D gratings being periodic in x
%  dimension of the Cartesian coordinates
%% input:
% no: number of Fourier harmonics
% kx0: incident plane wave wavevector x-projection (Bloch wavevector)
% kg: wavelength-to-period ratio (grating vector)
% kh: grating depth multiplied by the vacuum wavenumber
% eps1: permittivity of the substrate
% eps2: permittivity of the superstrate
% FE: Fourier matrix of the grating profile
% pol: polarization, either "TE" or "TM"
%% output:
% SM: scattering matrix of size (no,no,2,2)
% block SM(:,:,1,1) corresponds to refelection from substrate to substrate
% block SM(:,:,2,2) corresponds to refelection from superstrate to superstrate
% block SM(:,:,2,1) corresponds to transmission from substrate to superstrate
% block SM(:,:,1,2) corresponds to transmission from superstrate to substrate
% central harmonic index is ind_0 = ceil(no/2)
% for example, an ampitude of the transmitted wave to i-th diffraction order
%  from the substrate to the superstrate under the plane wave illumination
%  with unit amplitude is SM(ind_0+i, ind_0, 2, 1)
%% implementation
function SM = fmm(no, kx0, kg, kh, eps1, eps2, FE, pol)
		% wavevector projections
	[kz1, kz2, kx] = fmm_kxz(no, kx0, 0, kg, eps1, eps2);
		% block indices
	ib1 = 1:no;
	ib2 = no+1:2*no;
	  % solve the eigenvalue problem:
	ME = toeplitz(FE(no:2*no-1,1),FE(no:-1:1,1)); % permittivity Toeplitz matrix	
	if (strcmp(pol,'TE')) % if TE polarization
		ME(1:no+1:end) = ME(1:no+1:end) - (kx.^2);
		[EV,MB] = eig(ME);
		beta = transpose(sqrt(diag(MB)));
			% check the branch of the square root
		ind = angle(beta) < -1e-7;
		beta(ind) = -beta(ind);
			% eigen vector of the magnetic field
		HV = -EV.*beta; % Hx
	else % if TM polarization
		MU = eye(no) / toeplitz(FE(no:2*no-1,2),FE(no:-1:1,2)); % inverce permittivity Toeplitz matrix
		ME = -(diag(kx) / ME).*kx;
		ME(1:no+1:end) = ME(1:no+1:end) + 1;
		[EV,MB] = eig(ME*MU);
		beta = transpose(sqrt(diag(MB)));
			% check the branch of the square root
		ind = angle(beta) < -1e-7;
		beta(ind) = -beta(ind);
			% eigen vector of the magnetic field
		HV = (MU*EV)./beta;
	end

	bexp = exp((1i*kh)*beta);
	  % apply the boundary conditions:
	TS = fmm_calc_T(no,EV,HV,kz1,eps1,pol); % susbtrate-grating T-matrix
	TC = fmm_calc_T(no,EV,HV,kz2,eps2,pol); % grating-cover T-matrix
		% combine T-matrices
	M1 = zeros(2*no,2*no);
	M2 = zeros(2*no,2*no);
	M1(ib1,ib1) = TS(ib2,ib1);
	M1(ib1,ib2) = TS(ib2,ib2).*bexp;
	M1(ib2,ib1) = TC(ib1,ib1).*bexp;
	M1(ib2,ib2) = TC(ib1,ib2);
	M2(ib1,ib1) = TS(ib1,ib1);
	M2(ib1,ib2) = TS(ib1,ib2).*bexp;
	M2(ib2,ib1) = TC(ib2,ib1).*bexp;
	M2(ib2,ib2) = TC(ib2,ib2);

	SM = M2S(M1/M2);
end




