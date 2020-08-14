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
%  in the case of the diffraction by 2D gratings being periodic in x and y
%  dimensions of the Cartesian coordinates
%% input:
% xno, yno: numbers of Fourier harmonics
% kx0, ky0: incident plane wave x and y projections (Bloch wavevector projections)
% kgx, kgy: wavelength-to-period ratios (grating vector projections)
% kh: grating depth multiplied by the vacuum wavenumber
% eps1: permittivity of the substrate
% eps2: permittivity of the superstrate
% FE: Fourier matrix of the grating profile
%% output:
% SM: scattering matrix of size (2*no,2*no,2,2) where no = xno*yno
%  is the total number of Fourier harmonics
% block SM(:,:,1,1) corresponds to refelection from substrate to substrate
% block SM(:,:,2,2) corresponds to refelection from superstrate to superstrate
% block SM(:,:,2,1) corresponds to transmission from substrate to superstrate
% block SM(:,:,1,2) corresponds to transmission from superstrate to substrate
% first (no) components in each of the two first dimensions of the S-matrix
%  correspond to the TE polarization, and indeces from (no+1) to (2*no) 
%  correspond to the TM polarization
% central harmonic index is ind_0 = (ceil(xno/2)-1)*yno+ceil(yno/2)
% for example, an ampitude of the TE-polarized transmitted wave to i-th diffraction order
%  from the substrate to the superstrate under the TM plane wave illumination
%  with unit amplitude is SM(ind_0+i, no+ind_0, 2, 1)
%% implementation
function [SM] = fmmtd(xno, yno, kx0, ky0, kgx, kgy, kh, eps1, eps2, FE)
	no = xno*yno;
	ib1 = 1:no;
	ib2 = no+1:2*no;
	ib3 = 2*no+1:3*no;
	ib4 = 3*no+1:4*no;

		% wavevector projections
	[kz1, kz2, kx, ky, kxy] = fmmtd_kxyz(xno, yno, kx0, ky0, kgx, kgy, eps1, eps2);
	
	ME = toeplitz2(FE{1,1},xno,yno);
	MU = toeplitz2(FE{1,2},xno,yno);

		% solve the eigenvalue problem:
	EV = zeros(2*no,2*no);
	HV = zeros(2*no,2*no);
		% matrix for the electric field
	EV(ib1,ib1) = ((kx').*MU).*ky;
	EV(ib1,ib2) = -((kx').*MU).*kx;
	EV(2*no*no+1:2*no+1:end) = EV(2*no*no+1:2*no+1:end) + 1;
	EV(ib2,ib1) = ((ky').*MU).*ky;
	EV(no+1:2*no+1:2*no*no) = EV(no+1:2*no+1:2*no*no) - 1;
	EV(ib2,ib2) = -((ky').*MU).*kx;
		% matrix for the magnetic field
	HV(2*no*no+no+1:2*no+1:end) = HV(2*no*no+no+1:2*no+1:end) + kx.*ky;
	HV(1:2*no+1:2*no*no) = -HV(2*no*no+no+1:2*no+1:end);
	HV(ib1,ib2) = -ME;
	HV(2*no*no+1:2*no+1:end) = HV(2*no*no+1:2*no+1:end) + kx.^2;
	HV(ib2,ib1) = ME;
	HV(no+1:2*no+1:2*no*no) = HV(no+1:2*no+1:2*no*no) - ky.^2;

		% solve the eigenvalue problem for the electric field
	[EV,MB] = eig(EV*HV);
	beta = transpose(sqrt(diag(MB))); % row of eigenvalues
	ind = angle(beta) < -1e-7; % check the branch of the square root
	beta(ind) = -beta(ind);

		% calculate the magnetic field amplitude eigen vectors
	HV = (HV*EV).*(1./beta);

		% calculate T-matrices
	TS = fmmtd_calc_T(no,EV,HV,kx,ky,kxy,kz1,eps1); % substrate-grating T-matrix
	TC = fmmtd_calc_T(no,EV,HV,kx,ky,kxy,kz2,eps2); % grating-cover T-matrix

		% combine T-matrices
	bexp = exp((1i*kh)*beta);
	M1 = zeros(4*no,4*no);
	M2 = zeros(4*no,4*no);

	M1([ib1,ib2],[ib1,ib2]) = TS([ib3,ib4],[ib1,ib2]);
	M1([ib3,ib4],[ib3,ib4]) = TC([ib1,ib2],[ib3,ib4]);
	M1(ib1,ib3) = TS(ib3,ib3).*bexp(ib1);
	M1(ib2,ib3) = TS(ib4,ib3).*bexp(ib1);
	M1(ib1,ib4) = TS(ib3,ib4).*bexp(ib2);
	M1(ib2,ib4) = TS(ib4,ib4).*bexp(ib2);
	M1(ib3,ib1) = TC(ib1,ib1).*bexp(ib1);
	M1(ib4,ib1) = TC(ib2,ib1).*bexp(ib1);
	M1(ib3,ib2) = TC(ib1,ib2).*bexp(ib2);
	M1(ib4,ib2) = TC(ib2,ib2).*bexp(ib2);

	M2([ib1,ib2],[ib1,ib2]) = TS([ib1,ib2],[ib1,ib2]);
	M2([ib3,ib4],[ib3,ib4]) = TC([ib3,ib4],[ib3,ib4]);
	M2(ib1,ib3) = TS(ib1,ib3).*bexp(ib1);
	M2(ib2,ib3) = TS(ib2,ib3).*bexp(ib1);
	M2(ib1,ib4) = TS(ib1,ib4).*bexp(ib2);
	M2(ib2,ib4) = TS(ib2,ib4).*bexp(ib2);
	M2(ib3,ib1) = TC(ib3,ib1).*bexp(ib1);
	M2(ib4,ib1) = TC(ib4,ib1).*bexp(ib1);
	M2(ib3,ib2) = TC(ib3,ib2).*bexp(ib2);
	M2(ib4,ib2) = TC(ib4,ib2).*bexp(ib2);

		% S-matrix
	SM = M2S(M1/M2);
end
%
% END
%