%{
Copyright © 2022 Alexey A. Shcherbakov. All rights reserved.

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
%% nonlinear Fourier Modal Method
function [Vdif, Veff, eps_pnl, ME] = fmmtd_nonlinear(xno, yno, kx0, ky0, ...
																kgx, kgy, kh, eps_s, eps_c, eps_m, Vinc, ...
																cx, cy, dx, dy, eps_p, eps_p0, ...
																nlcf_p, weight, tol)
%	eps_p : initial values for pixel permittivities
%	eps_p0 : linear values of pixel permittivities
% weight : coefficient to update permittivity
		% auxiliary indexing parameters
	no = xno*yno;
	ib1 = 1:no;
	ib2 = no+1:2*no;
	ixy = (ceil(xno/2)-1)*yno+ceil(yno/2);

		%% calculate zero order approximation:
		% initialize the vector of diffraction amplitudes:
	FE = fmmtd_feps_binary_m(xno, yno, cx, cy, dx, dy, eps_p, eps_m);
	[SM, beta, Ms, EV, HV] = fmmtd(xno, yno, kx0, ky0, kgx, kgy, kh, eps_s, eps_c, FE);
	Vdif = smatrix_diffract(SM, Vinc);
	Veff = fmmtd_efficiency(xno, yno, Vinc, Vdif, kx0, ky0, kgx, kgy, eps_s, eps_c);
	R2 = Veff(ixy,1);% + Veff(no+ixy,2); % refelected power carried by zeroth TE and TM orders
	R1 = 1.1;

	eps_pnl = eps_p; % current array of nonlinear pixel permittivities

		%% loop for calculation of the self-consistent nonlinear permittivity
	while abs(R2-R1) > tol % check convergence via the refelection coefficient
			% calculate local field:
		[~, ~, kx, ky] = fmmtd_kxyz(xno, yno, kx0, ky0, kgx, kgy, eps_s, eps_c);
		[MEx, MEy] = local_field_td(xno, yno, Vinc, kx, ky, beta, Ms, EV, HV, kh, 1, FE);
		ME = abs(MEx.^2) + abs(MEy.^2);
			% calculate new permittivities
		eps_pnl0 = eps_pnl;
		eps_pnl = eps_p0;

		ixm = 1:xno;
		iym = 1:yno;
		for ip = 1:numel(eps_p) % loop over all nonlinear pixels
			mask = 0*(ixm + iym');
			mask(abs((yno-iym+1)/yno-0.5-cy(ip)) <= 0.5*dy(ip), ...
						abs((xno-ixm+1)/xno-0.5-cx(ip)) <= 0.5*dx(ip)) = 1;
			MEm = abs(ME .* mask);
			eps_pnl(ip) = eps_pnl(ip) + nlcf_p(ip)*mean(MEm(MEm > 1e-10));
		end

		eps_pnl = eps_pnl0 + weight*(eps_pnl - eps_pnl0);
		
			%% calculate diffraction for the new epsilon
		FE = fmmtd_feps_binary_m(xno, yno, cx, cy, dx, dy, eps_pnl, eps_m);
		[SM, beta, Ms, EV, HV] = fmmtd(xno, yno, kx0, ky0, kgx, kgy, kh, eps_s, eps_c, FE);
		Vdif = smatrix_diffract(SM, Vinc);
		Veff = fmmtd_efficiency(xno, yno, Vinc, Vdif, kx0, ky0, kgx, kgy, eps_s, eps_c);

		R1 = R2;
		R2 = Veff(ixy,2);

		fprintf("dR = %.5f, R = %.5f, de = %5f\n", abs(R1 - R2), R1, max(eps_pnl-eps_p));
	end
end


