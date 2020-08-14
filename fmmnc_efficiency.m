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
% calculate a matrix of diffraction efficiencies in case of the
% non-collinear diffraction by 1D gratings being periodic in x dimension
%% input:
% no: number of Fourier harmonics
% V_inc: incident field amplitude matrix of size (2*no,2)
% V_dif: diffracted field amplitude matrix of size (2*no,2)
% kx0: incident plane wave wavevector x-projection (Bloch wavevector)
% ky0: incident plane wave wavevector y-projection
% kg: wavelength-to-period ratio (grating vector)
% eps1, eps2: substrate and superstrate permittivities
%% output:
% V_eff: efficiency matrix of size (2*no,2) if the if the incident field has
%  propagating harmonics, otherwise (if the incident field is purely evanescent)
%  the matrix of partial powers carried by each diffraction order
% first index of V_inc, V_dif, V_eff indicates diffraction harmonics
%   with indices 1:no being TE orders and no+1:2*no being TM orders
%  (0-th order index is ind_0 = ceil(no/2))
% second index of V_inc, V_dif, V_eff indicates whether the diffraction orders
%  are in the substrate (V(:,1)) or in the superstrate (V(:,2))
%% implementation
function [V_eff] = fmmnc_efficiency(no, V_inc, V_dif, kx0, ky0, kg, eps1, eps2)
	[kz1, kz2] = fmm_kxz(no, kx0, ky0, kg, eps1, eps2);
	kz1 = transpose(kz1);
	kz2 = transpose(kz2);
		
	ib1 = 1:no;
	ib2 = no+1:2*no;

	V_eff = zeros(2*no,2);
	P_inc = sum( abs(V_inc(ib1,1).^2).*real(kz1) + abs(V_inc(ib1,2).^2).*real(kz2) ) ...
				+ sum( abs(V_inc(ib2,1).^2).*real(kz1/eps1) + abs(V_inc(ib2,2).^2).*real(kz2/eps2) );
	V_eff(ib1,1) = abs(V_dif(ib1,1).^2).*real(kz1);
	V_eff(ib1,2) = abs(V_dif(ib1,2).^2).*real(kz2);
	V_eff(ib2,1) = abs(V_dif(ib2,1).^2).*real(kz1/eps1);
	V_eff(ib2,2) = abs(V_dif(ib2,2).^2).*real(kz2/eps2);

	if abs(P_inc) > 1e-15
		V_eff = V_eff/P_inc;
	else
		V_eff = 0.5*V_eff;
	end
end
%
% END
%