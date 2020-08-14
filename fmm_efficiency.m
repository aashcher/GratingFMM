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
% collinear diffraction by 1D gratings
%% input:
% no: number of Fourier harmonics
% V_inc: incident field amplitude matrix of size (no,2)
% V_dif: diffracted field amplitude matrix of size (no,2)
% kx0: incident plane wave wavevector x-projection (Bloch wavevector)
% kg: wavelength-to-period ratio (grating vector)
% eps1, eps2: substrate and superstrate permittivities
% pol: polarization (either "TE" or "TM")
%% output:
% V_eff: efficiency matrix of size (no,2) if the if the incident field has
%  propagating harmonics, otherwise (if the incident field is purely evanescent)
%  the matrix of partial powers carried by each diffraction order
% first index of V_inc, V_dif, V_eff indicates diffraction harmonics
%  (0-th order index is ind_0 = ceil(no/2))
% second index of V_inc, V_dif, V_eff indicates whether the diffraction orders
%  are in the substrate (V(:,1)) or in the superstrate (V(:,2))
%% implementation
function [V_eff] = fmm_efficiency(no, V_inc, V_dif, kx0, kg, eps1, eps2, pol)
	[kz1, kz2] = fmm_kxz(no, kx0, 0, kg, eps1, eps2);
	kz1 = transpose(kz1);
	kz2 = transpose(kz2);
	if (strcmp(pol,'TM'))
		kz1 = kz1/eps1;
		kz2 = kz2/eps2;
	end

	V_eff = zeros(no,2);
	P_inc = sum( abs(V_inc(:,1).^2).*real(kz1) + abs(V_inc(:,2).^2).*real(kz2) );
	V_eff(:,1) = abs(V_dif(:,1).^2).*real(kz1);
	V_eff(:,2) = abs(V_dif(:,2).^2).*real(kz2);

	if (abs(P_inc) > 1e-15)
		V_eff = V_eff/P_inc;
	end
end
%
% END
%