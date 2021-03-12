%{
Copyright © 2021 Alexey A. Shcherbakov. All rights reserved.

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
% calculate the power balance (check the energy conservation law) in case
% of the diffraction by 2D gratings with generally nonorthogonal directions
% of periodicity
%% input:
% no1, no2: numbers of Fourier harmonics
%  total number of Fourier harmonics is no = xno*yno;
% V_inc: incident field amplitude matrix of size (2*no, 2)
% V_dif: diffracted field amplitude matrix of size (2*no, 2)
%  first index of V_inc, V_dif indicates diffraction harmonics
%   with indices 1:no being TE orders and no+1:2*no being TM orders
%   (0-th order index is ind_0 = (ceil(xno/2)-1)*yno+ceil(yno/2))
%  second index of V_inc, V_dif, V_eff indicates whether the diffraction orders
%   are in the substrate (V(:,1)) or in the superstrate (V(:,2))
% kx0, ky0: incident plane wave wavevector x and y projections (Bloch wavevector projections)
% kg1, kg2: wavelength-to-period ratios
% psi: angle between periodicity directions (first direction is supposed 
%  to be parallel to axis X)
% eps1, eps2: substrate and superstrate permittivities
%% output:
% if the incident field has propagating harmonics the function returns the
%  normalized difference between the incident and diffractied field total
%  power, otherwise (if the incident field is purely evanescent) it returns
%  the total power carried by propagating diffraction orders
%% implementation
function [b] = fmmtdno_balance(no1, no2, VI, VD, kx0, ky0, kg1, kg2, psi, eps1, eps2)
	no = no1*no2;
	ib1 = 1:no;
	ib2 = no+1:2*no;
	[kz1, kz2] = fmmtdno_kxyz(no1, no2, kx0, ky0, kg1, kg2, psi, eps1, eps2);
	kz1 = transpose(kz1);
	kz2 = transpose(kz2);

	Pi = sum( abs((VI(ib1,1)).^2).*real(kz1) + abs((VI(ib1,2)).^2).*real(kz2) ) ...
		+ sum( abs((VI(ib2,1)).^2).*real(kz1/eps1) + abs((VI(ib2,2)).^2).*real(kz2/eps2) );
	Pd = sum( abs((VD(ib1,1)).^2).*real(kz1) + abs((VD(ib1,2)).^2).*real(kz2) ) ...
		+ sum( abs((VD(ib2,1)).^2).*real(kz1/eps1) + abs((VD(ib2,2)).^2).*real(kz2/eps2) );
		
	if (abs(Pi) > 1e-15)
		b = abs(Pd/Pi-1);
	else
		b = 0.5*Pd;
	end
end


%%% END OF FILE %%%