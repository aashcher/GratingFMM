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
% Fourier matrix of the permittivity and the inverse permittivity of a
% 1D lamellar grating
%% input:
% no - number of Fourier harmonics
% alp - lamellae width relative to the grating period (0 < alp < 1)
% eps1 - lamellae permittvity
% eps2 - surrounding medium permittvity
%% output:
% FE: Fourier matrix of size (2*no-1,2), which contains
% Fourier components of the permittivity (FE(:,1)) and the inverse
% permittivity (FE(:,2))
%% implementation
function [FE] = calc_emn_lam(no, alp, eps1, eps2)
	FE = zeros(2*no-1,2);
	te1 = eps1 - eps2;
	te2 = 1/eps1 - 1/eps2;
		% zero harmonic
	FE(no,1) = eps1*alp + eps2*(1-alp);
	FE(no,2) = alp/eps1 + (1-alp)/eps2;
		% non-zero harmonics:
	ind = linspace(1,no-1,no-1);
	ifun = transpose(sin(ind*pi*alp)./(ind*pi));
	FE(no+1:2*no-1,1) = te1*ifun;
	FE(no+1:2*no-1,2) = te2*ifun;
	FE(no-1:-1:1,1) = FE(no+1:2*no-1,1);
	FE(no-1:-1:1,2) = FE(no+1:2*no-1,2);
end
%
% end of calc_emn_lam
%