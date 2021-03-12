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
%% demonstration script for the collinear 1D grating Fourier Modal Method calculations
clc;
format long;
%% initialization
wl = 1; % wavelength in micrometers
wv = 2*pi/wl; % wavevector
pol = 'TE'; % polarization, "TE" or "TM"
  % grating parameters
gp = 1.5; % grating period
gh = 0.4; % grating depth
ff = 0.4; % filling factor
	% permittivities
eps_sub = 1.5^2; % substrate permittivity
eps_gr = 2.1^2; %grating permittivity
eps_sup = 1; % superstrate permittivity
	% method parameters
no = 11; % number of Fourier modes
ind0 = ceil(no/2); % index of the zero harmonic (0th order diffraction)
	% incidence
theta = 10; % angle of incidence
kx0 = sin(theta*pi/180); % incidence wavevector projection
V_inc = zeros(no,2); % matrix of incident field amplitudes
	% second index indicates wether the amplitudes are in the substrate (1)
	% or in the superstrate (2)
V_inc(ind0,2) = 1; % plane wave coming from the superstrate
%V_inc(ind0,1) = 1; % plane wave coming from the substrate

%% scattering matrix calculation
	% calculate Fourier image matrix of the dielectric permittivity function
	% for a lamellar grating
FM = calc_emn_lam(no, ff, eps_gr, eps_sup); % lamellar grating
	% scattering matrix of the grating of size (no,no,2,2)
	% block SM(:,:,1,1) corresponds to refelection from substrate to substrate
	% block SM(:,:,2,2) corresponds to refelection from superstrate to superstrate
	% block SM(:,:,2,1) corresponds to transmission from substrate to superstrate
	% block SM(:,:,1,2) corresponds to transmission from superstrate to substrate
SM = fmm(no, kx0, wl/gp, wv*gh, eps_sub, eps_sup, FM, pol);

%% diffraction of a plane wave example
V_dif = zeros(no,2); % allocate a vector of diffracted field amplitudes
	% apply the calculated scattering matrix to the incident vector:
V_dif(:,1) = SM(:,:,1,1)*V_inc(:,1) + SM(:,:,1,2)*V_inc(:,2); % diffraction to the substrate
V_dif(:,2) = SM(:,:,2,1)*V_inc(:,1) + SM(:,:,2,2)*V_inc(:,2); % diffraction to the superstrate
  % check the power conservation:
b = fmm_balance(no, V_inc, V_dif, kx0, wl/gp, eps_sub, eps_sup, pol);
fprintf("power balance: %f\n",b); % precision of the power conservation
	% calculate the vector of diffraction efficiencies
V_eff = fmm_efficiency(no, V_inc, V_dif, kx0, wl/gp, eps_sub, eps_sup, pol);
fprintf("diffraction efficiencies:\n");
for i = -2:2
	fprintf(" %d order in the substrate: %f\n", i, V_eff(ind0+i,1));
end
for i = -2:2
	fprintf(" %d order in the superstrate: %f\n", i, V_eff(ind0+i,2));
end

%%% END OF FILE
