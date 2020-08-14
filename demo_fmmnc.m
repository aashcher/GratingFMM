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
%% demonstration script for the non-collinear 1D grating Fourier Modal Method calculations
clc;
format long;
%% initialization
wl = 1; % wavelength in micrometers
wv = 2*pi/wl; % wavevector
	% grating parameters
gp = 1.5; % grating period
gh = 0.5; % grating depth
	% dimensionless parameters
kg = wl/gp;
kh = wv*gh;
	% permittivities
eps_sub = 1.5^2; % substrate permittivity
eps_gr = 3.17^2; % grating permittivity
eps_sup = 1; % superstrate permittivity
	% method parameters
no = 15; % number of Fourier modes
ind0 = ceil(no/2); % index of the zero harmonic (0th order diffraction)
	% incidence
theta = 10; % angle of incidence (angle between the direction of incidence
						% and the vertical axis being perpendicular to the grating plane)
phi = 0; % turn angle in the grating plane,
					% phi = 0 corresponds to the case of collinear diffraction
	% incidence wavevector projections:
kx0 = sin(theta*pi/180)*cos(phi*pi/180);
ky0 = sin(theta*pi/180)*sin(phi*pi/180);
V_inc = zeros(2*no,2); % matrix of incident field amplitudes
	% first index indicates different Fourier harmonics, first (no) correspond
	% to the TE polarization; second (no) correspond to the TM polarization
	% second index indicates wether the amplitudes are in the substrate (1) in the superstrate (2)
V_inc(1*no+ind0,2) = 1; % "TM" polarized plane wave (0-th harmonic) incoming from the superstrate

%% scattering matrix calculation
	% calculate Fourier image matrix of the dielectric permittivity function
	% for a lamellar grating with filling factor 0.4
FM = calc_emn_lam(no,0.4,eps_gr,eps_sup);
	% scattering matrix of the grating
SM = fmmnc(no,kx0,ky0,kg,kh,eps_sub,eps_sup,FM);

%% diffraction of a plane wave example
V_sca = zeros(2*no,2); % allocate a vector of diffracted field amplitudes
	% apply the calculated scattering matrix to the incident vector:
V_sca(:,1) = SM(:,:,1,1)*V_inc(:,1) + SM(:,:,1,2)*V_inc(:,2); % diffraction to the substrate
V_sca(:,2) = SM(:,:,2,1)*V_inc(:,1) + SM(:,:,2,2)*V_inc(:,2); % diffraction to the superstrate
  % check the power conservation
b = fmmnc_balance(no,V_inc,V_sca,kx0,ky0,kg,eps_sub,eps_sup);
disp(b); % precicision of the power conservation
	% calculate the vector of diffraction efficiencies:
V_eff = fmmnc_efficiency(no,V_inc,V_sca,kx0,ky0,kg,eps_sub,eps_sup);
disp(V_eff(0*no+ind0,1)); % zero order power transmission to TE
disp(V_eff(0*no+ind0,2)); % zero order power reflection to TE
disp(V_eff(1*no+ind0,1)); % zero order power transmission to TM
disp(V_eff(1*no+ind0,2)); % zero order power reflection to TM

