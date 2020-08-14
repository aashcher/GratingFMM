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
%% demonstration script for the 2D grating Fourier Modal Method calculations
clc;
%format long;
%% initialization
wl = 1; % wavelength in micrometers
phi = 0;
theta = 0.001; % for normal incidence take a small non-zero value
kx0 = sin(theta*pi/180)*cos(phi*pi/180); % incidence wavevector horizontal projection (dimensionless) (Bloch wavevector)
ky0 = sin(theta*pi/180)*sin(phi*pi/180); % incidence wavevector horizontal projection (dimensionless) (Bloch wavevector)
gpx = 0.72; % grating period in x dimension
gpy = 0.72; % grating period it y dimension
gh = 0.5; % grating depth

wv = 2*pi/wl; % wavevector

  % dimensionless variables
kgx = wl/gpx;
kgy = wl/gpy;
kh = wv*gh;

xno = 15; % number of Fourier modes
yno = 15; % number of Fourier modes
no = xno*yno;
ixy = (ceil(xno/2)-1)*yno+ceil(yno/2);
eps_sub = 1.5; % substrate permittivity
eps_gr = 3.17^2; % get_epsAu_Drude(wl); %grating permittivity
eps_sup = 1; % superstrate permittivity

	%% S-matrix calculation
	% calculate Fourier image matrix of the dielectric permittivity function
	% for a 2D lamellar grating with filling factors 0.5,0.5
FE = calc_emntd_lam(xno,yno,0.5,0.5,eps_gr,eps_sup);
	% scattering matrix of the grating
SM = fmmtd(xno,yno,kx0,ky0,kgx,kgy,kh,eps_sub,eps_sup,FE);

%% diffraction of a plane wave example
	% incident field amplitude vector
V_inc = zeros(2*no,2);
V_inc((ceil(xno/2)-1)*yno+ceil(yno/2),2) = 1; % TE polarized plane wave (0-th harmonic) coming from the superstrate
	% define diffracted amplitude vector
V_dif = zeros(2*no,2);
	% calculate diffraction vector
V_dif(:,1) = SM(:,:,1,1)*V_inc(:,1) + SM(:,:,1,2)*V_inc(:,2); % diffraction to the substrate
V_dif(:,2) = SM(:,:,2,1)*V_inc(:,1) + SM(:,:,2,2)*V_inc(:,2); % diffraction to the superstrate
V_eff = fmmtd_efficiency(xno,yno,V_inc,V_dif,kx0,ky0,kgx,kgy,eps_sub,eps_sup);
disp("efficiency:");
disp(V_eff(ixy,1)); % 0th order power transmission coefficient
disp(V_eff(ixy,2)); % 0th order power reflection coefficient
	% calculate the power balance
b = fmmtd_balance(xno,yno,V_inc,V_dif,kx0,ky0,kgx,kgy,eps_sub,eps_sup);
disp("balance:");
disp(b);
