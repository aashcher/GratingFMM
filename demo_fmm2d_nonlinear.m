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
%%
clc;
clear;
%% initialization
wl = 1.2; % wavelength in micrometers
wv = 2*pi/wl; % wavevector

gpx = 0.8; % grating period along x
gpy = 0.8; % grating period along y
gh = 0.4; % grating depth
kgx = wl/gpx;
kgy = wl/gpy;
kh = wv*gh;

xno = 15; % number of Fourier harmonics along x
yno = 15; % number of Fourier harmonics along y
no = xno*yno;
ixy = (ceil(xno/2)-1)*yno+ceil(yno/2);

eps_substrate = 1.44^2;
eps_air = 1;

	%% checkerboard grating
cx = [-0.25, 0.25]; % pixel centers relative to the period along x
cy = [-0.25, 0.25]; % pixel centers relative to the period along y
dx = [0.5, 0.5]; % pixel widths relative to the period along x
dy = [0.5, 0.5]; % pixel widths relative to the period along y
eps_pix = [3.45^2, 1, 3.45^2, 1]; % pixel permittivities
nonlinear_mask = [1, 0, 1, 0]; % pixel mask

	% incident plane wave
 % angles of incidence:
phi = 0;
theta = 0.001;
 % incident wavevector projections:
kx0 = sin(theta*pi/180)*cos(phi*pi/180);
ky0 = sin(theta*pi/180)*sin(phi*pi/180);
 % incident amplitude vectors:
Vinc_TE = zeros(2*no,2);
Vinc_TM = zeros(2*no,2);
Vinc_TE(ixy,2) = 1; % (0,0)-th harmonic - plane wave from the top
Vinc_TM(ixy+no,2) = 1; % (0,0)-th harmonic - plane wave from the top

	%% calculate nonlinear power dependence:
reflected_power = [];
nonlinear_coefficient = 1e-1; % use the coefficient to adapt demensionality
input_power = [0.1:0.5:100.1,99.6:-0.5:0.1]; % increasing and decreasing power
current_eps = eps_pix; % initial grating permittivity
for i_pnl = 1:numel(input_power)
	power = input_power(i_pnl); % nonlinear power coefficient
	[Vdif, Veff, current_eps] = fmmtd_nonlinear(xno, yno, kx0, ky0, kgx, kgy, kh, ...
																eps_substrate, eps_air, Vinc_TE, ...
																cx, cy, dx, dy, current_eps, eps_pix, ...
																(power*nonlinear_coefficient)*nonlinear_mask, 0.1, 1e-5);
	reflected_power = [reflected_power; [power*nonlinear_coefficient, Veff(ixy,2)]];
end

plot(reflected_power(:,1), reflected_power(:,2));
xlabel('power, arb. units');
ylabel('transmission');

return;
