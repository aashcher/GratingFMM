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
%% demonstration script for the 2D grating Fourier Modal Method calculations
clc;
%format long;
%% initialization
wl = 1; % wavelength
wv = 2*pi/wl; % wavevector
	% angles for the incident wavevector in spherical coordinates:
phi = 0;
theta = 10.001; % for normal incidence take a small non-zero value
	% in-plane projections for the incident wavevector:
kx0 = sin(theta*pi/180)*cos(phi*pi/180);
ky0 = sin(theta*pi/180)*sin(phi*pi/180);
gpx = 1.5; % grating period in x dimension
gpy = 1.5; % grating period it y dimension
gh = 0.6; % grating depth
ffx = 0.5;
ffy = 0.4;

  % dimensionless variables:
kgx = wl/gpx;
kgy = wl/gpy;
kh = wv*gh;
	% numbers of Fourier orders in the reciprocal space:
xno = 15;
yno = 15;
no = xno*yno; % total number of Fourier orders
ixy = (ceil(xno/2)-1)*yno+ceil(yno/2); % zero harmonic index
eps_sub = 1; % substrate permittivity
eps_gr = 2.1^2; % grating permittivity
eps_sup = 1.5^2; % superstrate permittivity

	%% S-matrix calculation
fprintf("ORTHOGONAL PERIOD GRATING:\n");

	% calculate Fourier image matrix of the dielectric permittivity function
	% for a 2D photonic crystal slab with rectangular period and rectangular
	% rods:
FE = calc_emntd_lam(xno, yno, ffx, ffy, eps_gr, eps_sup);
	% OR
	% calculate Fourier image matrix of the dielectric permittivity function
	% for a 2D photonic crystal slab with rectangular period and circular
	% rods:
%r = 0.3*gpx; % rod radius
%FE = calc_emntd_cyl(xno, yno, r/gpx, r/gpy, eps_gr, eps_sup);
	% OR
	% calculate Fourier image matrix of the dielectric permittivity function
	% for a 2D grating with a pixelated period, e.g. checkerboard structure:
	% conditions: sum(dx) == 1; 0.5*(dx(i)+dx(i+1)) == cx(i+1)-cx(i); and the same for Y
%dx = [0.5, 0.5]; % array of pixel widths divided x-period along X
%dy = [0.5, 0.5]; % array of pixel widths divided y-period alobg Y
%cx = [-0.25, 0.25]; % array of pixel center positions along X relative to x-period
%cy = [-0.25, 0.25]; % array of pixel center positions along Y relative to y-period
%eps = [eps_gr, eps_sup, eps_sup, eps_gr]; % array of pixel permittivities, checkerboard
%eps = [eps_gr, eps_sup, eps_sup, eps_sup]; % OR array corresponding to calc_emntd_lam with ffx = ffy = 0.5
%FE = calc_emntd_bin(xno, yno, cx, cy, dx, dy, eps);

	% scattering matrix of the grating
SM = fmmtd(xno, yno, kx0, ky0, kgx, kgy, kh, eps_sub, eps_sup, FE);

%% diffraction of a plane wave example
	% incident field amplitude vector
V_inc = zeros(2*no,2);
V_inc(0*no+ixy,2) = 1; % TE polarized plane wave (0-th harmonic) coming from the superstrate
	% V_inc(1*no+ixy,2) = 1; % OR TM polarized plane wave coming from the superstrate
	% V_inc(0*no+ixy,1) = 1; % OR TE polarized plane wave coming from the bubstrate
	% V_inc(1*no+ixy,1) = 1; % OR TM polarized plane wave coming from the bubstrate
	% OR one can define a plane wave amplitude vector for a given beam
	
	% define diffracted amplitude vector
V_dif = zeros(2*no,2);
	% calculate diffraction vector
V_dif(:,1) = SM(:,:,1,1)*V_inc(:,1) + SM(:,:,1,2)*V_inc(:,2); % diffraction to the substrate
V_dif(:,2) = SM(:,:,2,1)*V_inc(:,1) + SM(:,:,2,2)*V_inc(:,2); % diffraction to the superstrate
	% calculate the power balance and diffraction efficiencies:
b = fmmtd_balance(xno,yno,V_inc,V_dif,kx0,ky0,kgx,kgy,eps_sub,eps_sup);
V_eff = fmmtd_efficiency(xno, yno, V_inc, V_dif, kx0, ky0, kgx, kgy, eps_sub, eps_sup);

fprintf("power balance: %e\n",b);

fprintf("diffraction efficiencies:\n");
fprintf(" 0th TE order efficiency in substrate: %f\n", V_eff(0*no+ixy,1));
fprintf(" 0th TM order efficiency in substrate: %f\n", V_eff(1*no+ixy,1));
fprintf(" 0th TE order efficiency in superstrate: %f\n", V_eff(0*no+ixy,2));
fprintf(" 0th TM order efficiency in superstrate: %f\n", V_eff(1*no+ixy,2));

%% diffraction by hexagonal cylinder gratings
fprintf("\nNON-ORTHOGONAL PERIOD GRATING:\n");

no1 = 15;
no2 = 15;
no = no1*no2;

f_ind = @(i,j) (ceil(no1/2)-1+i)*no2+ceil(no2/2)+j;
ixy = f_ind(0,0);

rc = 0.5; % cylinder radius
gp1 = 1.5; % period 1 (along X)
gp2 = 1.5; % period 2 directed at an angle psi to X
psi = pi / 3; % hexagonal grating

kg1 = wl/gp1;
kg2 = wl/gp2;

FEh = calc_emntdno_cyl(no1, no2, rc/gp1, rc/gp2, psi, eps_gr, eps_sup); % Fourier matrices
SMh = fmmtdno(no1, no2, kx0, ky0, kg1, kg2, psi, wv*gh, eps_sub, eps_sup, FEh); % S-matrix
V_inc = zeros(2*no,2);
V_inc(0*no+ixy,2) = 1; % incident field amplitude vector
	% calculate diffraction vector
V_dif = zeros(2*no,2);
V_dif(:,1) = SMh(:,:,1,1)*V_inc(:,1) + SMh(:,:,1,2)*V_inc(:,2); % diffraction to the substrate
V_dif(:,2) = SMh(:,:,2,1)*V_inc(:,1) + SMh(:,:,2,2)*V_inc(:,2); % diffraction to the superstrate
	% calculate the power balance and diffraction efficiencies:
b = fmmtdno_balance(no1, no2, V_inc, V_dif, kx0, ky0, kg1, kg2, psi, eps_sub, eps_sup);
fprintf("power balance: %e\n",b);

V_eff = fmmtdno_efficiency(no1, no2, V_inc, V_dif, kx0, ky0, kg1, kg2, psi, eps_sub, eps_sup);
fprintf("diffraction efficiencies in the substrate:\n");
fprintf(" (0,0) order efficiency: %f\n", V_eff(f_ind(0,0),1));
fprintf(" (1,0) order efficiency: %f\n", V_eff(f_ind(1,0),1));
fprintf(" (-1,0) order efficiency: %f\n", V_eff(f_ind(-1,0),1));
fprintf(" (0,1) order efficiency: %f\n", V_eff(f_ind(0,1),1));
fprintf(" (0,-1) order efficiency: %f\n", V_eff(f_ind(0,-1),1));
fprintf(" (1,1) order efficiency: %f\n", V_eff(f_ind(1,1),1));
fprintf(" (-1,-1) order efficiency: %f\n", V_eff(f_ind(-1,-1),1));

%%% END OF FILE %%%




