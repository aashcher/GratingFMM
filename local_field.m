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
%% calculate field in the grating region
function [MFx, MFy, MFz] = local_field(Vinc, beta, Ms, Ev, Hv, pol, kh, nz, FE, ieps_mask)
	no = size(Vinc,1);
		% mode amplitudes
	Vmod = reshape(Ms\reshape(Vinc,[],1),[no,2]);
	MFx = zeros(no, nz);
	MFy = zeros(no, nz);
	if strcmp(pol, 'TE')
		for i = 1:nz
			bexp1 = exp((1i*((i-0.5)/nz)*kh)*transpose(beta));
			bexp2 = exp((1i*((nz-i+0.5)/nz)*kh)*transpose(beta));
			Vp = Vmod(:,1).*bexp1;
			Vm = Vmod(:,2).*bexp2;
			MFy(:,i) = fftshift(fft(Ev*(Vp + Vm))); % E_y
			MFx(:,i) = fftshift(fft(Hv*(Vp - Vm))); % H_x
		end
	elseif strcmp(pol, 'TM')
		MU = eye(no) / toeplitz(FE(no:2*no-1,2),FE(no:-1:1,2));
		for i = 1:nz
			bexp1 = exp((1i*((i-0.5)/nz)*kh)*transpose(beta));
			bexp2 = exp((1i*((nz-i+0.5)/nz)*kh)*transpose(beta));
			Vp = Vmod(:,1).*bexp1;
			Vm = Vmod(:,2).*bexp2;
			MFy(:,i) = fftshift(fft(Hv*(Vp - Vm))); % H_y
			MFx(:,i) = ieps_mask.*fftshift(fft(MU*Ev*(Vp + Vm))); % E_x
			%MFy(:,i) = fftshift(fft(Ev*(Vp + Vm))); % E_x
		end
	end
end


