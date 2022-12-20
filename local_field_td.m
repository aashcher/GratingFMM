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
function [MEx, MEy, MEz, MHx, MHy, MHz] = local_field_td(xno, yno, Vinc, kx, ky, ...
																beta, Ms, Ev, Hv, kh, nx, ny, nz, FE)
	no = xno*yno;
	ib1 = 1:no;
	ib2 = no+1:2*no;
	kx = transpose(kx);
	ky = transpose(ky);
	beta = transpose(beta);

	ME = toeplitz2(FE{1},xno,yno);
	MU = toeplitz2(FE{2},xno,yno);

		% mode amplitudes:
	Vmod = reshape(Ms\reshape(Vinc,[],1), [2*no,2]);
		% allocate field matrices:
	MEx = zeros(ny, nx, nz);
	MEy = zeros(ny, nx, nz);
	MEz = zeros(ny, nx, nz);
	MHx = zeros(ny, nx, nz);
	MHy = zeros(ny, nx, nz);
	MHz = zeros(ny, nx, nz);

		% loop over slices
	for i = 1:nz
		bexp1 = exp((1i*((i-0.5)/nz)*kh)*beta);
		bexp2 = exp((1i*((nz-i+0.5)/nz)*kh)*beta);

		Vp = Vmod(:,1).*bexp1;
		Vm = Vmod(:,2).*bexp2;

		Ex = Ev(ib1,:)*(Vp + Vm);
		Ey = Ev(ib2,:)*(Vp + Vm);
		Hx = Hv(ib1,:)*(Vp + Vm);
		Hy = Hv(ib2,:)*(Vp + Vm);
%		Ez = ME \ ( (1./beta(ib1)) .* (kx.*(ME*(Ev(ib1,:)*Vp))+ky.*(ME*(Ev(ib2,:)*Vp))) ...
%							- (1./beta(ib1)) .* (kx.*(ME*(Ev(ib1,:)*Vm))+ky.*(ME*(Ev(ib2,:)*Vm))));
%		Ez = -ME \ ( (1./beta(ib1)) .* (kx.*(ME*Ex) + ky.*(ME*Ey)) );
%		Ez = -(1./beta(ib1)) .* ( ME \ (kx.*(ME*Ex) + ky.*(ME*Ey)) );

%		Ez = ( beta(ib2).*Ey + Hx ) ./ ky;
%		Hz = ( beta(ib2).*Hy - ME*Ex ) ./ ky;
		Ez = -ME \ (kx.*Hy - ky.*Hx);
		Hz = ky.*Ex - kx.*Ey;

		nf = xno*yno;
		MEx(1:yno,1:xno,i) = reshape(Ex,[yno,xno]); MEx(:,:,i) = nf*fftshift(ifft2(MEx(:,:,i)));
		MEy(1:yno,1:xno,i) = reshape(Ey,[yno,xno]); MEy(:,:,i) = nf*fftshift(ifft2(MEy(:,:,i)));
		MEz(1:yno,1:xno,i) = reshape(Ez,[yno,xno]); MEz(:,:,i) = nf*fftshift(ifft2(MEz(:,:,i)));
		MHx(1:yno,1:xno,i) = reshape(Hx,[yno,xno]); MHx(:,:,i) = nf*fftshift(ifft2(MHx(:,:,i)));
		MHy(1:yno,1:xno,i) = reshape(Hy,[yno,xno]); MHy(:,:,i) = nf*fftshift(ifft2(MHy(:,:,i)));
		MHz(1:yno,1:xno,i) = reshape(Hz,[yno,xno]); MHz(:,:,i) = nf*fftshift(ifft2(MHz(:,:,i)));
	end
end


