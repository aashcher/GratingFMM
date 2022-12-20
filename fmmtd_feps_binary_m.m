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
	% binary grating in a medium (pixel permittivitities being different from 
	% the medium permittivity are given only)
%% implementation:
function FE = fmmtd_feps_binary_m(xno, yno, cx, cy, dx, dy, eps, eps_m)
	np = numel(cx); % number of pixels

		% check for array compatibility
	if (numel(cx)~=numel(dx)) || (numel(cy)~=numel(dy)) || (numel(eps)~=np) ...
			|| (numel(cx)~=numel(cy))
		error("incorrect binary(m) grating definition (1)");
	end

		% check for filling-in and intersections
	if sum(dx.*dy) > 1
		error("incorrect binary grating definition (2)");
	end

	ix = linspace(1,xno-1,xno-1);
	iy = linspace(1,yno-1,yno-1);
	[IX,IY] = meshgrid(ix,iy);
	
	M = zeros(2*yno-1,2*xno-1);
	FE = cellmat(2,1,2*yno-1,2*xno-1);

	for ip = 1:np % loop over nonzero pixels
		fx = (sin(ix*pi*dx(ip))./(pi*ix)).*exp((-2*pi*1i*cx(ip))*ix);
		fy = (sin(iy*pi*dy(ip))./(pi*iy)).*exp((-2*pi*1i*cy(ip))*iy);
		FX = (sin(IX*pi*dx(ip))./(pi*IX)).*exp((-2*pi*1i*cx(ip))*IX);
		FY = (sin(IY*pi*dy(ip))./(pi*IY)).*exp((-2*pi*1i*cy(ip))*IY);

		M = zeros(2*yno-1,2*xno-1);
		
		M(yno,xno) = dx(ip)*dy(ip);

		M(yno+1:2*yno-1,xno) = dx(ip)*fy;
		M(yno-1:-1:1,xno) = conj(M(yno+1:2*yno-1,xno));
		M(yno,xno+1:2*xno-1) = dy(ip)*fx;
		M(yno,xno-1:-1:1) = conj(M(yno,xno+1:2*xno-1));

		M(yno+1:2*yno-1,xno+1:2*xno-1) = FX.*FY;
		M(yno+1:2*yno-1,xno-1:-1:1) = conj(FX).*FY;
		M(yno-1:-1:1,xno+1:2*xno-1) = FX.*conj(FY);
		M(yno-1:-1:1,xno-1:-1:1) = conj(FX.*FY);

		FE{1} = FE{1} + (eps(ip) - eps_m)*M;
		FE{2} = FE{2} + (1/eps(ip) - 1/eps_m)*M;
	end
	FE{1}(yno,xno) = FE{1}(yno,xno) + eps_m;
	FE{2}(yno,xno) = FE{2}(yno,xno) + 1/eps_m;
end
