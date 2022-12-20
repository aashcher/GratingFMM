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
function Vdif = smatrix_diffract(S, Vinc)
	Vdif = Vinc;
	Vdif(:,1) = S(:,:,1,1)*Vinc(:,1) + S(:,:,1,2)*Vinc(:,2);
	Vdif(:,2) = S(:,:,2,1)*Vinc(:,1) + S(:,:,2,2)*Vinc(:,2);
end