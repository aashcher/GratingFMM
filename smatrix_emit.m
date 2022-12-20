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
	% (+,-) -> (-below,+above)
%%
function Vemi = smatrix_emit(S1, S2, Vloc)
	no = size(S1,1);
	Vemi = Vloc;
	if (numel(size(S1)) == 4) && (numel(size(S2)) == 4)
		Vemi(:,1) = ( S1(:,:,1,2) / (eye(no) - S2(:,:,1,1)*S1(:,:,2,2)) ) * ...
								( S2(:,:,1,1) * Vloc(:,1) + Vloc(:,2) );
		Vemi(:,2) = ( S2(:,:,2,1) / (eye(no) - S1(:,:,2,2)*S2(:,:,1,1)) ) * ...
								( Vloc(:,1) + S1(:,:,2,2) * Vloc(:,2) );
	elseif (numel(size(S1)) == 3) && (numel(size(S2)) == 4)
		Vemi(:,1) = transpose( transpose(S1(:,1,2)) / (eye(no) - S2(:,:,1,1).*transpose(S1(:,2,2))) ) .* ...
								( S2(:,:,1,1) * Vloc(:,1) + Vloc(:,2) );
		Vemi(:,2) = ( S2(:,:,2,1) / (eye(no) - S1(:,2,2).*S2(:,:,1,1)) ) * ...
								( Vloc(:,1) + S1(:,2,2) .* Vloc(:,2) );
	elseif (numel(size(S1)) == 4) && (numel(size(S2)) == 3)
		Vemi(:,1) = ( S1(:,:,1,2) / (eye(no) - S2(:,1,1).*S1(:,:,2,2)) ) * ...
								( S2(:,1,1) .* Vloc(:,1) + Vloc(:,2) );
		Vemi(:,2) = transpose( transpose(S2(:,2,1)) / (eye(no) - S1(:,:,2,2).*transpose(S2(:,1,1))) ) .* ...
								( Vloc(:,1) + S1(:,:,2,2) * Vloc(:,2) );
	elseif (numel(size(S1)) == 3) && (numel(size(S2)) == 3)
		Vemi(:,1) = ( S1(:,1,2) ./ (1 - S2(:,1,1).*S1(:,2,2)) ) .* ...
								( S2(:,1,1) .* Vloc(:,1) + Vloc(:,2) );
		Vemi(:,2) = ( S2(:,2,1) ./ (1 - S1(:,2,2).*S2(:,1,1)) ) .* ...
								( Vloc(:,1) + S1(:,2,2) .* Vloc(:,2) );
	else
		error('smatrix_emit: unknown matrix dimensions');
	end
end