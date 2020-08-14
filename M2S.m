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
%% convert a matrix to S-matrix format
function S = M2S(M)
	no = size(M,1)/2;
	S = zeros(no,no,2,2);
	S(:,:,1,1) = M(1:no,1:no);
	S(:,:,1,2) = M(1:no,(no+1):(2*no));
	S(:,:,2,1) = M((no+1):(2*no),1:no);
	S(:,:,2,2) = M((no+1):(2*no),(no+1):(2*no));
end