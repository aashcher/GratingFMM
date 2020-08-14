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
%% description:
% multiplication of two S-matrices
%% input:
% SM1, SM2: S-matrices of size (n,n,2,2) or (n,2,2)
% the multiplication order is meaningful!
%% output:
% SM: S-matrix of size (n,n,2,2) or (n,2,2)
%% implementation:
function [SM] = mul_SM(SM1, SM2)
	if (numel(size(SM1)) == 4) && (numel(size(SM2)) == 4) % both S-matrices are full
		if ~isequal(size(SM1),size(SM2))
			error("mul_SM: different size input");
		end
		n = size(SM1,1);
		SM = zeros(n,n,2,2);

		TM = -SM2(:,:,1,1)*SM1(:,:,2,2);
		TM(1:n+1:end) = TM(1:n+1:end) + 1;
		TM = SM1(:,:,1,2)/TM;
		SM(:,:,1,2) = TM*SM2(:,:,1,2);
		SM(:,:,1,1) = SM1(:,:,1,1) + TM*SM2(:,:,1,1)*SM1(:,:,2,1);

		TM = -SM1(:,:,2,2)*SM2(:,:,1,1);
		TM(1:n+1:end) = TM(1:n+1:end) + 1;
		TM = SM2(:,:,2,1)/TM;
		SM(:,:,2,1) = TM*SM1(:,:,2,1);
		SM(:,:,2,2) = SM2(:,:,2,2) + TM*SM1(:,:,2,2)*SM2(:,:,1,2);

	elseif (numel(size(SM1)) == 3) && (numel(size(SM2)) == 4) % first S-matrix is diagonal
		if size(SM1,1) ~= size(SM2,1)
			error("mul_SM: different size input");
		end
		n = size(SM1,1);
		SM = zeros(n,n,2,2);

		TM = -SM2(:,:,1,1).*transpose(SM1(:,2,2));
		TM(1:n+1:end) = TM(1:n+1:end) + 1;
		TM = diag(SM1(:,1,2))/TM;
		SM(:,:,1,2) = TM*SM2(:,:,1,2);
		SM(:,:,1,1) = diag(SM1(:,1,1)) + TM*(SM2(:,:,1,1).*transpose(SM1(:,2,1)));

		TM = -SM1(:,2,2).*SM2(:,:,1,1);
		TM(1:n+1:end) = TM(1:n+1:end) + 1;
		TM = SM2(:,:,2,1)/TM;
		SM(:,:,2,1) = TM.*transpose(SM1(:,2,1));
		SM(:,:,2,2) = SM2(:,:,2,2) + TM*(SM1(:,2,2).*SM2(:,:,1,2));

	elseif (numel(size(SM1)) == 4) && (numel(size(SM2)) == 3) % second S-matrix is diagonal
		if size(SM1,2) ~= size(SM2,1)
			error("mul_SM_SMD: different size input");
		end
		n = size(SM1,2);
		SM = zeros(n,n,2,2);

		TM = -SM2(:,1,1).*SM1(:,:,2,2);
		TM(1:n+1:end) = TM(1:n+1:end) + 1;
		TM = SM1(:,:,1,2)/TM;
		SM(:,:,1,2) = TM.*transpose(SM2(:,1,2));
		SM(:,:,1,1) = SM1(:,:,1,1) + TM*(SM2(:,1,1).*SM1(:,:,2,1));

		TM = -SM1(:,:,2,2).*transpose(SM2(:,1,1));
		TM(1:n+1:end) = TM(1:n+1:end) + 1;
		TM = diag(SM2(:,2,1))/TM;
		SM(:,:,2,1) = TM*SM1(:,:,2,1);
		SM(:,:,2,2) = diag(SM2(:,2,2)) + TM*(SM1(:,:,2,2).*transpose(SM2(:,1,2)));

	elseif (numel(size(SM1)) == 3) && (numel(size(SM2)) == 3) % both S-matrices are diagonal
		if ~isequal(size(SM1),size(SM2))
			error("mul_SM: different size input");
		end
		n = size(SM1,1);
		SM = zeros(n,2,2);

		TM = SM1(:,1,2)./(1 - SM1(:,2,2).*SM2(:,1,1));
		SM(:,1,2) = SM2(:,1,2).*TM;
		SM(:,1,1) = SM1(:,1,1) + TM.*SM2(:,1,1).*SM1(:,2,1);

		TM = SM2(:,2,1)./(1 - SM2(:,1,1).*SM1(:,2,2));
		SM(:,2,1) = SM1(:,2,1).*TM;
		SM(:,2,2) = SM2(:,2,2) + SM2(:,1,2).*SM1(:,2,2).*TM;

	else
		error("mul_SM: incorrect input size");

	end
end
%
% end of mul_SM
%