function A = removeBorders3(A, sizebord, nrows, ncols, v)

% puts to value v pixels on the borders (of size sizebord) of the 3D image A 

nim = length(A)/(ncols*nrows); 

A(1 : sizebord*nrows*ncols) = v; % put to v first "sizebord" slices
A(end-(sizebord*nrows*ncols)+1 : end) = v; % put to v last "sizebord" slices


npx2d = nrows*ncols; % number of pixels in every slice
for k = 2 : nim - sizebord
	
	idx_begIm = ((k-1)*npx2d) + 1;
	idx_endIm = k*npx2d;
	
% 	[idx_begIm , idx_begIm + (sizebord*nrows) - 1]
	A(idx_begIm : idx_begIm + (sizebord*nrows) - 1) = v; % put to v first "sizebord" columns of image k
	
% 	[idx_endIm - (sizebord*nrows) + 1, idx_endIm]
	A(idx_endIm - (sizebord*nrows) + 1 : idx_endIm) = v; % put to v last "sizebord" columns of image k


	for i = 2 : ncols - sizebord
        idx_begCol = idx_begIm + ((i-1)*nrows);
		idx_endCol = idx_begIm + i*nrows - 1;
		A(idx_begCol : idx_begCol + sizebord - 1) = v;
		A(idx_endCol - sizebord + 1  : idx_endCol) = v; 
	end
end

end 