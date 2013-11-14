function d = tangent_d(A, B, ~)

% Modified to return Tangent distance (in place of L2)
%-----------------------------------------------------
% tangent_d - computes tangent distance matrix
%
% E = L2_distance(A,B)
%
%    A - (DxM) matrix 
%    B - (DxN) matrix
% 
% Returns:
%    E - (MxN) tangent distances between vectors in A and B
%
% 
%
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    d = tangent_d(A,B);
%
%
  

if (nargin < 2)
   error('Not enough input arguments');
end

if (size(A,1) ~= size(B,1))
   error('A and B should be of same dimensionality');
end
if ~(isreal(A) && isreal(B))
   disp('Warning: running distance.m with imaginary numbers.  Results may be off.'); 
end

disp('Computing the tangent distances...');
d = zeros(size(A,2),size(B,2));
for i =1:size(A,2)
	if(mod(i,10) == 0)
		fprintf(1,'.');
	end
	for j = 1: size(B,2)
		d(i,j) = td(A(:,i),B(:,j));
	end
end