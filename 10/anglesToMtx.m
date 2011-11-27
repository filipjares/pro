function mtx = anglesToMtx(A, B, C)
%ANGLESTOMTX Converts yaw, pitch and roll to rotation matrix.
%
% mtx = anglesToMtx(A)
% mtx = anglesToMtx(A, B, C)
%
% Arguments:
%		A: if one argumet is used, then it is [yaw; pitch; roll], otherwise yaw
%		B: pitch
%		C: roll
% Return value:
%		mtx: The rotation matrix

if (nargin == 1)
	A = A*pi/180;
	mtx = mtxRotateZ(A(3)) * mtxRotateY(A(2)) * mtxRotateX(A(1));
else
	A = A*pi/180;
	B = B*pi/180;
	C = C*pi/180;
	mtx = mtxRotateZ(C) * mtxRotateY(B) * mtxRotateX(A);
end
	
end

function R = mtxRotateX(alpha)
	R = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
end

function R = mtxRotateY(alpha)
	R = [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
end

function R = mtxRotateZ(alpha)
	R = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
end
