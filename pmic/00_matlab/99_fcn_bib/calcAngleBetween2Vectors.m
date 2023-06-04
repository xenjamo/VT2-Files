function [ang, vn] = calcAngleBetween2Vectors(v1, v2)
% positive angle around vn rotates v1 into v2 whereas |vn| = 1

% https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
vn = cross(v1, v2);
vn = vn / norm(vn);
ang = atan2( cross(v1, v2).' * vn, v1.' * v2 );

end

