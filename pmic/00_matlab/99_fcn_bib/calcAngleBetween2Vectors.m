function [ang, vn] = calcAngleBetween2Vectors(v1, v2)
% positive angle around vn rotates v1 into v2 whereas |vn| = 1

% https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
vn = cross(v1, v2);
vn_norm = norm(vn);
if vn_norm > 1e-6
    vn = vn / vn_norm;
    ang = atan2( cross(v1, v2).' * vn, v1.' * v2 );
else
    ang = 0.0;
end

end

