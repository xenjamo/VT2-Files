function x_ecefr = transformWGS84ToECEF_R(lat, lon, h)

    a = 6378137.0;         % WGS | 84 Earth semimajor axis
    e = 8.1819191 * 1e-2;  % eccentricity

    N = a ./ sqrt(1 - e^2*sin(lat).^2);
    x_ecefr = [(h + N) .* cos(lat) .* cos(lon), ...
               (h + N) .* cos(lat) .* sin(lon), ...
               (h + (1 - e^2)*N) .* sin(lat)];
end