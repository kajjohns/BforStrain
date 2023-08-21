function llh = local2llh(xy, origin)
% LOCAL2LLH Converts from local coordinates to longitude and latitude.
%   llh = local2llh(xy, origin)
%
%   Converts from local coordinates to longitude and latitude given the [lon, lat] 
%   of an origin. 'origin' should be in decimal degrees. Note that heights are ignored, 
%   and xy is in km. Output is [lon, lat, height] in decimal degrees. This is an 
%   iterative solution for the inverse of a polyconic projection.
%
%   Inputs:
%   - xy: A 2xN matrix of local coordinates [x; y] in kilometers.
%   - origin: A 1x2 vector [lon, lat] of the origin in decimal degrees.
%
%   Output:
%   - llh: A 3xN matrix of [lon, lat, height] in decimal degrees.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 12, 2023  Nicolas Castro        Improved convergence check, array
%                                       prealocation and precomputation of
%                                       common terms
%   Aug 23, 2001  Jessica Murray        Clarification to help.
%
%   Apr 4, 2001   Peter Cervelli        Added failsafe to avoid
%                                       infinite loop because of
%                                       covergence failure.
%   Sep 7, 2000   Peter Cervelli		Original Code
%
%-------------------------------------------------------------

% Set ellipsoid constants (WGS84)
a = 6378137.0; % semi-major axis
e = 0.0818191908426; % eccentricity

% Convert to radians / meters
xy = xy * 1000;
origin = origin * pi / 180;

% Precompute some common terms
e2 = e^2;
e4 = e^4;
e6 = e^6;

% Iterate to perform inverse projection
M0 = a * ((1 - e2/4 - 3*e4/64 - 5*e6/256) * origin(2) - ...
          (3*e2/8 + 3*e4/32 + 45*e6/1024) * sin(2 * origin(2)) + ...
          (15*e4/256 + 45*e6/1024) * sin(4 * origin(2)) - ...
          (35*e6/3072) * sin(6 * origin(2)));

z = xy(2, :) ~= -M0;

A = (M0 + xy(2, z)) / a;
B = xy(1, z).^2 / a^2 + A.^2;

llh = zeros(3, size(xy, 2));
llh(2, z) = A;

% Iterate for inverse projection
for iter = 1:100
    C = sqrt(1 - e2 * sin(llh(2, z)).^2) .* tan(llh(2, z));

    M = a * ((1 - e2/4 - 3*e4/64 - 5*e6/256) * llh(2, z) - ...
             (3*e2/8 + 3*e4/32 + 45*e6/1024) * sin(2 * llh(2, z)) + ...
             (15*e4/256 + 45*e6/1024) * sin(4 * llh(2, z)) - ...
             (35*e6/3072) * sin(6 * llh(2, z)));

    Mn = 1 - e2/4 - 3*e4/64 - 5*e6/256 - ...
         2 * (3*e2/8 + 3*e4/32 + 45*e6/1024) * cos(2 * llh(2, z)) + ...
         4 * (15*e4/256 + 45*e6/1024) * cos(4 * llh(2, z)) + ...
         -6 * (35*e6/3072) * cos(6 * llh(2, z));

    Ma = M / a;
   
    delta = -(A .* (C .* Ma + 1) - Ma - 0.5 * (Ma.^2 + B) .* C) ./ ...
            (e2 * sin(2 * llh(2, z)) .* (Ma.^2 + B - 2 * A .* Ma) ./ (4 * C) + (A - Ma) .* (C .* Mn - 2./sin(2 * llh(2, z))) - Mn);

    llh(2, z) = llh(2, z) + delta;

    if all(abs(delta) < 1e-8)
        break; % Convergence reached
    end
end

% Calculate longitude and latitude for the entire array
llh(1, :) = asin(xy(1, :) .* C ./ a) ./ sin(llh(2, :)) + origin(1);

% Handle special case of latitude = 0
llh(1, ~z) = xy(1, ~z) / a + origin(1);
llh(2, ~z) = 0;

% Convert back to decimal degrees
llh = llh * 180 / pi;
end
