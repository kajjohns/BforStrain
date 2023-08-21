function xy = llh2local(llh, origin)
% LLH2LOCAL Converts from longitude and latitude to local coordinates.
%   xy = llh2local(llh, origin)
%
%   Converts from longitude and latitude to local coordinates given an origin. 
%   'llh' (lon; lat; height) and 'origin' should be in decimal degrees. Note that 
%   heights are ignored, and xy is in km.
%
%   Inputs:
%   - llh: A 3xN matrix of [lon; lat; height] in decimal degrees.
%   - origin: A 1x2 vector [lon, lat] of the origin in decimal degrees.
%
%   Output:
%   - xy: A 2xN matrix of local coordinates [x; y] in kilometers.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
% 
%   Sept 7, 2000  Peter Cervelli		Original Code
%   Oct 20, 2000  Jessica Murray        Changed name from DM_llh2local to 
%                                       llh2local for use with non-DM functions;
%                                       Added to help message to clarify order
%                                       of 'llh' (i.e., lon, lat, height).
%   Dec. 6, 2000  Jessica Murray        Clarified help to show that llh 
%                                       is a column vector
%   Aug 12, 2023  Nicolas Castro        Updated WGS84 excentricity value,
%                                       improved variable naming, code
%                                       readability and handling of special
%                                       cases to avoid numerical issues.
%                                       Improved computation of the reduced
%                                       and prime vertical radii of curvature 
%                                       by using the average latitude and its 
%                                       sine/cosine values.
%
%-------------------------------------------------------------

% Set ellipsoid constants (WGS84)
a = 6378137.0; % semi-major axis
e = 0.0818191908426; % eccentricity

% Convert to radians
llh = llh * pi / 180;
origin = origin * pi / 180;

% Calculate differences in longitudes and latitudes
dlat = llh(2, :) - origin(2);
dlon = llh(1, :) - origin(1);

% Ensure longitude differences are within [-180, 180]
dlon = mod(dlon + pi, 2*pi) - pi;

% Handle special case of latitude = 0
z = abs(llh(2, :)) > 1e-10; % Use a small tolerance for latitude comparison
dlat(~z) = 0; % Set latitude difference to zero for points at the equator

% Compute reduced latitudes and their sines and cosines
lat_avg = (llh(2, :) + origin(2)) / 2;
lat_avg_sin = sin(lat_avg);
lat_avg_cos = cos(lat_avg);

% Compute the meridian radius of curvature
M = a * (1 - e^2) ./ (sqrt(1 - e^2 * lat_avg_sin.^2)).^3;

% Compute the prime vertical radius of curvature
N = a ./ sqrt(1 - e^2 * lat_avg_sin.^2);

% Compute local coordinates
xy = zeros(2, size(llh, 2));
xy(1, :) = dlon .* N .* lat_avg_cos;
xy(2, :) = dlat .* M;

% Convert to km
xy = xy / 1000;
end
