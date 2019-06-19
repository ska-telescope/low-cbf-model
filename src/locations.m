function [P] = getDelayFunctions(stations,sky)
% Stations and sky are read from configuration files.
% stations should have fields:
%  -latitude
%  -longitude
%  -altitude
%  -offsets
% sky should have fields:
%  - 
%
doplot = 0;
% MWA coordinates are
% Latitude: -26.70331940°, Longitude: 116.67081524°
latitude = stations.latitude * pi/180;
longitude = stations.longitude * pi/180;
altitude = stations.altitude; % of the telescope (about 400m for SKA)

% Relative coordinates of the stations
% distance in meters east and north
%stationLoc = [0, 0; 1000, 0; -1000, 0; 1000, 1000];
stationLoc = stations.offsets;
stations = size(stationLoc,1);

% Locations of objects in the sky (right ascension and declination, degrees)
% Note that a right ascension of 0 places the object on the horizon at t=0 for the center of the telescope.
% A right ascension of -90 places the object directly overhead at t = 0.
%objectLocations = [0, -20; 90, -30; 180, 0];
objectLocations = [-45, -63.2967];
objects = size(objectLocations,1);

% Distance of the telescope center from the center of the earth.
a = 6378.1370e3;
b = 6356.7523e3;
rCenter = altitude + sqrt( ((a^2 * cos(latitude))^2 + (b^2 * sin(latitude))^2) / ((a * cos(latitude))^2 + (b * sin(latitude))^2));
% Distance of the telescope center from the Earths axis
rAxis = rCenter * cos(latitude);

% For each station and each object, generate offset d0, amplitude A and phase P such that the
% delay of the object at the station is given by 
%  D = d0 + A * sin(R*w*t + P)
% where w is sidereal rate and t is time.
d0 = zeros(stations,objects);
A = zeros(stations,objects);
P = zeros(stations,objects);
R = zeros(stations,objects);

for  s = 1:stations
    % Convert east,north locations to actual latitude, longitude
    long = stationLoc(s,1)/rAxis;  % Only use the offset for the longitude; This sets the zero for time such that a right ascension of 0 is on the horizon at t = 0
    lat = latitude - stationLoc(s,2)/rCenter;
    
    % Distance from the earths axis for this station
    axialRadius = rCenter * cos(lat);
    
    for object = 1:objects
        ra = objectLocations(object,1) * pi/180;  % right ascension in radians
        dec = objectLocations(object,2) * pi/180; % declination in radians
        % Unit vector pointing at the object
        v1 = [sin(dec), cos(ra)*cos(dec), sin(ra) * cos(dec)];
        % Three components of the dot product with a vector pointing at the station
        % the vertical component of the dot product is constant
        d0(s,object) = rCenter * sin(lat) * v1(1);
        % The other two components vary sinusoidally with time.
        % modelled as A*sin(time + P)
        A(s,object) = axialRadius * cos(dec);   % component in the horizontal plane [which = sqrt((axialRadius * v1(2))^2 + (axialRadius * v1(3))^2)]
        % Phase offset P. 
        % To place an object with right ascension(ra)=0 on the horizon at t=0,
        % the angle between the vector pointing at the telescope and the object is (pi/2 + ra - long - time)
        % so the dot product is cos(pi/2 + ra - long - time), 
        % We use sin(time + offset), so using cos(x) = sin(pi/2 - x) gives
        % cos(pi/2 + ra - long - time) = sin(long - ra + time) 
        P(s,object) = long - ra;
        R(s,object) = objectRate;
    end
    
end


if (doplot)
    figure(1);
    clf;
    hold on;
    grid on;
    for s = 1:stations
        for object = 1:objects
            plot(d0(s,object) + A(s,object) * sin((0:(2*pi/100):2*pi) + P(s,object)),'r.-');
        end
    end
end