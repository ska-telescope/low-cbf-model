function [delayFunctions] = getDelayFunctions(stations,sky,observationTime)
% Generates delay functions for an object in the sky, for all the stations.
% Delay functions are sinusoids, with DC offset, amplitude, phase offset and rate.
%
% Stations should be a struct with fields:
%  .latitude  - latitude of the center point of the telescope, in degrees.
%  .longitude - longitude of the center point of the telescope, in degrees.
%  .altitude  - altitude of the telescope.
%  .offsets   - offset of the stations from the center point, Nx2 array with distance East and North in meters for each station.
%
% Sky should be a struct with the fields:
%  .ascension   - right ascension, in degrees.
%  .declination - declination, in degrees.
%  .rate        - Objects move across the sky at this multiple of the sidereal rate
%  .index       - to index this object in the sky
%  .subIndex    - The center of field for the object has subIndex = 0. This is used for telescope pointing.
%                  subIndex 1 and higher is used for generating the sky.
%
% The output is a structure with fields
%  .offset
%  .amplitude
%  .rate
%  .phase
%  .poly
% Each of these is an array with length equal to the number of stations, except for .poly which is a matrix 
% with dimensions (number of stations)x4, representing a 3rd order polynomial approximation centered on the 
% observation time. The polynomial fit is done to the data in a region of +/- 10 minutes.
%
doplot = 1;
c = 299792458;
siderealSeconds = 23*60*60 + 56*60 + 4.0905; % seconds per sidereal day (23 hours, 56 minutes, 4.0905 seconds)
siderealRate = 2*pi/(siderealSeconds);
% Note : MWA coordinates are
% Latitude: -26.70331940�, Longitude: 116.67081524�
latitude = stations.latitude * pi/180;
% We don't actually use the longitude of the station center, 
% since time is assumed such that a right ascension of -90 degrees is directly overhead at the center of the telescope at t=0.
longitude = stations.longitude * pi/180;  
altitude = stations.altitude; % of the telescope (about 400m for SKA)

% Relative coordinates of the stations
% distance in meters east and north
%stationLoc = [0, 0; 1000, 0; -1000, 0; 1000, 1000];
stationLoc = stations.offsets;
stations = size(stationLoc,1);

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
delayFunctions.offset = zeros(stations,1);
delayFunctions.amplitude = zeros(stations,1);
delayFunctions.phase = zeros(stations,1);
delayFunctions.rate = zeros(stations,1);
delayFunctions.poly = zeros(stations,4);

for s = 1:stations
    % Convert east,north locations to actual latitude, longitude
    long = stationLoc(s,1)/rAxis;  % Only use the offset for the longitude; This sets the zero for time such that a right ascension of 0 is on the horizon at t = 0
    lat = latitude - stationLoc(s,2)/rCenter;
    
    % Distance from the earths axis for this station
    axialRadius = rCenter * cos(lat);
    
    ra = sky.ascension * pi/180;  % right ascension in radians
    dec = sky.declination * pi/180; % declination in radians
    % Unit vector pointing at the object
    v1 = [sin(dec), cos(ra)*cos(dec), sin(ra) * cos(dec)];
    % Three components of the dot product with a vector pointing at the station
    % the vertical component of the dot product is constant
    delayFunctions.offset(s) = (1/c) * rCenter * sin(lat) * v1(1);
    % The other two components vary sinusoidally with time.
    % modelled as A*sin(time + P)
    delayFunctions.amplitude(s) = (1/c) * axialRadius * cos(dec);   % component in the horizontal plane [which = sqrt((axialRadius * v1(2))^2 + (axialRadius * v1(3))^2)]
    % Phase offset P.
    % To place an object with right ascension(ra)=0 on the horizon at t=0,
    % the angle between the vector pointing at the telescope and the object is (pi/2 + ra - long - time)
    % so the dot product is cos(pi/2 + ra - long - time),
    % We use sin(time + offset), so using cos(x) = sin(pi/2 - x) gives
    % cos(pi/2 + ra - long - time) = sin(long - ra + time)
    delayFunctions.phase(s) = long - ra;
    delayFunctions.rate(s) = sky.rate * siderealRate;
end

% Generate 3rd order polynomial approximations.
% Note : CSP interface spec states "each polynomial is provided as an array of at most 4 coefficients"
polyWindow = 10; % seconds
t = (observationTime - (polyWindow/2)):1:(observationTime + (polyWindow/2));

% t_fit is the time interval recentered so observationTime = 0.
t_fit = t + observationTime;

for s = 1:stations
    d = delayFunctions.offset(s) + delayFunctions.amplitude(s) * sin(delayFunctions.rate(s) * t + delayFunctions.phase(s));
    p = polyfit(t_fit,d,3);
    delayFunctions.poly(s,1:4) = p;
end

if (doplot)
    clist = 'rgbcmyk';
    % Plot the full delay functions
    figure(1);
    clf;
    hold on;
    grid on;
    for s = 1:stations
        plot(delayFunctions.offset(s) + delayFunctions.amplitude(s) * sin(delayFunctions.rate(s) * (0:100:siderealSeconds) + delayFunctions.phase(s)),[clist(mod(s,7) + 1) '.-']);
    end
    title('Delay functions (full 24 hours)');
    
    % Plot the portion of the delay functions estimated by the polynomial approximation
    figure(2);
    clf;
    hold on;
    grid on;
    tplot = 10 * (((-polyWindow/2):0.1:(polyWindow/2)));
    for s = 1:stations
        plot(tplot,delayFunctions.offset(s) + delayFunctions.amplitude(s) * sin(delayFunctions.rate(s) * tplot + delayFunctions.phase(s)),[clist(mod(s,7) + 1) '.-']);
        plot(tplot,delayFunctions.poly(s,1)*(tplot.^3) + delayFunctions.poly(s,2)*(tplot.^2) + delayFunctions.poly(s,3)*(tplot) + delayFunctions.poly(s,4),[clist(mod(s,7) + 1) 'o-']);
    end
    title(['Delay Functions (centered on observation time = ' num2str(observationTime) '), . true, o fit'])
    
    % Plot the difference between the polynomial estimate and the true delay 
    figure(3);
    clf;
    hold on;
    grid on;
    for s = 1:stations
        % True Delay
        d1 = delayFunctions.offset(s) + delayFunctions.amplitude(s) * sin(delayFunctions.rate(s) * tplot + delayFunctions.phase(s));
        % Polynomial Approximation
        d2_4 = delayFunctions.poly(s,1)*(tplot.^3) + delayFunctions.poly(s,2)*(tplot.^2) + delayFunctions.poly(s,3)*(tplot) + delayFunctions.poly(s,4);
        %d2_3 = delayFunctions.poly(s,2)*(tplot.^2) + delayFunctions.poly(s,3)*(tplot) + delayFunctions.poly(s,4);
        %d2_2 = delayFunctions.poly(s,3)*(tplot) + delayFunctions.poly(s,4);
        plot(tplot,d1 - d2_4,[clist(mod(s,7) + 1) '.-']);
        %plot(tplot,d1-d2_3,[clist(mod(s,7) + 1) '*-']);
        %plot(tplot,d1-d2_2,[clist(mod(s,7) + 1) 's-']);
    end
    title('True delay - polynomial approximation');
    
end

