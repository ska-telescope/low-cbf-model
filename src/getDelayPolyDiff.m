function [poly] = getDelayPolyDiff(sky1,sky2,observationTime,polyWindow)
% Gets a 3rd order polynomial approximation to the difference in the functions given by sky1 and sky2 (i.e. sky1 - sky2)
% sky1 and sky2 specify delays as per "getDelayFunctions", as sinusoids with
%  .offset
%  .amplitude
%  .rate
%  .phase
% observationTime is the center time to evaluate the delays over.
% polyWindow is the length of time to evaluate the delays over.
%

t = (observationTime - (polyWindow/2)):(polyWindow/10):(observationTime + (polyWindow/2));

% t_fit is the time interval recentered so observationTime = 0.
t_fit = t + observationTime;
stations = length(sky1.offset);
poly = zeros(stations,4);
%keyboard
for s = 1:stations
    d1 = sky1.offset(s) + sky1.amplitude(s) * sin(sky1.rate(s) * t + sky1.phase(s));
    d2 = sky2.offset(s) + sky2.amplitude(s) * sin(sky2.rate(s) * t + sky2.phase(s));
    p = polyfit(t_fit,d1 - d2,3);
    poly(s,1:4) = p;
end
