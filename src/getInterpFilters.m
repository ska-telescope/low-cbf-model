function [filters] = getInterpFilters(taps,N,BWFrac)
% Get a set of bandlimited filters for interpolating.
% Generates a set of filters, where each filter interpolates to a
% particular offset from 0 to 1.
% Inputs :
%  taps : Number of filter taps to use
%  N : Number of filters to generate - 1. Interpolation offsets for the filters are
%      0:(1/N):1, so one extra filter is generated for the end point of the interval
%  BWFrac : fraction of the full bandwidth used for interpolation.
%           set to 1 for full bandwidth (i.e. no attenuation of the signal)
%
doplot = 0;
sincPoints = (-1 * floor((taps-1)/2)):((-1 * floor((taps-1)/2)) + taps - 1);
windowFull = (1/2) * (1 + cos(-pi :(2*pi/(taps*N)):pi));
filters = zeros(N+1,taps);
for f = 1:(N+1)
    %disp('---')
    %disp('f = ');
    %disp(sinc((sincPoints  - (f-1)/N) * BWFrac));
    %disp('window = ');
    %disp(windowFull((N+2-f):N:(N+2-f + N*taps - 1)));
    filters(f,:) = sinc((sincPoints  - (f-1)/N) * BWFrac) .* windowFull((N+2-f):N:(N+2-f + N*taps - 1));
    filters(f,:) = filters(f,:) / sum(filters(f,:));
end

if (doplot)
    figure(2);
    clf;
    hold on;
    grid on;
    for f = 1:(N+1)
       plot(abs(fft([filters(f,:) zeros(1,1000)])),'r.-')
    end
end


