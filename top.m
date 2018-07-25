%% Generate station data

% get 
filterTaps = 32;
NFilters = 8;
BWFrac = 1;
tic
[filters] = getInterpFilters(filterTaps,NFilters,BWFrac);
toc

% testSig = sin(0:(pi/4):1000);
% 
% for n = 1:length(testSig)
%     
% end

din = sin(0:pi/3:(1e6 * pi/3 + 500));
din_shift = sin(pi/12:pi/3:(1e6 * pi/3 + 500));
Ts = 1080; % ns
CF = 0;
% delay = [offset, amplitude, frequency, phase] in ns and radians/sec
delay = [Ts*40 + Ts/4 + Ts/2048, 0, 0, 0];
%delay = [Ts*40, 0, 0, 0];
resampledPoints = 1e6;

% unit vector pointing at the sky

tic
[resampled] = resampleNU(din,Ts,CF,delay,resampledPoints);
toc

figure(3);
clf;
hold on;
grid on;
%plot(resampled - din_shift(40:49).','b.-');
plot(resampled,'r.-');
plot(din_shift(40:49),'go-');
title('red resampled');
