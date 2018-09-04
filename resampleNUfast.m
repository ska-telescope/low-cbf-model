function [resampled] = resampleNUfast(din,Ts,CF,delay,resampledPoints)
% Non-Uniform Resampling.
% Resamples the input signal using a non-uniform rate specified by delay.
% run time optimised version.
% Inputs:
%  din : complex baseband input data.
%  Ts : Sampling period in nanoseconds. The same value is used for the input data (din) and the output (resampled), except for the change created by the resampling.
%  CF : center frequency in Hz for the input data (impacts the doppler shift)
%  delay : Time offset for sampling the signal, specified as a sinusoid with an array of 4 numbers
%          [offset, amplitude, frequency, phase, slope]
%          delay = offset + amplitude * sin(frequency * t + phase) + slope*t
%           - offset and amplitude are in ns.
%           - slope is in ns/s
%           - phase is the phase at t=0, in radians.
%           - frequency is in radians per second. e.g. for sidereal rate frequency = 2*pi/(24*60*60 - 235.9) 
%             (note a sidereal day is 3 minutes 55.9 seconds = 235.9 seconds short)
%          Note delay should always be greater than 0, as the first sample taken from din is at delay.
%  resampledPoints : Number of samples to calculate. Error will occur if there are not enough input samples.
%  
% Outputs:
%  resampled : resampled version of the input data.

% Get interpolation filters
filterTaps = 32; % Pick an even number for the code below to work. For filterTaps = 32, interpolation occurs between the samples aligned to filtertap 16 and 17 (using matlab indexing from 1)
NFilters = 512;  % 32 taps x 512 filters will fit in core i7 L2 cache - should be faster than if it doesn't... actually doesn't make much difference because the same filter gets used for at least 1000 successive points.
BWFrac = 1;
persistent filters;
persistent filtersR; % time reversed version of filters, used for convolution.
if isempty(filters)
    [filters] = getInterpFilters(filterTaps,NFilters,BWFrac);
    filtersR = filters(:,filterTaps:-1:1);
end
ds = size(din);
if (ds(2) ~= 1)
    warning('data input must be a row vector');
    din = din.';
end

resampled = zeros(resampledPoints,1);

p_all = (0:(resampledPoints-1)).';
DelayOffset_all = (delay(1) + delay(2) * sin(delay(3) * (p_all)*Ts*1e-9 + delay(4)) + (delay(5) * p_all * Ts * 1e-9))/Ts;
DelayOffsetInt_all = floor(DelayOffset_all);
DelayOffsetFrac_all = DelayOffset_all - DelayOffsetInt_all;
filterSelect_all = round(DelayOffsetFrac_all * NFilters) + 1;
changePoints = find(diff(filterSelect_all)) + 1;  % Indexes at which the interpolation filter changes.

%keyboard
if isempty(changePoints)
    % Same filter throughout
    DelayOffsetInt_part = DelayOffsetInt_all(1);
    filterSelect_part = filterSelect_all(1);
    
    SampleOffsetStart = DelayOffsetInt_part(1) - filterTaps/2 + 1;
    SampleOffsetEnd = SampleOffsetStart + resampledPoints + filterTaps - 2;
    resampled = conv(din(SampleOffsetStart:SampleOffsetEnd),filtersR(filterSelect_part(1),:),'valid');
    
else
    % Multiple different filter are required.
    DelayOffsetInt_part = DelayOffsetInt_all(changePoints - 1);
    filterSelect_part = filterSelect_all(changePoints - 1);
    
    % First block is point 1 to point changePoints(1) - 1, using filterSelect_all(changePoints(1) - 1)
    SampleOffsetStart = DelayOffsetInt_part(1) - filterTaps/2 + 1;
    SampleOffsetEnd = SampleOffsetStart + changePoints(1) - 1 + filterTaps - 2;
    resampled(1:(changePoints(1) - 1)) = conv(din(SampleOffsetStart:SampleOffsetEnd),filtersR(filterSelect_part(1),:),'valid');
    
    for p = 2:length(changePoints)
        % Everything from changePoints(p-1) to (changePoints(p) - 1) uses the filter filterSelect_all(changePoints(p) - 1)
        SampleOffsetStart = DelayOffsetInt_part(p) + changePoints(p-1) - filterTaps/2 + 1 - 1;
        SampleOffsetEnd = DelayOffsetInt_part(p) + changePoints(p) - 1 + filterTaps/2 - 1 - 1;
        totalPoints = (changePoints(p) - 1) - changePoints(p-1);
        %keyboard
        resampled(changePoints(p-1):(changePoints(p-1) + totalPoints - 1)) = conv(din(SampleOffsetStart:SampleOffsetEnd),filtersR(filterSelect_part(p),:),'valid');
    end
end
resampled = resampled .* exp(1i * 2*pi*DelayOffsetFrac_all * Ts * 1e-9 * CF);

%toc
%tic

% resampled2 = zeros(resampledPoints,1);
% fsel = zeros(resampledPoints,1);
% 
% for p = 1:resampledPoints
%     % Delay in fractions of a sample 
%     DelayOffset = (delay(1) + delay(2) * sin(delay(3) * (p-1)*Ts*1e-9 + delay(4)))/Ts;
%     DelayOffsetInt = floor(DelayOffset);
%     DelayOffsetFrac = DelayOffset - DelayOffsetInt;
%     % Base delay offset in samples
%     SampleOffset = (p-1) + DelayOffsetInt;
%     
%     % Interpolate 
%     filterSelect = round(DelayOffsetFrac * NFilters);
%     fsel(p) = filterSelect;
%     resampled2(p) = filters(filterSelect+1,:) * din((SampleOffset - filterTaps/2 + 1):(SampleOffset + filterTaps/2));
%     %resampled(p) = sum(filters(filterSelect+1,:) .* din((SampleOffset - filterTaps/2 + 1):(SampleOffset + filterTaps/2)));
%     
%     % Phase Rotation
%     % amount of time past the nominal sampling point is (DelayOffsetFrac * Ts)*1e-9 ns
%     % then multiply by the center frequency (CF) to get the number of rotations.
%     Phase = 2*pi*DelayOffsetFrac * Ts * 1e-9 * CF;
% %    Phase = 0;
%     resampled2(p) = resampled2(p) * exp(1i * Phase);    
% end

%toc

%keyboard

