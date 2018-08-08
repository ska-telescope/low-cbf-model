function [resampled] = resampleNU(din,Ts,CF,delay,resampledPoints)
% Non-Uniform Resampling.
% Resamples the input signal using a non-uniform rate specified by delay.
% Inputs:
%  din : complex baseband input data.
%  Ts : Sampling period in nanoseconds. The same value is used for the input data (din) and the output (resampled), except for the change created by the resampling.
%  CF : center frequency in Hz for the input data (impacts the doppler shift)
%  delay : Time offset for sampling the signal, specified as a sinusoid with an array of 4 numbers
%          [offset, amplitude, frequency, phase]
%          delay = offset + amplitude * sin(frequency * t + phase)
%           - offset and amplitude are in ns.
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
NFilters = 512;  % 32 taps x 512 filters will fit in core i7 L2 cache - should be faster than if it doesn't...
BWFrac = 1;
[filters] = getInterpFilters(filterTaps,NFilters,BWFrac);
ds = size(din);
if (ds(2) ~= 1)
    warning('data input must be a row vector');
    din = din.';
end

resampled = zeros(resampledPoints,1);
%Phase = 1.21341;

for p = 1:resampledPoints
    % Delay in fractions of a sample 
    DelayOffset = (delay(1) + delay(2) * sin(delay(3) * (p-1)*Ts*1e-9 + delay(4)))/Ts;
    DelayOffsetInt = floor(DelayOffset);
    DelayOffsetFrac = DelayOffset - DelayOffsetInt;
    % Base delay offset in samples
    SampleOffset = (p-1) + DelayOffsetInt;
    
    % Interpolate 
    filterSelect = round(DelayOffsetFrac * NFilters);
    resampled(p) = filters(filterSelect+1,:) * din((SampleOffset - filterTaps/2 + 1):(SampleOffset + filterTaps/2));
    %resampled(p) = sum(filters(filterSelect+1,:) .* din((SampleOffset - filterTaps/2 + 1):(SampleOffset + filterTaps/2)));
    
    % Phase Rotation
    % amount of time past the nominal sampling point is (DelayOffsetFrac * Ts)*1e-9 ns
    % then multiply by the center frequency (CF) to get the number of rotations.
    Phase = 2*pi*DelayOffsetFrac * Ts * 1e-9 * CF;
    resampled(p) = resampled(p) * exp(1i * Phase);
    
end

