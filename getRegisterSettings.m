function [LRUreg] = getRegisterSettings(arrayConfigFull,modelConfig)
%
% Generate the register settings in the firmware from the LMC configuration
% Returns LRUreg, an array of structs, where each struct has the registers for a single LRU ( = 1 FPGA).
% There is a field for each register setting. These are named according the module name in the firmware,
% Also see the module level documentation on confluence.
%
%  (1) LFAASPEADDecoder
%      - .virtualChannel 
%          * This is a lookup table with 384 entries used to assign the virtual channel number.
%  (2) LocalDoppler
%      - .phaseOffset
%      - .phaseIncrement
%          * Doppler corrections are derived from the delay polynomials
%  (3) CoarseCornerTurner
%      - .
%          * Sample delays for each station & channel
%  (4) CorrelatorFilterbank
%      - FIR Filter Coefficients
%  (5) PSSFilterbank
%      - FIR Filter Coefficients
%  (6) PSTFilterbank
%      - FIR Filter Coefficients
%  (4) Fine Delay
%      - phase slopes 
% 
%  (5) PSS beamformer
%      - Weights
%
%  (6) PST beamformer
%      - Weights
%
%  (7) Packetiser config (SPEAD, PSS, VDIF)
%      - MAC/IP/UDP addresses etc.
%

for LRU = 1:modelConfig.LRUs
    
    % .virtualChannel
    % Each entry in the virtual channel table has 
    %   bits(8:0) = frequency_id     (valid range 65 to 448)
    %   bits(12:9) = beam_id         (valid range 1 to 8)
    %   bits(15:13) = substation_id  (valid range 1 to 4)
    %   bits(20:16) = subarray_id    (valid range 1 to 16)
    %   bits(30:21) = station_id     (valid range 1 to 512)
    % There are two blocks of entries in the table, one for each of the two stations that come to this FPGA
    % Entries 0 to 383 relate to the first station, entries 512 to 895 are for the second station.
    LRUreg(LRU).virtualChannel = zeros(1024,1,'uint32');
    
    keyboard
    
end

