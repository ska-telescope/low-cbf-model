function [LRUreg, globalreg] = getRegisterSettings(arrayConfigFull,modelConfig)
%
% Generate the register settings in the firmware from the LMC configuration
% Returns:
%  LRUreg - an array of structs, where each struct has the registers for a single LRU ( = 1 FPGA).
%  globalReg - Register settings (or possibly ROM values) that are common to all FPGAs, e.g. filterbank taps.
%
% ------------------------------------------------------------------------------
% LRUreg
%
% There is a field for each register setting. These are named according the module and register name in the firmware.
% Also see the module level documentation on confluence.
% In some cases the register settings can change over time. To accomodate this, revised register settings are included
% in the output structure in the column index. e.g. localDoppler.startPhase refers to a block of 1024 registers. If
% there are 4 sets of values, then localDoppler.startPhase will be a 1024x4 matrix.
%
%  (1) LFAASPEADDecoder
%      - .virtualChannel 
%          * This is a lookup table with 384 entries used to assign the virtual channel number.
%  (2) localDoppler (LD)
%      Doppler corrections are derived from the delay polynomials. Multiple sets of register settings may be 
%      created for these settings. The time at which the registers switch over is stored in .validFrom
%      - .LDstationID0
%      - .LDstationID1
%      - .LDcountOffset0
%      - .LDcountOffset1
%      - .LDstartPhase
%      - .LDphaseStep
%      - .LDvalidFrom 
%         * Time in seconds at which each set of register settings becomes valid.
%           (not a register in the firmware)
%          
%  (3) coarseCornerTurner (CCT)
%      - .
%          * Sample delays for each station & channel
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
% ------------------------------------------------------------------------------
% globalreg
% 
%  (1) correlatorFilterbank
%      - .correlatorFilterbank  : FIR filter coefficients for the correlator filterbank
%  (2) PSSFilterbank
%      - .PSSFilterbank : FIR Filter Coefficients
%  (3) PSTFilterbank
%      - .PSTFilterbank : FIR Filter Coefficients
% 
%
%


%% Assign default values to the registers.
defaultReg.zxy = [1, 1, 1];
defaultReg.virtualChannel = zeros(1024,1,'uint32');
% LD : Registers for Local Doppler module
defaultReg.LDStationID0 = 0;
defaultReg.LDStationID1 = 0;
% Work out how many delay updates we need
if (modelConfig.delayUpdatePeriod == 0)
    updates = 1; % i.e. just the initial value.
else
    totalRuntime = modelConfig.runtime * 1080e-9 * 408 * 2048; 
    updates = floor(totalRuntime / modelConfig.delayUpdatePeriod) + 1;
end
defaultReg.LDcountOffset0 = zeros(1,updates);
defaultReg.LDcountOffset1 = zeros(1,updates);
defaultReg.LDstartPhase = zeros(1536,updates);
defaultReg.LDphaseStep = zeros(1536,updates);
defaultReg.LDvalidFrom = zeros(1,updates);
% CCT : Registers for Coarse Corner Turn module
defaultReg.CCTdelayTable = zeros(768,modelConfig.runtime,'int32');

LRUreg(1:modelConfig.LRUs) = defaultReg;

fullVirtualChannelTable = zeros(modelConfig.LRUs*1024,11,'int16');  % Use int16 to limit how big this gets.

%% Set register values.
for LRU = 1:modelConfig.LRUs
    
    % Not a register as such, but defines the location of the LRU in the array
    % The loop variable "LRU" first counts in the Z direction, then X, then Y.
    if (modelConfig.configuration == 0)  %  0 = PISA,     3 LRUs,  6 stations,   Z connect = 3, X = 1, Y = 1
        LRUreg(LRU).zxy = [LRU, 1, 1];
    elseif (modelConfig.configuration == 1) % 1 = AA1,     12 LRUs,  24 stations,  Z connect = 2, X = 6, Y = 1
        LRUreg(LRU).zxy = [mod(LRU,2) + 1, floor((LRU-1)/2) + 1, 1];
    elseif (modelConfig.configuration == 1) % 2 = AA2,     36 LRUs,  72 stations,  Z connect = 6, X = 6, Y = 1
        LRUreg(LRU).zxy = [mod(LRU-1,6) + 1, floor((LRU-1)/6) + 1, 1];
    elseif (modelConfig.configuration == 1) % 3 = AA3-ITF, 48 LRUs,  96 stations,  Z connect = 4, X = 6, Y = 2
        LRUreg(LRU).zxy = [mod(LRU-1,4) + 1, mod(floor((LRU-1)/4),24) + 1, floor((LRU-1)/24) + 1];
    elseif (modelConfig.configuration == 1) % 4 = AA3-CPF, 144 LRUs, 256 stations, Z connect = 8, X = 3, Y = 6
        LRUreg(LRU).zxy = [mod(LRU-1,8) + 1, mod(floor((LRU-1)/8),24) + 1, floor((LRU-1)/24) + 1];
    elseif (modelConfig.configuration == 1) % 5 = AA4,     288 LRUs, 512 stations, Z connect = 8, X = 6, Y = 6
        LRUreg(LRU).zxy = [mod(LRU-1,8) + 1, mod(floor((LRU-1)/8),48) + 1, floor((LRU-1)/48) + 1];
    end
    %% Registers for LFAASPEADDecoder module
    % .virtualChannel
    % Each entry in the virtual channel table has 
    %   bits(8:0) = frequency_id     (valid range 65 to 448)
    %   bits(12:9) = beam_id         (valid range 1 to 8)
    %   bits(15:13) = substation_id  (valid range 1 to 4)
    %   bits(20:16) = subarray_id    (valid range 1 to 16)
    %   bits(30:21) = station_id     (valid range 1 to 512)
    %   bit(31) = '1' to indicates that this entry is invalid.
    % There are two blocks of entries in the table, one for each of the two stations that come to this FPGA
    % Entries 0 to 383 relate to the first station, entries 512 to 895 are for the second station.
    
    LRUreg(LRU).virtualChannel = zeros(1024,1,'uint32');
    
    for station = 0:1 % Two stations per FPGA; first uses entries 0-383 in the table, second uses entries 512-895.
        % Get the station ID that goes to this FPGA
        if (((LRU-1)*2 + station + 1) <= length(modelConfig.stationMap))
            station_id = modelConfig.stationMap((LRU-1)*2 + station + 1);
            % Get all the logical IDs associated with this station ID
            logicalIDs = [];
            substationIDs = [];
            for stype = 1:length(arrayConfigFull.stations)
                IDindex = find(arrayConfigFull.stations(stype).stationIDs == station_id);
                logicalIDs = [logicalIDs arrayConfigFull.stations(stype).logicalIDs(IDindex)];
                substationIDs = [substationIDs arrayConfigFull.stations(stype).substationIDs(IDindex)];
            end
            
            % Search the subarrays for these logical IDs - Note each logical ID can only be a member of a single subArray
            subArrayIDs = zeros(length(logicalIDs),1);
            subArrayIndexes = zeros(length(logicalIDs),1);
            NsubArrays = length(arrayConfigFull.subArray);
            for logicalIDIndex = 1:length(logicalIDs)
                logicalID = logicalIDs(logicalIDIndex);
                for subArray = 1:NsubArrays
                    if any(arrayConfigFull.subArray(subArray).logicalIDs == logicalID)
                        subArrayIDs(logicalIDIndex) = arrayConfigFull.subArray(subArray).index;
                        subArrayIndexes(logicalIDIndex) = subArray;
                    end
                end
            end

            % At this point, we have :
            %  station_id - single value for this particular station
            %   logicalIDs      - List of logical IDs for this station_id. A single value for 300 MHz stations or up to 4 values when multiple substations are used.
            %   substationIDs   - substationIDs corresponding to the logicalIDs
            %   subArrayIDs     - subArrayIDs corresponding to the logical IDs
            %   subArrayIndexes - Index into arrayConfigFull.subArray
            %         
            % Find the configured station beams and frequency channels for each of the subArrays in the list subArrayIDs

            virtualChannelCount = 0;
            % Columns in virtualChannelTable
            %  1 => frequency_id, 2 => beam_id, 3 => substation_id, 4 => subarray_id, 5 => station_id, 
            %  6 => index X into arrayConfigFull.stationBeam(X), 7 => index X into arrayConfigFull.subArray(X), 8 => logical ID, 9 => table entry is unused
            % 10 => LRU this stationID came in on.
            % 11 => Virtual Channel.
            virtualChannelTable = zeros(512,11);

            for subArray = 1:length(subArrayIDs)
                currentSubArray = subArrayIDs(subArray);

                NstationBeams = length(arrayConfigFull.stationBeam);
                for b1 = 1:NstationBeams
                    if (arrayConfigFull.stationBeam(b1).subArray == currentSubArray)
                        % This stationBeam (arrayConfigFull.stationBeam(b1)) uses the subarray that this station is part of; 
                        % Add all the channels to the list for this station
                        for c1 = 1:length(arrayConfigFull.stationBeam(b1).channels)
                            virtualChannelCount = virtualChannelCount + 1;
                            % These entries get combined to form the virtual channel table in the LFAASPEADDecoder module
                            virtualChannelTable(virtualChannelCount,1) = arrayConfigFull.stationBeam(b1).channels(c1);
                            virtualChannelTable(virtualChannelCount,2) = arrayConfigFull.stationBeam(b1).index;
                            virtualChannelTable(virtualChannelCount,3) = substationIDs(subArray);
                            virtualChannelTable(virtualChannelCount,4) = currentSubArray;
                            virtualChannelTable(virtualChannelCount,5) = station_id;
                            % The remaining entries are for convenience use later in this function.
                            virtualChannelTable(virtualChannelCount,6) = b1;  % Just for convenience; used later in this function to look up the correct stationBeam parameters for this entry.
                            virtualChannelTable(virtualChannelCount,7) = subArrayIndexes(subArray);  % Also just for convenience; used to look up the subarray this virtual channel is part of.
                            virtualChannelTable(virtualChannelCount,8) = logicalIDs(subArray);  % logical ID
                            virtualChannelTable(virtualChannelCount,9) = 0; % this entry is not invalid.
                            virtualChannelTable(virtualChannelCount,10) = LRU;
                            virtualChannelTable(virtualChannelCount,11) = virtualChannelCount-1;
                        end
                    end
                end

            end

            % Convert virtualChannelTable matrix into the register settings
            virtualChannelTemp = zeros(512,1,'uint32');
            for vc = 1:virtualChannelCount
                virtualChannelTemp(vc) = uint32(virtualChannelTable(vc,1) + ...
                                                virtualChannelTable(vc,2) * 2^9 + ...
                                                virtualChannelTable(vc,3) * 2^13 + ...
                                                virtualChannelTable(vc,4) * 2^16 + ...
                                                virtualChannelTable(vc,5) * 2^21);
            end
            virtualChannelTemp((virtualChannelCount+1):512) = uint32(2^31);  % just set the top bit to indicate the entry is invalid.
            virtualChannelTable((virtualChannelCount+1):512,9) = 1; % This entry is invalid.
        else
            % No LFAA data into this FPGA.
            virtualChannelTemp = zeros(512,1,'uint32');
            virtualChannelTable = zeros(512,9);
            virtualChannelTemp(1:512) = uint32(2^31);  % Just set the top bit to indicate the entry is invalid.
            virtualChannelTable(1:512,9) = 1;          % All entries are invalid.
            station_id = -1;                           % No LFAA station sends data to this FPGA
        end
        
        if (station == 0)
            LRUreg(LRU).virtualChannel(1:512) = virtualChannelTemp;
            LRUreg(LRU).stationID(1) = station_id;       % record which station we expect to get data from.
            %LRUreg(LRU).virtualChannelTable0 = virtualChannelTable;  % store for easier reference in this function.
        else
            LRUreg(LRU).virtualChannel(513:1024) = virtualChannelTemp;
            LRUreg(LRU).stationID(2) = station_id;
            %LRUreg(LRU).virtualChannelTable1 = virtualChannelTable;  % store for easier reference in this function.
        end
        fullVirtualChannelTable(((LRU-1)*1024 + 512 * station + 1):((LRU-1)*1024 + 512 * station + 512),:) = int16(virtualChannelTable);
    end
    
    %% Registers for LocalDoppler module (LD)
    % .LDStationID0
    % .LDStationID1
    % .LDcountOffset0
    % .LDcountOffset1
    % .LDstartPhase
    % .LDphaseStep
    % .LDvalidFrom     
    
    LRUreg(LRU).LDStationID0 = modelConfig.stationMap((LRU-1)*2 + 1);
    LRUreg(LRU).LDStationID1 = modelConfig.stationMap((LRU-1)*2 + 2);
    
    %LRUreg(LRU).LDstartPhase
    %LRUreg(LRU).LDphaseIncrement
    % Step through the virtual channel tables and derive the start phase and phase increment from the delay Polynomials.
    for station = 0:1
        virtualChannelTable = double(fullVirtualChannelTable(((LRU-1)*1024 + 512 * station + 1):((LRU-1)*1024 + 512 * station + 512),:));
        for virtualChannel = 1:384
            
            if (virtualChannelTable(virtualChannel,9) == 0) % Not invalid, i.e. valid.
                logicalIDIndex = find(arrayConfigFull.subArray(virtualChannelTable(virtualChannel,7)).logicalIDs == virtualChannelTable(virtualChannel,8));
                delayPoly = arrayConfigFull.stationBeam(virtualChannelTable(virtualChannel,6)).delayPolynomial(logicalIDIndex,:);  % The delay polynomial for this virtual channel 
                HVOffset = arrayConfigFull.stationBeam(virtualChannelTable(virtualChannel,6)).HVOffset(logicalIDIndex);  % Offsets between the horizontal and vertical polarisations.
            else
                delayPoly = [0 0 0 0];  % Not a valid virtual channel, just set the polynomial to zeros.
                HVOffset = 0;
            end
            
            % Work out how many delay updates we need
            if (modelConfig.delayUpdatePeriod == 0)
                updates = 1; % i.e. just the initial value.
            else
                totalRuntime = modelConfig.runtime * 1080e-9 * 408 * 2048; 
                updates = floor(totalRuntime / modelConfig.delayUpdatePeriod) + 1;
            end
            
            % The slope at time t for the delay polynomial is -
            %   Delay = delayPoly(1) * t^3 + delayPoly(2) * t^2 + delayPoly(3) * t + delayPoly(4)
            %  where t is in seconds.
            %   So d(Delay)/dt = 3 * delayPoly(1) * t^2 + 2 * delayPoly(2) * t + delayPoly(3)
            %  e.g. at t = 0, d(Delay)/dt = delayPoly(3)
            % We also need the sky frequency from virtualChannelTable(virtualChannel,1)
            currentTime = 0;
            for update = 1:updates
                % Convert currentTime to sit at a whole number of LFAA blocks
                currentTimeRounded = floor(currentTime / (2048 * 1080e-9)) * 2048 * 1080e-9;
                % Generate a linear approximation to the delay polynomial 
                initialDelay = delayPoly(1) * currentTimeRounded^3 + delayPoly(2) * currentTimeRounded^2 + delayPoly(3) * currentTimeRounded + delayPoly(4);
                DdelayDt = 3 * delayPoly(1) * currentTimeRounded^2 + 2 * delayPoly(2) * currentTimeRounded + delayPoly(3); % delay in seconds/second
                skyFrequency = virtualChannelTable(virtualChannel,1) * 781.25e3;  % frequency in Hz
                skyPeriod = 1/skyFrequency;
                
                if (virtualChannelTable(virtualChannel,9) == 0) % Not invalid, i.e. valid.
                    LRUreg(LRU).LDcountOffset0(1,update) = round(currentTimeRounded / (2048 * 1080e-9));
                    LRUreg(LRU).LDcountOffset1(1,update) = round(currentTimeRounded / (2048 * 1080e-9));
                    
                    % Note (DdelayDt / skyPeriod) = Doppler revolutions per second.
                    % First 384 = horizontal polarisation
                    LRUreg(LRU).LDstartPhase(station*768 + virtualChannel,update) = round(2^32 * rem(initialDelay,skyPeriod)/skyPeriod);
                    LRUreg(LRU).LDphaseStep(station*768 + virtualChannel,update) = round((2^35)*(1080e-9)*(DdelayDt / skyPeriod));
                    % Second 384 = vertical polarisation
                    LRUreg(LRU).LDstartPhase(station*768 + 384 + virtualChannel,update) = round(2^32 * rem(initialDelay + HVOffset,skyPeriod)/skyPeriod);
                    LRUreg(LRU).LDphaseStep(station*768 + 384 + virtualChannel,update) = round((2^35)*(1080e-9)*(DdelayDt / skyPeriod));
                else % This virtual channel is not used.
                    LRUreg(LRU).LDcountOffset0(1,update) = 0;
                    LRUreg(LRU).LDcountOffset1(1,update) = 0;
                    LRUreg(LRU).LDstartPhase(station*768 + virtualChannel,update) = 0;
                    LRUreg(LRU).LDphaseStep(station*768 + virtualChannel,update) = 0;
                    LRUreg(LRU).LDstartPhase(station*768 + 384 + virtualChannel,update) = 0;
                    LRUreg(LRU).LDphaseStep(station*768 + 384 + virtualChannel,update) = 0;
                end
                LRUreg(LRU).LDvalidFrom(1,update) = currentTimeRounded;
                
                currentTime = currentTime + modelConfig.delayUpdatePeriod;
            end
            
            % Check the maximum size of the phase discontinuity at the update boundaries
            if (updates > 1)
                for update = 2:updates
                    predictedPhase = (2^(-32)) * LRUreg(LRU).LDstartPhase(station*768 + virtualChannel,update-1) + (2^(-35)) * 2048 * LRUreg(LRU).LDphaseStep(station*768 + virtualChannel,update) * (LRUreg(LRU).LDcountOffset0(1,update) - LRUreg(LRU).LDcountOffset0(1,update-1));
                    actualPhase = (2^(-32)) * LRUreg(LRU).LDstartPhase(station*768 + virtualChannel,update);
                    phaseError = predictedPhase - actualPhase - round(predictedPhase - actualPhase); % It is possible to be out by whole rotations due to dropping whole rotations from the 32 bit representation in the LDstartPhase register.
                    if (phaseError > 0.001)
                        disp(['Warning : Linear approximation phase error is ' num2str(phaseError) ' rotations. Use a shorter update period for the local Doppler parameters.']);
                        keyboard
                    end
                end
            end
            
        end
    end
end

%% Restart the loop through FPGAs, as we need all the setting for the input processing before dealing with the next stage, starting at the Coarse Corner Turn.

for LRU = 1:modelConfig.LRUs
    %% Registers for CoarseCornerTurner module (CCT)
    %  Delay for each coarse channel, based on the virtual channel number.
    %
    % For modelConfig.configuration of
    %  - 0 (PISA,     3 LRUs,  6 stations,   Z connect = 3) : Each corner turner has all 6 stations for 384/3 = 128 of the virtual channels. Virtual channels are in steps of 3 (i.e. 0, 3, 6... or 1, 4, 7... or 2, 5, 8...)
    %  - 1 (AA1,     12 LRUs,  24 stations,  Z connect = 2) : Each corner turner has 4 stations for 384/2 = 192 of the virtual channels. Virtual channels are in steps of 2.
    %  - 2 (AA2,     36 LRUs,  72 stations,  Z connect = 2) : Each corner turner has 4 stations for 384/2 = 192 of the virtual channels. Virtual channels are in steps of 2.
    %  - 3 (AA3-ITF, 48 LRUs,  96 stations,  Z connect = 4) : Each corner turner has 8 stations for 384/4 = 96 of the virtual channels. Virtual channels are in steps of 4.
    %  - 4 (AA3-CPF, 144 LRUs, 256 stations, Z connect = 8) : Each corner turner has 16 stations for 384/8 = 48 of the virtual channels. Virtual channels are in steps of 8.
    %  - 5 (AA4,     288 LRUs, 512 stations, Z connect = 8) : Each corner turner has 16 stations for 384/8 = 48 of the virtual channels. Virtual channels are in steps of 8.
    % 
    % For all cases of modelConfig.configuration, there are 768 distinct data sets that the coarse corner turn deals with.
    % Registers for each LRU:
    %  - .CCTdelayTable : Table with 1536 entries, coarse delay and fine delay parameters for each station/channel combination.
    %     There are two copies of this table in the firmware to enable updating the table without generating glitches. In the model,
    %     this is a matrix with one column for each 0.9 second integration period.
    %  - .CCTstationList : Table listing all the stations used in this LRU (This will be the same for all LRUs with the same X,Y coordinates).
    %  - .CCTvirtualChannelOffset : The virtual channel offset is the first virtual channel sent to this FPGA. This is defined by the Z coordinate.
    
    % First, get stationList, the stations that are routed to this coarseCornerTurn module.
    stationList = [];
    for LRU2 = 1:modelConfig.LRUs
        if (LRUreg(LRU2).zxy(2:3) == LRUreg(LRU).zxy(2:3))  % All LRUs on the same Z connect
            stationList = [stationList LRUreg(LRU2).stationID];
        end
    end
    LRUreg(LRU).CCTstationList = stationList;
    LRUreg(LRU).virtualChannelOffset = LRUreg(LRU).zxy(1) - 1;
    if (modelConfig.configuration == 0)
        virtualChannelStepsize = 3;  % For PISA, every 3rd virtual channel goes to a given LRU
    elseif (modelConfig.configuration == 1)
        virtualChannelStepsize = 2;
    elseif (modelConfig.configuration == 2)
        virtualChannelStepsize = 6;
    elseif (modelConfig.configuration == 3)
        virtualChannelStepsize = 4;
    else
        virtualChannelStepsize = 8;
    end
    
    % Now get the actual coarse delays
    if isempty(LRUreg(LRU).CCTstationList)
        LRUreg(LRU).CCTenable = 0;
        LRUreg(LRU).CCTcoarseDelay = zeros(768,modelConfig.runtime,'uint32');  % Note modelConfig.runtime is the number of 0.9 second integration periods.
    else
        LRUreg(LRU).CCTenable = 1;
        LRUreg(LRU).CCTcoarseDelay = zeros(768,modelConfig.runtime,'uint32');
        allVirtualChannels = LRUreg(LRU).virtualChannelOffset:virtualChannelStepsize:383;
        for virtualChannelCount = 1:128
            virtualChannel = allVirtualChannels(virtualChannelCount);
            for stationIndex = 1:length(stationList)
                stationID = stationList(stationIndex);
                % Look up this stationID and virtual Channel in fullVirtualChannelTable
                %keyboard
                fullIndex = find((fullVirtualChannelTable(:,5) == stationID) & (fullVirtualChannelTable(:,11) == virtualChannel));
                % Get the delay polynomial for this station and beam
                if (fullVirtualChannelTable(fullIndex,9) == 0) % Not invalid, i.e. valid.
                    logicalIDIndex = find(arrayConfigFull.subArray(fullVirtualChannelTable(fullIndex,7)).logicalIDs == fullVirtualChannelTable(fullIndex,8));
                    delayPoly = arrayConfigFull.stationBeam(fullVirtualChannelTable(fullIndex,6)).delayPolynomial(logicalIDIndex,:);  % The delay polynomial for this virtual channel 
                    HVOffset = arrayConfigFull.stationBeam(fullVirtualChannelTable(fullIndex,6)).HVOffset(logicalIDIndex);
                else
                    delayPoly = [0 0 0 0];  % Not a valid virtual channel, just set the polynomial to zeros.
                    HVOffset = 0;
                end
                for frame = 1:modelConfig.runtime
                    currentTime = (frame - 1) * 1080e-9 * 2048 * 408;
                    % Evaluate the delay polynomial and its derivative
                    delaySeconds = delayPoly(1) * currentTime^3 + delayPoly(2) * currentTime^2 + delayPoly(3) * currentTime + delayPoly(4);
                    DdelayDt = 3 * delayPoly(1) * currentTime^2 + 2 * delayPoly(2) * currentTime + delayPoly(3);
                    delaySamples = round(delaySeconds/(1080e-9));
                    delayFracH = delaySeconds - delaySamples * 1080e-9; % in seconds
                    delayFracV = delayFracH + HVOffset;
                    deltaPH = round((1/1080e-9) * (1/2) * delayFracH/(2^(-12)) * 2^3);  % See the coarse corner turn confluence page for an explanation.
                    deltaPV = round((1/1080e-9) * (1/2) * delayFracV/(2^(-12)) * 2^3);
                    deltadeltaP = round((1/1080e-9) * (1/2) * (DdelayDt * 64*1080e-9)/(2^(-15)) * 2^15);
                    if ((delaySamples < 0) || (delaySamples > 4095))
                        error(['Calculated Coarse Delay is ' num2str(delaySamples) ' LFAA samples. It must be between 0 and 4095']);
                    end
                    % Two words in the table for each channel and station.
                    LRUreg(LRU).CCTdelayTable(2*((virtualChannelCount-1)*6 + (stationIndex-1)) + 1,frame) = 2^16 * deltaPH + delaySamples;
                    LRUreg(LRU).CCTdelayTable(2*((virtualChannelCount-1)*6 + (stationIndex-1)) + 2,frame) = 2^16 * deltadeltaP + deltaPV;
                end
            end
        end
    end
    
end

%% Registers for the filterbanks.
% The filter taps for the filterbanks are fixed for all FPGAs, so they are stored in globalreg instead of LRUreg.
% generate_MaxFlt creates filters with a DC response of 0.5
% i.e. sum(globalreg.PSSFilterbankTaps(X:64:end)) = 0.5
% The maximum value for the filter taps is about 0.55
% The firmware has 18 bits to represent the filter taps (i.e. range = -131072 to 131071).
% So scale by 2^17 (=131072)
globalreg.PSSFilterbankTaps = round(2^17 * generate_MaxFlt(64,12));   % PSS, 64 point FFT, 12 taps.
globalreg.PSTFilterbankTaps = round(2^17 * generate_MaxFlt(256,12));  % PST, 256 point FFT, 12 taps
globalreg.correlatorFilterbankTaps = round(2^17 * generate_MaxFlt(4096,12)); % Correlator, 4096 point FFT, 12 taps.

%% Registers for PSS beamforming
% Each PSS beam has a set of complex weights for the stations involved.
% up to 512 stations contribute to each beam. 
% For each beam we need a complex weight for each station and each fine channel.
% So the maximum is (500 beams) * (512 stations) * (54*384 fine channels) * (2bytes/weight (at least)) / (288 LRUs) = 36 Mbytes of weights per FPGA
% Obviously that won't work.
% Instead, we have a phase slope across the entire (300MHz) band - same data needs to be stored for all LRUs, so
% (500 beams) * (512 stations) * (4 bytes/weight at least - 2 for slope + 2 for linear change in slope) = 1 Mbyte = 32 URAMs.
%  Data required for each station and beam :
%   * Phase slope
%   * Linear change in the phase slope
%   * valid - i.e. If this station contributes to this beam or not.
% We also need to specify which stations and coarse channels contribute to each beam - do this with a bit vector for each beam and each LRU
% For AA4, each LRU has 72 fine channels from all 512 stations (9 fine channels from each of 8 coarse channels).
%
% Registers (Same for all LRUs)
%  globalreg.PSSPhase : array 512x500xT (station x beam x updates), with phase across a single 14.49kHz channel, in units of 2^(-25) rotations. Signed 24 bit value.
%  globalreg.PSSPhaseDelta : array 512x500xT (station x beam x updates), with the phase change delta per 64 LFAA samples, in units of 2^(-42) rotations. Signed 24 bit value.
%  globalreg.PSSUsed : array 512x500 (station x beam), indicates if the station is used for the beam.
globalreg.PSSPhase = zeros(512,500,modelConfig.runtime,'int32');
globalreg.PSSPhaseDelta = zeros(512,500,modelConfig.runtime,'int32');
globalreg.PSSUsed = zeros(512,500,'int32');
% Get a list of stations (indexed by stationIDs, not logicalIDs)
stationList = unique(arrayConfigFull.stationsFull.stationIDs);
for stationIndex = 1:length(stationList)
    stationID = stationList(stationIndex);
    
    % Step through the PSS beams
    NBeams = length(arrayConfigFull.PSSBeam);
    for beam = 1:NBeams
        % -------------
        % Work out if this beam uses this stationID, and if so, which polynomial it is.
        % Get the list of logical IDs for the stations in the subarray for this beam.
        logicalIDs = arrayConfigFull.subArray(arrayConfigFull.PSSBeam(beam).subArray).logicalIDs;
        % Check if any of these logicalIDs match this stationID
        matchesFound = 0;
        for i1 = 1:length(logicalIDs)
            if (stationID == arrayConfigFull.stationsFull.stationIDs(find(arrayConfigFull.stationsFull.logicalIDs == logicalIDs(i1))))
                matchesFound = matchesFound + 1;
                matchIndex = i1;
            end
        end
        if (matchesFound > 1)
            error('More than one matching logical ID for the station ID in PSS beamforming');
        end
        % -------------
        if (matchesFound == 1)
            % Get the Polynomial
            delayPoly = arrayConfigFull.PSSBeam(beam).delayPolynomial(matchIndex,:);
            % Evaluate the polynomial and its derivative at the PSS update interval
            % The PSS update interval is currently set here to once per 0.9 second frame - in the real system this is more likely to be about once per 10 seconds.
            for frame = 1:modelConfig.runtime
                currentTime = (frame - 1) * 1080e-9 * 2048 * 408;
                delaySeconds = delayPoly(1) * currentTime^3 + delayPoly(2) * currentTime^2 + delayPoly(3) * currentTime + delayPoly(4);
                DdelayDt = 3 * delayPoly(1) * currentTime^2 + 2 * delayPoly(2) * currentTime + delayPoly(3);
                phaseRegister = ((1/1080e-9)/64) * delaySeconds / (2^(-25));   % Phase across one PSS fine channel, in units of 2^(-25) rotations.
                deltaPhaseRegister = ((1/1080e-9)/64) * DdelayDt * 64 * 1080e-9 / (2^(-42)); % phase step, units of 2^(-42) rotations per 64 LFAA samples.
                % Check values are reasonable
                if ((abs(phaseRegister) > 2^23) || (abs(deltaPhaseRegister) > 2^23))
                    warning('PSS phase corrections do not fit in 3 bytes');
                end
                globalreg.PSSPhase(stationID,beam,frame) = round(phaseRegister);
                globalreg.PSSPhaseDelta(stationID,beam,frame) = round(deltaPhaseRegister);
                globalreg.PSSUsed(stationID,beam) = 1;
                % Check the error in the linear approximation from the first to the last frame
                if (frame > 1)
                    predictedPhase = double(globalreg.PSSPhase(stationID,beam,1)) + 2^(-17) * double(globalreg.PSSPhaseDelta(stationID,beam,1)) * (frame - 1) * 2048 * 408 / 64; 
                    actualPhase = globalreg.PSSPhase(stationID,beam,frame);
                    if ((predictedPhase - actualPhase) > 2)
                        warning('Linear prediction for PSS phase is poor');
                        keyboard
                    end
                end
                %keyboard
            end
        else
            % No match found, i.e. this station does not contribute to this beam.
            globalreg.PSSPhase(stationID,beam,:) = 0;
            globalreg.PSSPhaseDelta(stationID,beam,:) = 0;
            globalreg.PSSUsed(stationID,beam) = 0;
        end
    end
    
end

%% ----------------------------------------------------------------------------------
% Supporting functions
%
function W = generate_MaxFlt(nbuff, nTap)
% Generate maximal flat filter coefficients
% Adapted from John Buntons code.

%{
% Author: J Bunton, 22 August 2015
Filter Response to meet correlator requirements
For a monochromatic signal, total power (all channels) remains constant
independent of frequency
Starting point are the maximally flat filters
this is improved with some simple optimisation

Calculation should be done only once per Simulink simulation  

nbuff = typically 4096
nTaps = # of taps, typically 8 or 12 

%} 

%{ 
SVN keywords
 $Rev::                                                                                            $: Revision of last commit
 $Author::                                                                                         $: Author of last commit
 $Date::                                                                                           $: Date of last commit
 $LastChangedDate::                                                                                $: Date of last change
 $HeadURL: svn://codehostingbt.aut.ac.nz/svn/LOWCBF/Modelling/CSP_DSP/CSP_Dataflow/generate_MaxFl#$: Repo location
%}


%% Filter design coefficients 
% disp('generate_MaxFlt')

nTap2 = 2*nTap;  % say around 8 or 12 
nTap2p1 = nTap2+1; 

imp=maxflat(nTap2,'sym',.5*nTap2/nTap2p1);
imp=interpft(imp,nTap2)*nTap2p1/nTap2; %Take to 2*ntap (24) tap filter (2 channel, 12tap FIR)

% plot(db(fft(imp)),'o-')

% Interate to improve (hard coded 10 times) 
for k=1:10

    impf=fft(imp);
    imph=imp.*cos( ((1:length(imp))-1)*pi);
    impfh=fft(imph);
    errorf =(impf.*conj(impf)+impfh.*conj(impfh));
    errorf=errorf/errorf(1);
    errorf=1-errorf;
    error=fftshift((ifft(errorf)));
    imp=imp+error/2.0; %(2.5*( abs(impf)+abs(impfh) ));

end

cor=imp;
corh=cor.*cos( ((1:length(imp))-1)*pi);

corf=freqz(cor,2048)*2048;
corfh=freqz(corh,2000)*2000;
ampf=corf.*conj(corf)+corfh.*conj(corfh);
error = fftshift((ifft(1-ampf)));


%{
Variable nbuff*nTap is the length of the filter. 
Typically: 4096 freq channels * 8 taps) 

12x32 is a 12 tap FIR section, 32 channel filterbank
%} 
%
W=interpft(cor,nbuff*nTap);  %change this line to alter length of filter.
W = W(:); % force column 

return 


