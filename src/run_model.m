function run_model(rundir)
    % Run the low.CBF packetised model
    % 
    % ---------------------------------------------------------------------
    % Inputs
    %  The configuration for the run is read from the run directory provided -
    %   'lmcConfig.txt' - json file with the arrayConfigFull structure, as created by create_config.m
    %   'registerSettings.txt' - json file with register settings for the run, as created by create_config.m
    %   'modelConfig.txt' - User defined json file describing the way the model is created, 
    %     e.g. which array configuration to use, amount of time to run for, which data products to keep.
    %
    % ----------------------------------------------------------------------
    % General structure of the code is similar to the firmware. 
    % The processing is split into stages based on when data is put on the interconnect between the FPGAs.
    % 
    % 1. Stage1 processing:
    %     For each FPGA
    %        load the input data from LFAA.mat
    %        Processing:
    %           * LFAA SPEAD ingest
    %           * Local Doppler
    %           * Coarse splitter
    %        save the results to stage1.mat
    %
    % 2. Stage2 Processing:
    %     For each FPGA
    %        Take input data from stage1
    %        Processing:
    %           * Coarse corner turn
    %           * filterbanks
    %           * fine delays
    %           * RFI
    %        save results to stage2.mat
    %
    % 3. Stage3 Processing:
    %     For each FPGA
    %        Take input data from stage2
    %        Processing:
    %          * Fine corner turn
    %          * Correlator
    %          * PSS & PST Beamformers
    %        save results to stage3.mat
    %        
    % 4. Stage4 Processing
    %     For each FPGA
    %        Take input data from stage3
    %        Processing:
    %          * PSS corner turn
    %          * packetising for SDP, PSS and PST
    % -----------------------------------------------------------------------

    %% read the configuration files
    % Model setup
    modelConfig = parseModelConfig(rundir);

    % Firmware register settings
    fid = fopen([rundir '/registerSettings.txt']);
    regjson = fread(fid,inf);
    fclose(fid);
    regjson = char(regjson');
    registers = jsondecode(regjson);

    % lmc Configuration
    fid = fopen([rundir '/lmcConfig.txt']);
    lmcjson = fread(fid,inf);
    fclose(fid);
    lmcjson = char(lmcjson');
    arrayConfigFull = jsondecode(lmcjson);

    %% -----------------------------------------------------------------------
    % Stage1 Processing - LFAA ingest and Local Doppler

    % Load the data file
    fname = [rundir '/LFAA'];
    if exist(fullfile(cd,rundir,'LFAA.mat'),'file')
        load(fname);  % Should load a variable "fpga"
    else
        error(['Cannot find ' fname '. Run create_config.m to generate it']);
    end

    for lru = 1:modelConfig.LRUs

        % Step through the packets, gather statistics, and form the station processing output packets.
        % The only processing is the Doppler correction.
        % The output of the station processing contains the same amount of data as the input, except that
        % the header on each packet is much smaller.

        % First, generate the output data structure
        % Output data has fields
        %  .headers : 4 x <number of input packets>. Four things in the header are as defined on the confluence page -
        %             - virtual channel
        %             - channel frequency
        %             - station ID
        %             - packet count
        %  .data    : NxM, as per the input LFAA data. N is the number of data samples, M is the number of non-zero channels.
        %  .dataPointers : <number of packets> x 2. As per fpga(lru).dataPointers. 
        %                  First column is the row index into .data,
        %                  Second column is the column index into .data. 
        %                  Note this is identical to the input structure, fpga(lru).dataPointers

        hsize = size(fpga(lru).headers); % hsize(1) should be 114 (number of bytes in an LFAA header, hsize(2) is the number of headers.
        dsize = size(fpga(lru).data);
        NChannels = dsize(2); % Number of non-zero channels.

        stage1(lru).headers = zeros(4,hsize(2));
        stage1(lru).data = fpga(lru).data;  % start with a copy of the input data, then apply the Doppler correction
        channelDone = zeros(NChannels,1);  % Used to keep track of when we have processed the Doppler correction for each channel
        stage1(lru).dataPointers = fpga(lru).dataPointers;

        for packet = 1:hsize(2)
            fullHeader = fpga(lru).headers(:,packet);
            % Get the relevant fields from the header and construct the header used on the internal packet.
            % bytes 1:12 - MAC addresses, bytes 13:14 - ethertype, bytes 15:34 - IPv4 Header, bytes 35:42 - UDP header
            % bytes 43:114 - SPEAD header
            SPEADHeader = fullHeader(43:114);
            % Extract the relevant fields from the SPEAD header
            frequency_id = double(SPEADHeader(55)) * 256 + double(SPEADHeader(56));
            beam_id = double(SPEADHeader(53)) * 256 + double(SPEADHeader(54));
            substation_id = double(SPEADHeader(59));
            subarray_id = double(SPEADHeader(60));
            station_id = double(SPEADHeader(61)) * 256 + double(SPEADHeader(62));
            packet_count = double(SPEADHeader(13)) * 2^24 + double(SPEADHeader(14)) * 2^16 + double(SPEADHeader(15)) * 256 + double(SPEADHeader(16));
            % Compile into a 32 bit number to compare with the entries in the virtual channel table
            virtualChannelContent = station_id * 2^21 + subarray_id * 2^16 + substation_id * 2^13 + beam_id * 2^9 + frequency_id;

            if (station_id == registers.LRU(lru).stationID(1))
                virtualChannel = find(registers.LRU(lru).virtualChannel(1:384) == virtualChannelContent);
            elseif (station_id == registers.LRU(lru).stationID(2))
                virtualChannel = find(registers.LRU(lru).virtualChannel(513:896) == virtualChannelContent);
            else
                error('No match for the stationID in the virtual channel table');
            end

            if (length(virtualChannel) ~= 1)
                error('No Match or more than 1 match in the virtual channel table');
            end
            virtualChannel = virtualChannel - 1;

            % Generate the internal packet header
            stage1(lru).headers(1,packet) = virtualChannel;
            stage1(lru).headers(2,packet) = frequency_id;
            stage1(lru).headers(3,packet) = station_id;
            stage1(lru).headers(4,packet) = packet_count;


            % Only process the data part of the packet if the data is non-zero
            if (fpga(lru).dataPointers(packet,2) ~= 0)
                % To be more efficient, this processes all the data for the channel in one go,
                % rather than as individual packets.
                % Check if we have already done this channel, and if not, do it.
                if ~channelDone(fpga(lru).dataPointers(packet,2))
                    % Do the Doppler correction for the entire channel
                    if (registers.LRU(lru).LDStationID0 == stage1(lru).headers(3,packet))
                        countOffset = registers.LRU(lru).LDcountOffset0;
                        startPhaseH = registers.LRU(lru).LDstartPhase(stage1(lru).headers(1,packet) + 1,:);
                        startPhaseV = registers.LRU(lru).LDstartPhase(stage1(lru).headers(1,packet) + 1 + 384,:);
                        phaseStepH = registers.LRU(lru).LDphaseStep(stage1(lru).headers(1,packet) + 1,:);
                        phaseStepV = registers.LRU(lru).LDphaseStep(stage1(lru).headers(1,packet) + 1 + 384,:);
                    elseif (registers.LRU(lru).LDStationID1 == stage1(lru).headers(3,packet))
                        countOffset = registers.LRU(lru).LDcountOffset1;
                        startPhaseH = registers.LRU(lru).LDstartPhase(stage1(lru).headers(1,packet) + 1 + 768,:);
                        startPhaseV = registers.LRU(lru).LDstartPhase(stage1(lru).headers(1,packet) + 1 + 1152,:);
                        phaseStepH = registers.LRU(lru).LDphaseStep(stage1(lru).headers(1,packet) + 1 + 768,:);
                        phaseStepV = registers.LRU(lru).LDphaseStep(stage1(lru).headers(1,packet) + 1 + 1152,:);
                    else
                        error('No match found for the Station ID in local Doppler');
                    end

                    % Registers have updates for each frame; Step through each frame and process the data
                    for frame = 1:modelConfig.runtime
                        % The data for this channel is stage1.data(:,stage1.dataPointers(packet,2)). Samples alternate between vertical and horizontal polarisations.

                        %for sample = 1:(408 * 2048)  % 408 blocks of 2048 samples per 0.9 second frame.
                        %    phaseH = startPhaseH(frame) * 2^(-32) + (sample-1) * phaseStepH(frame) * 2^(-35);
                        %    phaseV = startPhaseV(frame) * 2^(-32) + (sample-1) * phaseStepV(frame) * 2^(-35);
                        %    stage1(lru).data((frame-1)*408*2048*2 + 2*(sample-1) + 1, stage1(lru).dataPointers(packet,2)) = double(stage1(lru).data((frame-1)*408*2048*2 + 2*(sample-1) + 1, stage1(lru).dataPointers(packet,2))) * exp(1i*2*pi*phaseV);
                        %    stage1(lru).data((frame-1)*408*2048*2 + 2*(sample-1) + 2, stage1(lru).dataPointers(packet,2)) = double(stage1(lru).data((frame-1)*408*2048*2 + 2*(sample-1) + 2, stage1(lru).dataPointers(packet,2))) * exp(1i*2*pi*phaseH);
                        %end
                        % Vectorised version of the above code (about 10x faster) :
                        for packet1 = 1:408  % 408 blocks of 2048 samples per 0.9 second frame.
                            phaseH = shiftdim(startPhaseH(frame) * 2^(-32) + (packet1 - 1) * 2048 * phaseStepH(frame) * 2^(-35) + (0:2047) * phaseStepH(frame) * 2^(-35));
                            phaseV = shiftdim(startPhaseV(frame) * 2^(-32) + (packet1 - 1) * 2048 * phaseStepV(frame) * 2^(-35) + (0:2047) * phaseStepV(frame) * 2^(-35));
                            stage1(lru).data((frame-1)*408*2048*2 + (packet1-1)*4096 + (0:2:4095) + 1, stage1(lru).dataPointers(packet,2)) = double(stage1(lru).data((frame-1)*408*2048*2 + (packet1-1)*4096 + (0:2:4095) + 1, stage1(lru).dataPointers(packet,2))) .* exp(1i*2*pi*phaseV);
                            stage1(lru).data((frame-1)*408*2048*2 + (packet1-1)*4096 + (0:2:4095) + 2, stage1(lru).dataPointers(packet,2)) = double(stage1(lru).data((frame-1)*408*2048*2 + (packet1-1)*4096 + (0:2:4095) + 2, stage1(lru).dataPointers(packet,2))) .* exp(1i*2*pi*phaseH);
                        end
                    end
                    %keyboard
                    % Mark this channel as processed
                    channelDone(fpga(lru).dataPointers(packet,2)) = 1;
                end
            end
        end

    end

    % Save the results to stage1.mat
    save([rundir '/stage1'], 'stage1');

    % Gather packets from all LRUs that go to a particular destination
    keyboard

    %% -----------------------------------------------------------------------
    % Stage2 Processing 
    %  - coarse corner turn, filterbanks, fine delay
    for lru = 1:modelConfig.LRUs
        
        % Get the combinations of virtual channels and stations that this LRU processes
        [cz,cx,cy] = getLRUCoordinates(modelConfig.configuration,lru); % coordinates of this LRU in the array
        if (modelConfig.configuration == 0) % PISA, 3x1, all 6 stations in each LRU
            virtualChannels = ((cz-1):3:383);
            stations = modelConfig.stationMap;
            CCTstationsPerLRU = 6;  % Number of stations per LRU for the coarse corner turn. This is used for addressing and indexing. Not all stations have to be used (i.e. there may be no data from some stations).
            zStart = 1;  % index of the first LRU with the same z coordinate as this LRU
            corSplit = 1;  % number of ways each coarse channel is split at the output of the correlator.
        elseif (modelConfig.configuration == 1) % AA1, 2x6x1, 4 stations/LRU
            virtualChannels = ((cz-1):2:383);
            stations = modelConfig.stationMap((floor((lru-1)/2) * 4 + 1):(floor((lru-1)/2) * 4 + 4));
            CCTstationsPerLRU = 4;
            zStart = floor((lru-1)/2) * 2 + 1;
            corSplit = 6;
        elseif (modelConfig.configuration == 2) % AA2, 2x6x3, 4 stations/LRU
            virtualChannels = ((cz-1):2:383);
            stations = modelConfig.stationMap((floor((lru-1)/2) * 4 + 1):(floor((lru-1)/2) * 4 + 4));
            CCTstationsPerLRU = 4;
            zStart = floor((lru-1)/2) * 2 + 1;
            corSplit = 18;
        elseif (modelConfig.configuration == 3) % AA3-ITF, 4x6x2, 8 stations/LRU
            virtualChannels = ((cz-1):4:383);
            stations = modelConfig.stationMap((floor((lru-1)/4) * 8 + 1):(floor((lru-1)/4) * 8 + 8));
            CCTstationsPerLRU = 8;
            zStart = floor((lru-1)/4) * 4 + 1;
            corSplit = 12;
        elseif (modelConfig.configuration == 4) % AA3-CPF, 8x6x3, 16 stations/LRU
            virtualChannels = ((cz-1):8:383);
            stations = modelConfig.stationMap((floor((lru-1)/8) * 16 + 1):(floor((lru-1)/8) * 16 + 16));
            CCTstationsPerLRU = 16;
            zStart = floor((lru-1)/8) * 8 + 1;
            corSplit = 18;
        else % AA4, 8x6x6, 16 stations/LRU
            virtualChannels = ((cz-1):8:383);
            stations = modelConfig.stationMap((floor((lru-1)/8) * 16 + 1):(floor((lru-1)/8) * 16 + 16));
            CCTstationsPerLRU = 16;
            zStart = floor((lru-1)/8) * 8 + 1;
            corSplit = 36;
        end

        % Coarse corner turn inputs and outputs - just collect the appropriate packets and save to a data structure.
        if (modelConfig.keepCoarseCornerTurner)
            % TBD 
        end

        % The stage2 output data structure for the correlator only contains data for non-zero packets.
        % (unlike the stage1 structure which contains headers for all packets - stage2 
        %  outputs can be much smaller packets and so would use a lot of memory).
        % Find the number of active stations & channels
        totalStationsChannels = 0;
        for channelCount = 1:length(virtualChannels)
            channel = virtualChannels(channelCount);
            for stationCount = 1:length(stations)
                station = stations(stationCount);
                packetIndex = find((channel == stage1(srcLRU).headers(1,:)) & (station == stage1(srcLRU).headers(3,:)),1);  % First packet for this station and channel. Note that all packets for this station+channel should reference the same data stream.
                if (stage1(srcLRU).dataPointers(packetIndex,2) ~= 0)
                    totalStationsChannels = totalStationsChannels + 1;
                end
            end
        end
        stage2Cor(lru).headers = zeros(4,corSplit * totalStationsChannels * modelConfig.runtime * 204);
        stage2Cor(lru).data = zeros(2*3456/corSplit,corSplit * totalStationsChannels * modelConfig.runtime * 204);
        corPacketCount = 0; % incremented as packets are put into the stage2Cor structure
        
        % Filterbanks
        % Step through the channels and stations that this LRU processes
        for channelCount = 1:length(virtualChannels)
            channel = virtualChannels(channelCount);
            for stationCount = 1:length(stations)
                station = stations(stationCount);
                CCToffset = (channelCount-1) * CCTstationsPerLRU + (stationCount-1) + 1;
                % Work out which LRU this combination of station and channel comes from.
                % Note that the order we step through stations & channels here is also the order in which they are stored in the coarse corner turn buffer,
                % and the order in which the delay information is stored in the registers (i.e. in registers.LRU(lru).CCTdelayTable())
                % Z coordinate is just stationCount/2, since "stations" is a list of all the stations in that go to this z coordinate (and 2 stations/LRU).
                srcLRU = floor((stationCount-1)/2) + (zStart-1) + 1;
                packetIndex = find((channel == stage1(srcLRU).headers(1,:)) & (station == stage1(srcLRU).headers(3,:)),1);  % First packet for this station and channel. Note that all packets for this station+channel should reference the same data stream.
                if (stage1(srcLRU).dataPointers(packetIndex,2) ~= 0)
                    % Data exists for this channel
                    % Process it in 0.9 second frames, taking into account the coarse delay, and using zeros for initialisation data for the first frame.
                    % data is stage1(srcLRU).data(:,stage1(srcLRU).dataPointers(packetIndex,2))
                    % Coarse delay is in the low 16 bits of registers.LRU(lru).CCTdelayTable((CCToffset-1)*2 + 1)
                    keyboard
                    for frame = 1:modelConfig.runtime
                        coarseDelay = mod(registers.LRU(lru).CCTdelayTable((CCToffset-1)*2 + 1,frame),65536);
                        hPolDeltaPCoarse = floor(registers.LRU(lru).CCTdelayTable((CCToffset-1)*2 + 1,frame)/65536);
                        vPolDeltaPCoarse = mod(registers.LRU(lru).CCTdelayTable((CCToffset-1)*2 + 2,frame),65536);
                        deltaDeltaPCoarse = floor(registers.LRU(lru).CCTdelayTable((CCToffset-1)*2 + 2,frame)/65536);
                        % Correlator filterbank. Data per frame is 204 * 4096 samples, plus 11 * 4096 samples for the initialisation data.
                        if (frame == 1) % first frame, initialization for the filterbank is zeros.
                            vPolData = zeros((204+12)*4096,1);
                            hPolData = zeros((204+12)*4096,1);
                            vPolData((12*4096+1):end) = stage1(srcLRU).data(1:2:(2*204*4096),stage1(srcLRU).dataPointers(packetIndex,2));
                            hPolData((12*4096+1):end) = stage1(srcLRU).data(2:2:(2*204*4096),stage1(srcLRU).dataPointers(packetIndex,2));
                        else  % Part way through the frame, so we have initialisation data.
                            sampleStart = (frame-1) * 4096 * 204 * 2 - 12 * 4096 * 2;   % x2 because polarisations are interleaved.
                            vPolData = stage1(srcLRU).data((sampleStart+1):2:(sampleStart + 2*(204+12)*4096),stage1(srcLRU).dataPointers(packetIndex,2));
                            hPolData = stage1(srcLRU).data((sampleStart+2):2:(sampleStart + 2*(204+12)*4096),stage1(srcLRU).dataPointers(packetIndex,2));
                        end
                        vPolCorrelatorFB = filterbank(vPolData((coarseDelay+1):(coarseDelay+215*4096)),registers.global.correlatorFilterbankTaps,4096,1,3456);
                        hPolCorrelatorFB = filterbank(hPolData((coarseDelay+1):(coarseDelay+215*4096)),registers.global.correlatorFilterbankTaps,4096,1,3456);
                        % Fine time delay.
                        % The fine delay in vPolDeltaPCoarse is the fine delay at the start of the frame.
                        % We need to adjust it by deltaDeltaPCoarse to the group of 64 samples closest to the peak of the filter.
                        filterpeakStart = -12 * 4096 + coarseDelay + registers.global.correlatorFilterbankMaxIndex;
                        vPolCorrelatorFD = zeros(3456,204);
                        hPolCorrelatorFD = zeros(3456,204);
                        for block = 1:204
                            filterpeak = round((filterpeakStart + (block - 1)*4096)/64); % Number of 64 sample steps to where the filter peak is.
                            vPolDeltaPCoarseFinal = vPolDeltaPCoarse + filterpeak * deltaDeltaPCoarse * 2^(-15);
                            hPolDeltaPCoarseFinal = hPolDeltaPCoarse + filterpeak * deltaDeltaPCoarse * 2^(-15);
                            vPolPhaseCorrection = exp(1i * 2^(-15) * 2 * pi * ((-1728:1727).') * vPolDeltaPCoarseFinal/2048);
                            hPolPhaseCorrection = exp(1i * 2^(-15) * 2 * pi * ((-1728:1727).') * hPolDeltaPCoarseFinal/2048);
                            vPolCorrelatorFD(:,block) = vPolCorrelatorFB(:,block) .* vPolPhaseCorrection;
                            hPolCorrelatorFD(:,block) = hPolCorrelatorFB(:,block) .* hPolPhaseCorrection;
                        end
                        
                        keyboard
                        corSplitSize = 3456/corSplit;
                        for block = 1:204
                            for split = 1:corSplit
                                corPacketCount = corPacketCount + 1;
                                
%                                 stage2Cor(lru).headers(1,corPacketCount) = ;   % virtual channel (9 bits)
%                                 stage2Cor(lru).headers(2,corPacketCount) = ;   % true frequency (9 bits)
%                                 stage2Cor(lru).headers(3,corPacketCount) = ;   % station id (9 bits)
%                                 stage2Cor(lru).headers(4,corPacketCount) = ;   % The LFAA packet count that this came from, divided by 2 (31 bits). (Divide by 2 since there are 4096 samples per filterbank output, but 2048 samples per LFAA input packet).
%                                 stage2Cor(lru).headers(5,corPacketCount) = ;   % First fine channel in the packet. Actual fine channel is this value * 96. (6 bits required)
%                                 
%                                 
%                                 stage1(lru).headers(1,packet) = virtualChannel;
%             stage1(lru).headers(2,packet) = frequency_id;
%             stage1(lru).headers(3,packet) = station_id;
%             stage1(lru).headers(4,packet) = packet_count;
                                
                                stage2Cor(lru).data(1:2:end,corPacketCount) = vPolCorrelatorFD((((split-1)*corSplitSize+1):((split-1)*corSplitSize + corSplitSize)),block);
                                stage2Cor(lru).data(2:2:end,corPacketCount) = hPolCorrelatorFD((((split-1)*corSplitSize+1):((split-1)*corSplitSize + corSplitSize)),block);
                            end
                        end
                        % PSS and PST operate on the same data, just need to trim the initialisation portion to the right length
                        %  -- TBD --
                        
                        
                        
                    end  % end of loop through all frames
                end
                %keyboard
            end  % end of loop through stations.
        end % end of loop through channels.

    end

    %% -----------------------------------------------------------------------
    % Stage3 Processing
    %  Fine corner turn, correlators and beamformers.
    
    
    %% -----------------------------------------------------------------------
    % Stage4 Processing
    %  * PSS corner turn
    %  * Packetising for SDP, PSS and PST    

end
%% ----------------------------------------------------------------------- 
% Supporting Functions

function [dout] = filterbank(din,filterTaps,FFTsize,oversampling,keep)
    % Polyphase filterbank.
    % Inputs:
    %   din : input data, assumes that the leading data is initialisation data for the filterbank.
    %   filterTaps : polyphase filter taps.
    %   fftSize : length of the FFT.
    %   oversampling : oversampling ratio for the output.
    %   keep : number of frequency channels to keep (e.g. 3456 out of a 4096 point FFT).
    % dout : output data. Length will be different to din as the outputs from the initialisation data are discarded.
    inputSamples = length(din);
    totalTaps = length(filterTaps);
    blocksPerOutput = totalTaps/FFTsize; % Number of blocks of input samples per output sample (expected to be 12)
    outputFFTs = oversampling * (inputSamples/FFTsize - (blocksPerOutput - 1)); % Number of times the FFT is done
    sampleStep = round(FFTsize / oversampling);
    dout = zeros(keep,outputFFTs);
    for f1 = 1:outputFFTs
        filtered1 = din(((f1-1)*sampleStep + 1):((f1-1)*sampleStep + blocksPerOutput*FFTsize)) .* filterTaps;
        filtered2 = reshape(filtered1,[FFTsize,blocksPerOutput]).';
        filtered3 = sum(filtered2);
        fd = fftshift(fft(filtered3));
        dout(:,f1) = fd((FFTsize/2 - keep/2):(FFTsize/2 + keep/2 - 1)).';  % With the -1 added to the second limit, this keeps one more negative frequency than positive frequency (as suggested by JohnB).
    end
    
end


function [z,x,y] = getLRUCoordinates(configuration,lru)
    % Get the coordinates of the LRU in the array (convert linear order to z,x,y)
    % lru should count from 1 up
    if (configuration == 0) % PISA, 3x1
        z = lru;
        x = 1;
        y = 1;
    elseif (configuration == 1) % AA1, 2x6x1
        z = mod(lru-1,2) + 1;
        x = floor((lru-1)/2) + 1;
        y = 1;
    elseif (configuration == 2) % AA2, 2x6x3
        z = mod(lru-1,2) + 1;
        x = mod(floor((lru-1)/2),6) + 1;
        y = floor((lru-1)/12) + 1;
    elseif (configuration == 3) % AA3-ITF, 4x6x2
        z = mod(lru-1,4) + 1;
        x = mod(floor((lru-1)/4),6) + 1;
        y = floor((lru-1)/24) + 1;
    elseif (configuration == 4) % AA3-CPF, 8x6x3
        z = mod(lru-1,8) + 1;
        x = mod(floor((lru-1)/8),6) + 1;
        y = floor((lru-1)/48) + 1;        
    else % AA4, 8x6x6
        z = mod(lru-1,8) + 1;
        x = mod(floor((lru-1)/8),6) + 1;
        y = floor((lru-1)/48) + 1;        
    end        
        
end