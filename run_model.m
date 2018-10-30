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
    %        save results to stage2.mat
    %
    % 3. stage3 Processing:
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
            zStart = 1;  % index of the first LRU with the same z coordinate as this LRU
        elseif (modelConfig.configuration == 1) % AA1, 2x6x1, 4 stations/LRU
            virtualChannels = ((cz-1):2:383);
            stations = modelConfig.stationMap((floor((lru-1)/2) * 4 + 1):(floor((lru-1)/2) * 4 + 4));
            zStart = floor((lru-1)/2) * 2 + 1;
        elseif (modelConfig.configuration == 2) % AA2, 2x6x3, 4 stations/LRU
            virtualChannels = ((cz-1):2:383);
            stations = modelConfig.stationMap((floor((lru-1)/2) * 4 + 1):(floor((lru-1)/2) * 4 + 4));
            zStart = floor((lru-1)/2) * 2 + 1;
        elseif (modelConfig.configuration == 3) % AA3-ITF, 4x6x2, 8 stations/LRU
            virtualChannels = ((cz-1):4:383);
            stations = modelConfig.stationMap((floor((lru-1)/4) * 8 + 1):(floor((lru-1)/4) * 8 + 8));
            zStart = floor((lru-1)/4) * 4 + 1;
        elseif (modelConfig.configuration == 4) % AA3-CPF, 8x6x3, 16 stations/LRU
            virtualChannels = ((cz-1):8:383);
            stations = modelConfig.stationMap((floor((lru-1)/8) * 16 + 1):(floor((lru-1)/8) * 16 + 16));
            zStart = floor((lru-1)/8) * 8 + 1;
        else % AA4, 8x6x6, 16 stations/LRU
            virtualChannels = ((cz-1):8:383);
            stations = modelConfig.stationMap((floor((lru-1)/8) * 16 + 1):(floor((lru-1)/8) * 16 + 16));
            zStart = floor((lru-1)/8) * 8 + 1;
        end

        % Coarse corner turn inputs and outputs - just collect the appropriate packets and save to a data structure.
        if (modelConfig.keepCoarseCornerTurner)
            % TBD 
        end

        %keyboard
        % Correlator filterbank
        % Step through the channels and stations that this LRU processes
        for stationCount = 1:length(stations)
            station = stations(stationCount);
            for channelCount = 1:length(virtualChannels)
                channel = virtualChannels(channelCount);
                % Work out which LRU this combination of station and channel comes from.
                % Z coordinate is just stationCount/2, since "stations" is a list of all the stations in that go to this z coordinate (and 2 stations/LRU).
                srcLRU = floor((stationCount-1)/2) + (zStart-1) + 1;
                packetIndex = find(channel == stage1(srcLRU).headers(1,:),1);  % First packet for this station and channel
                if (stage1(srcLRU).dataPointers(packetIndex,2) ~= 0)
                    % Data exists for this channel
                    % Process it in 0.9 second frames, taking into account the coarse delay, and using zeros for initialisation data for the first frame.
                    % data is stage1(srcLRU).data(:,stage1(srcLRU).dataPointers(packetIndex,2))
                end
                keyboard
            end
        end

    end

    %% -----------------------------------------------------------------------
    % Stage3 Processing


end
%% ----------------------------------------------------------------------- 
% Supporting Functions

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