function [channelPower1Vpol,channelPower1Hpol, channelPower2Vpol, channelPower2Hpol, histogram] = get_LFAA_statistics(rundir,lru, maxPackets, histogramStation, histogramVirtualChannel)
    % Gets the LFAA signal statistics for a run.
    % Includes the power and voltage histograms
    %  rundir : directory for the data for the run
    %  lru    : which LRU to calculate statistics for
    %  maxPackets : Total number of input packets to process
    %  histogramStation : which of the two stations to log the voltage histogram for (1 or 2)
    %  histogramVirtualChannel : Which virtual channel to calculate the voltage histogram for.

    %% Load info for this run
    % Firmware register settings
    fid = fopen([rundir '/registerSettings.txt']);
    regjson = fread(fid,inf);
    fclose(fid);
    regjson = char(regjson');
    registers = jsondecode(regjson);
    
    % Load the data file
    fname = [rundir '/LFAA'];
    if exist(fullfile(cd,rundir,'LFAA.mat'),'file')
        load(fname);  % Should load a variable "fpga"
    else
        error(['Cannot find ' fname '. Run create_config.m to generate it']);
    end
    
    % Process the first "maxPackets" packets for lru "LRUSelect" to get power and voltage histograms
    d1 = fpga(lru);
    
    hsize = size(d1.headers); % hsize(1) should be 114 (number of bytes in an LFAA header, hsize(2) is the number of headers.
    if (maxPackets > hsize(1))
        error('requested processing more packets than there are in the data file')
    end
    
    dsize = size(d1.data);
    NChannels = dsize(2); % Number of non-zero channels.
    
    %% Setup the output data structures
    % station 1 channel power
    channelPower1Vpol = zeros(384,1);
    channelPower1Hpol = zeros(384,1);
    % station 2 channel power
    channelPower2Vpol = zeros(384,1);
    channelPower2Hpol = zeros(384,1);
    
    % Get the subset of samples that are used in the firmware for power measurements.
    PsamplesUsed = zeros(512,1);
    for k1 = 1:128
        PsamplesUsed((k1-1)*4 + 1) = (k1-1)*16 + 1;
        PsamplesUsed((k1-1)*4 + 2) = (k1-1)*16 + 6;
        PsamplesUsed((k1-1)*4 + 3) = (k1-1)*16 + 11;
        PsamplesUsed((k1-1)*4 + 4) = (k1-1)*16 + 16;
    end
    
    % Voltage histogram for a particular station and virtual channel
    histogram = zeros(1024,1);  % As per firmware, 1-256 = vpol real, 257-512 = vpol imaginary, 513-768 = hpol real, 769-1024 = hpol imaginary 
    
    %% 
    for packet = 1:maxPackets
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
            stationSel = 1;
        elseif (station_id == registers.LRU(lru).stationID(2))
            virtualChannel = find(registers.LRU(lru).virtualChannel(513:896) == virtualChannelContent);
            stationSel = 2;
        else
            error('No match for the stationID in the virtual channel table');
        end
        if (length(virtualChannel) ~= 1)
            error('No Match or more than 1 match in the virtual channel table');
        end
        virtualChannel = virtualChannel - 1;
        
        %% Get the voltage histogram
        if ((stationSel == histogramStation) && (virtualChannel == histogramVirtualChannel))
            % Only process the data part of the packet if the data is non-zero
            if (d1.dataPointers(packet,2) ~= 0)
                vpol = d1.data(d1.dataPointers(packet,1):2:(d1.dataPointers(packet,1) + 4095));
                hpol = d1.data((d1.dataPointers(packet,1)+1):2:(d1.dataPointers(packet,1) + 4095));
                for k1 = 1:2048
                    if (real(vpol(k1)) < 0)
                        VreIndex = 256 + double(real(vpol(k1)));
                    else 
                        VreIndex = double(real(vpol(k1)));
                    end
                    VreIndex = VreIndex + 1;

                    if (imag(vpol(k1)) < 0)
                        VimIndex = 256 + double(imag(vpol(k1)));
                    else 
                        VimIndex = double(imag(vpol(k1)));
                    end
                    VimIndex = VimIndex + 1;
                    
                    if (real(hpol(k1)) < 0)
                        HreIndex = 256 + double(real(hpol(k1)));
                    else 
                        HreIndex = double(real(hpol(k1)));
                    end
                    HreIndex = HreIndex + 1;

                    if (imag(hpol(k1)) < 0)
                        HimIndex = 256 + double(imag(hpol(k1)));
                    else 
                        HimIndex = double(imag(hpol(k1)));
                    end
                    HimIndex = HimIndex + 1;
                    
                    histogram(VreIndex) = histogram(VreIndex) + 1;
                    histogram(VimIndex + 256) = histogram(VimIndex + 256) + 1;
                    histogram(HreIndex + 512) = histogram(HreIndex + 512) + 1;
                    histogram(HimIndex + 768) = histogram(HimIndex + 768) + 1;
                    
                end
            end
        end
        
        %% Update the total power measurement
        % Only process the data part of the packet if the data is non-zero
        if (d1.dataPointers(packet,2) ~= 0)
            vpol = d1.data(d1.dataPointers(packet,1):2:(d1.dataPointers(packet,1) + 4095));
            hpol = d1.data((d1.dataPointers(packet,1)+1):2:(d1.dataPointers(packet,1) + 4095));
            
            % Select samples to match the firmware
            vpolSamples = abs(double(vpol(PsamplesUsed)));
            hpolSamples = abs(double(hpol(PsamplesUsed)));
            
            if (stationSel == 1)
                %keyboard
                channelPower1Vpol(virtualChannel+1) = channelPower1Vpol(virtualChannel+1) + sum(vpolSamples.^2);
                channelPower1Hpol(virtualChannel+1) = channelPower1Hpol(virtualChannel+1) + sum(hpolSamples.^2);
            else
                channelPower2Vpol(virtualChannel+1) = channelPower2Vpol(virtualChannel+1) + sum(vpolSamples.^2);
                channelPower2Hpol(virtualChannel+1) = channelPower2Hpol(virtualChannel+1) + sum(hpolSamples.^2);
            end
        end
        
    end
    
    figure(1);
    clf;
    hold on;
    grid on;
    hxvals = [-128:127];
    plot(hxvals,fftshift(histogram(1:256)),'r.-');
    plot(hxvals,fftshift(histogram(257:512)),'g.-');
    plot(hxvals,fftshift(histogram(513:768)),'b.-');
    plot(hxvals,fftshift(histogram(769:1024)),'c.-');
    legend('Vpol real','Vpol Imag','Hpol real','Hpol imag');
    title(['Histogram of raw samples, virtual channel = ' num2str(histogramVirtualChannel)])
    
    keyboard
    