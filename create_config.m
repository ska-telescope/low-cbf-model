function create_config(rundir)
% Create setup files for the LOW.CBF packetised model.
% Configuration is stored in directory rundir
%
%% The setup that is generated consists of the following files:
% - "lmc_config.txt"
%    This is in JSON format.
%    It describes the setup for LOW.CBF according to the parameters specified
%    in the LMC (local monitoring and control) interface specification.
%
% - "register_settings.txt"
%    Register settings used in the model. For example the 4th order delay 
%    polynomials in lmc_config.txt are converted to a linear interpolation
%    for use in the firmware.
%
% - "sky.mat"
%    matlab format data file with simulated LFAA data.
%    All numbers are 8 bits, as per the LFAA spec.
%    The sky is included in the setup because it can take a significant 
%    amount of time to generate, and the sky must match the delays/beams 
%    configured for low.cbf
%
%% Source files
% This script uses the following files:
% - "stations.txt"
%    Locations of the stations on earth. (Each station is an array of antennas.)
%    JSON format with two variables:
%    - center : Latitude and longitude of the midpoint of the
%      array. Longitude is included for completeness but is ignored because
%      the date/time of the observation is arbitrarily set such that objects
%      in sky_definition.txt with a right ascension of 0 are on the horizon 
%      at t=0 (observation time is also set in sky_definition).
%    - offsets : Nx2 array. Each row is the location of a station,
%      relative to station_center. The first number is the east/west distance,
%      with east positive. The second number is the north/south distance, with 
%      north positive. Distances are in meters.
%
% - "sky.txt"
%    JSON format defining the objects in the sky that the telescope is looking
%    at. The sky model supports objects consisting of a collection of
%    sinusoids with pseudo-random phases located at points on the celestial
%    sphere with fixed equatorial coordinates (i.e. no accounting for 
%    precession, proper motion etc.).
%    Each object in the sky has the following properties:
%      - ascension : single double precision number specifying right ascension
%                    in decimal degrees e.g. 15.12325125
%      - declination : single double precision number as per right ascension
%      - coarse : LFAA coarse channel, numbered 1 to 384
%      - fine : low.cbf center fine channel for the source, numbered 1 to 3456
%      - SNR : source SNR in dB. If the object is an image, then the SNR is
%        set for the entire image, i.e. the ratio of the total power across
%        all sinusoids generated to the total noise power.
%      - sinusoids : Number of sinusoids to generate for each source
%      - BW : bandwidth in fine channels. Sinusoids will be spread out over
%        the bandwidth specified. This can be a fractional number of fine
%        channels (so beat frequencies can be arbitrarily low).
%      - sky_image : either an empty string, or the name of an image file. The
%        image defines a region of the sky, which will be centered at the
%        right ascension and declination defined above. The image should have 
%        the same number of pixels in each dimension. The image should be 
%        greyscale. Each pixel in the image is treated as a point source,
%        with relative brightness defined by the greyscale value.
%        If no image is specified, then a single point source is generated.
%      - sky_image_dimension : Number of degrees on an edge of the sky_image
%       (if used)
%      - index : a unique number used to identify this object in the sky.
%    Time
%      - This is the time to start viewing the objects, measured in hours.
%        The "true" date/time is arbitrarily defined such that an object 
%        with a right ascension of 0 will be on the horizon if time is 0. 
%        So, e.g. if time is set to 6 hours, then objects with 
%        right ascension 0 will be directly overhead.
%
% - "ArrayConfig.txt"
%   JSON format file which defines arrays and subarrays in use, and links
%   objects in sky_definition to an array, subarray, PSS beam or PST beam.
%   Two structures are specified in this file:
%   1. "array_definition" : Each array is configured with parameters:
%     - Array id  : 1 to 256 (?)
%     - subarray id : 1 to 256 (?)
%     - stations : List of stations to use for this array. Station numbers
%       correspond to the index into the station_offsets array in
%       "station_locations.txt"
%     - coarse_channels : List of coarse channels to use with this
%       array/subarray. Numbered 1 to 384. Coarse channels are always in
%       groups of 8, so legal values are 1, 9, 17 etc.
%   2. "sky_array_map" : Defines sky objects viewed by each array.
%     Parameters are
%     - Array_id
%     - subarray_id
%     - type : either correlator, PSS or PST.
%     - index : 1 to 512 for PSS or 1 to 16 for PST. (ignored for
%       correlator).
%     - sky_object : defines which object from "sky_definition.txt" to
%       view, using the "index" field from the sky object definition. 
%
% - "ModelConfig.txt"
%   JSON format file with any configuration that doesn't fit in any other category.

%% Load sky definition
fid = fopen([rundir '/sky.txt']);
sky_json = fread(fid,inf);
fclose(fid);
sky_json = char(sky_json');
sky1 = jsondecode(sky_json);
% expandSky adds the field "subIndex", and expands entries in sky that contain an image into a set of objects at each non-zero pixel in the image.
sky = expandSky(sky1.sky,rundir);


%% Load Model setup
modelConfig = parseModelConfig(rundir);


%% Check if lmcConfig.txt is out of date, if so then generate it and save it.
lmcOutOfDate = 0;
if exist(fullfile(cd,rundir,'lmcConfig.txt'),'file')
    d0 = dir(fullfile(cd,rundir,'lmcConfig.txt'));
    d1 = dir(fullfile(cd,rundir,'modelConfig.txt'));
    d2 = dir(fullfile(cd,rundir,'arrayConfig.txt'));
    d3 = dir(fullfile(cd,rundir,'sky.txt'));
    d4 = dir(fullfile(cd,rundir,'stations.txt'));
    if (d0.datenum < max([d1.datenum, d2.datenum, d3.datenum, d4.datenum]))
        lmcOutOfDate = 1;
    end
else
    lmcOutOfDate = 1;
end

%% Check is LFAA.mat is out of date.
skyOutOfDate = 0;
if exist(fullfile(cd,rundir,'LFAA.mat'),'file')
    d0 = dir(fullfile(cd,rundir,'LFAA.mat'));
    d1 = dir(fullfile(cd,rundir,'modelConfig.txt'));
    d2 = dir(fullfile(cd,rundir,'arrayConfig.txt'));
    d3 = dir(fullfile(cd,rundir,'sky.txt'));
    d4 = dir(fullfile(cd,rundir,'stations.txt'));
    if (d0.datenum < max([d1.datenum, d2.datenum, d3.datenum, d4.datenum]))
        skyOutOfDate = 1;
    end
else
    skyOutOfDate = 1;
end

%% Check if the register settings are out of date.
regOutOfDate = 0;
if exist(fullfile(cd,rundir,'registerSettings.txt'),'file')
    d0 = dir(fullfile(cd,rundir,'registerSettings.txt'));
    d1 = dir(fullfile(cd,rundir,'modelConfig.txt'));
    d2 = dir(fullfile(cd,rundir,'arrayConfig.txt'));
    d3 = dir(fullfile(cd,rundir,'sky.txt'));
    d4 = dir(fullfile(cd,rundir,'stations.txt'));
    if (d0.datenum < max([d1.datenum, d2.datenum, d3.datenum, d4.datenum]))
        regOutOfDate = 1;
    end
else
    regOutOfDate = 1;
end

%% Generate the full LMC Configuration if required, or load it if it is already up to date.

if (lmcOutOfDate)
    disp('Generating lmcConfig as it is out of date.')
    arrayConfigFull = parseArrayConfig(rundir,sky,modelConfig);
    %% -----------------------------------------------------------------------
    % Save lmcConfig as a json file in the run directory
    lmcdata = jsonencode(arrayConfigFull);
    fid = fopen([rundir '/lmcConfig.txt'],'w');
    fprintf(fid,lmcdata);
    fclose(fid);
else
    % lmcConfig was up to date, load from file.
    disp('LMC data is up to date, loading from lmcConfig.txt');
    fid = fopen([rundir '/lmcConfig.txt']);
    lmcjson = fread(fid,inf);
    fclose(fid);
    lmcjson = char(lmcjson');
    arrayConfigFull = jsondecode(lmcjson);
end

%% -----------------------------------------------------------------------
% 
if (skyOutOfDate)
    disp('Generating LFAA data as it is out of date.')
    % Step through each stationBeam and create data.
    % for each stationbeam
    %   for each station
    %     for each coarse channel
    %
    %
    % Data is stored as 
    % fpga(1:N) : Array with one entry for each FPGA. Each entry is a structure with fields as below, where M is the total number of packets to send.
    %   .headers : List of packet headers to send. Dimensions 114 x M. Data type is uint8.
    %   .txTimes : Time to send each packet in nanoseconds. Dimensions M x 1. Doubles.
    %   .dataPointers : offset into .data field to take the data part of the packet from. Dimensions Mx2.
    %                   At each row, first entry is the row offset, second entry is the column offset.
    %                   column offset value of 1 or more selects a column from data.
    %                   column offset value of 0 indicates filling the data part with zeros.
    %   .data : complex int8 array with data, dimensions (2*M*2048) x (2*[active channels]). 
    %           Dual polarisation with even and odd row indexes selecting the polarisation.
    %           Each column is for one station.
    %
    
    % First, get a list of all the stations
    stationList = [];
    for s1 = 1:length(arrayConfigFull.stations)
        stationList = [stationList arrayConfigFull.stations(s1).stationIDs];
    end
    stationList = unique(stationList);
    
    % Create the data as it comes from the sky
    % i.e. without taking into account the delays at each station
    %keyboard
    for s1 = 1:length(sky)
        % Add two fields to each entry in the sky structure:
        %  LFAAChannels : a list of coarse channels that are non-zero. N integers in the range 65-448.
        %  rawData : Mx(2N) data structure with complex single precision data, where M is the number of samples generated and N is the number of non-zero coarse channels.
        %            The factor of 2 in 2N for dual polarisation.
        %            Data is stored in single precision to reduce the storage required. Data is only generated if sky(s1).power > 0, 
        %            but the signal level of the generated data is NOT modulated by sky(s1).power. Power needs to be taken into account
        %            when the data received at a given station is calculated later. Note that the received data at a particular station is the sum of
        %            the all the different sky.data fields with a given sky.index (and different sky.subIndex)
        if (sky(s1).power > 0)
            % start and end frequency for the signal, in terms of fine channels.
            % The Correlator filterbank createse 4096 fine channels. Due to the 32/27 oversampling, 
            % there are 27/32 * 4096 = 3456 fine channel between coarse channels.
            frequencyStart = sky(s1).coarse * 3456 + sky(s1).fine - sky(s1).BW/2;  % Note that frequency start is in units of correlator fine channels, so the actual frequency is frequencyStart * (1/1080e-9)/4096
            frequencyEnd = sky(s1).coarse * 3456 + sky(s1).fine + sky(s1).BW/2;
            coarseStart = round(frequencyStart / 3456);
            fineStart = frequencyStart - coarseStart * 3456;
            coarseEnd = round(frequencyEnd/3456);
            fineEnd = frequencyEnd - coarseEnd * 3456;
            % runtime is in integration periods of 0.9sec = 408 LFAA packets = 408 * 2048 coarse channel samples.
            Nsamples = modelConfig.runtime * 408 * 2048;
            Nsamples = 32768 * (1 + ceil(Nsamples/32768)); % Need a few extra samples for resampling. Round up to a multiple of 32768. 
            if (sky(s1).sinusoids == 0)
                % put in the list of coarse channels
                sky(s1).LFAAChannels = (coarseStart:coarseEnd);
                sky(s1).rawData = zeros(Nsamples,2*length(sky(s1).LFAAChannels),'single');
                coarseCount = 0;
                % Signal is filtered noise
                for coarse = coarseStart:coarseEnd
                    
                    if (coarse == coarseStart)
                        currentFineStart = fineStart;
                    else
                        currentFineStart = -3456/2 - 10;  % When using the whole channel, pad it out a bit (i.e. the -10) so that after doppler correction the whole channel is still covered.
                    end
                    if (coarse == coarseEnd)
                        currentFineEnd = fineEnd;
                    else
                        currentFineEnd = 3456/2 + 10;  % When using the whole channel, pad it out a bit (i.e. the +10) so that after doppler correction the whole channel is still covered.
                    end
                    % Get some random data and filter it to the bandwidth required for this coarse channel.
                    d1 = (1/sqrt(2)) * (randn(Nsamples,2) + 1i * randn(Nsamples,2));  % Two polarisations. Each will have variance 1.
                    d1fft = fft(d1);
                    % Blank out frequencies before currentFineStart
                    if (round(Nsamples + 1 + Nsamples*currentFineStart/4096) <= Nsamples)
                        % currentFineStart was negative
                        d1fft((Nsamples/2 + 1):round(Nsamples + 1 + Nsamples*currentFineStart/4096),:) = 0;
                    else
                        % currentFineStart is positive, blank out all negative frequencies
                        d1fft((Nsamples/2+1):Nsamples,:) = 0;
                        % Blank out positive frequencies up to currentFineStart
                        d1fft(1:(Nsamples*currentFineStart/4096),:) = 0;
                    end
                    % Blank out frequencies after currentFineEnd
                    if (currentFineEnd < 0)
                        d1fft(round(Nsamples + 1 + Nsamples*currentFineEnd/4096):Nsamples,:) = 0;
                        d1fft(1:(Nsamples/2),:) = 0;
                    else
                        % currentFineEnd is positive
                        d1fft(round(Nsamples*currentFineEnd/4096):(Nsamples/2),:) = 0;
                    end
                    
                    d1 = ifft(d1fft);
                    % Force each polarisation in d1 to have variance 1.
                    d1(:,1) = d1(:,1)/std(d1(:,1));
                    d1(:,2) = d1(:,2)/std(d1(:,2));
                    
                    % Make the polarization match p = sky(s1).polarisation
                    % p is stokes parameters p = [I Q U V]
                    % Let x,y be the horizontal and vertical components of the polarisation
                    % So I = E(|x|^2) + E(|y|^2)
                    %    Q = E(|x|^2) - E(|y|^2)
                    %  U-iV = 2 * E(xy')
                    %
                    % Now given two complex signals a, b with E(|a|^2) = 1, E(|b|^2) = 1
                    % Then set 
                    %   x = na
                    %   y = nka + mb
                    % This will give the stokes parameters [I, Q, U, V] providing
                    %   n = sqrt((I+Q)/2)
                    %   k = (1/n^2) * (U-iV)/2
                    %   m = (I-Q)/2 - n^2*k^2
                    % 
                    n = sqrt((sky(s1).polarisation(1) + sky(s1).polarisation(2))/2);
                    if (n < 0.000001)
                        k = 0;
                    else
                        k = (1/n^2) * (sky(s1).polarisation(3) - 1i * sky(s1).polarisation(4))/2;
                    end
                    m = sqrt((sky(s1).polarisation(1) - sky(s1).polarisation(2))/2 - n^2*(abs(k)^2));
                    
                    d2 = zeros(size(d1));
                    d2(:,1) = n*d1(:,1);
                    d2(:,2) = k*d2(:,1) + m * d1(:,2);
                    
                    sky(s1).rawData(:,(coarseCount+1:coarseCount+2)) = single(d2);
                    coarseCount = coarseCount + 2; % Which column in sky(s1).rawData to put the signal into
                    % Check
                    %IQUV(1) = var(d2(:,1)) + var(d2(:,2));
                    %IQUV(2) = var(d2(:,1)) - var(d2(:,2));
                    %IQUV(3) = 2*mean(d2(:,1).*conj(d2(:,2)));
                    %disp(IQUV)
                    %keyboard
                    
                end
            else
                % Signal is a set of sinusoids
                % put in the list of coarse channels
                sky(s1).LFAAChannels = (coarseStart:coarseEnd);
                sky(s1).rawData = zeros(Nsamples,2*length(sky(s1).LFAAChannels),'single');
                coarseCount = 0;
                % Signal is a set of sinusoids
                % sky(s1).sinusoids = 1 => 1 sinusoid at the center frequency
                % sky(s1).sinusoids = 2 => frequencies are [center - BW/2,center + BW/2]
                % etc.
                % Get a list of frequencies for the sinusoids
                if (sky(s1).sinusoids == 1)
                    flist = sky(s1).coarse * 3456 + sky(s1).fine;
                else
                    flist = [frequencyStart:(frequencyEnd-frequencyStart)/(sky(s1).sinusoids - 1):frequencyEnd];
                end
                for coarse = coarseStart:coarseEnd
                    % Get the frequencies in units of fine channels for the start and end of this coarse channel,
                    % and find anything in flist in that range
                    fstart = round(coarse * 3456 - 1728);
                    fend = round(coarse * 3456 + 1728);
                    currentflist = flist(find((flist > fstart) && (flist < fend)));
                    for f1 = 1:length(currentflist)
                        thisF = currentflist(f1) - coarse * 3456; % the frequency in fine channels, offset from the center frequency of this coarse channel
                        d1 = exp(1i*(0:(Nsamples-1)) * (1080e-9) * 2 * pi * (thisF * (1/1080e-9) / 4096));
                        % Make the polarization match p = sky(s1).polarisation
                        % p is stokes parameters p = [I Q U V]
                        % Let x,y be the horizontal and vertical components of the polarisation
                        % So I = E(|x|^2) + E(|y|^2)
                        %    Q = E(|x|^2) - E(|y|^2)
                        %  U-iV = 2 * E(xy')
                        %
                        % We can't actually make a finite set of sinusoids uncorrelated, so we just make the correlation term (U-iV) correct.
                        % This means that we only generate a single source signal and small values for the correlation forces one of the
                        % polarisations to have lower power.
                        % So, given a complex signal a with E(|a|^2) = 1
                        % Then set 
                        %   x = na    (with n real)
                        %   y = ma    (with m complex)
                        % This will give the stokes parameters [I, Q, U, V]
                        %   n = sqrt((I+Q)/2)
                        %  |m| = sqrt((I-Q)/2)
                        %  angle(m) = angle(U-iV)
                        % (Providing the Stokes parameters imply the signal is fully polarised)
                        n = sqrt((sky(s1).polarisation(1) + sky(s1).polarisation(2))/2);
                        mabs = sqrt((sky(s1).polarisation(1) - sky(s1).polarisation(2))/2);
                        mphase = angle((sky(s1).polarisation(3) - sky(s1).polarisation(4))/2);
                        m = mabs * exp(1i*mphase);
                        d2 = zeros(Nsamples,2);
                        d2(:,1) = n*d1;
                        d2(:,2) = m*d1;
                        sky(s1).rawData(:,(coarseCount+1:coarseCount+2)) = sky(s1).rawData(:,(coarseCount+1:coarseCount+2)) + single(d2);
                        coarseCount = coarseCount + 2; % Which column in sky(s1).rawData to put the signal into
                    end
                    % Check
                    %IQUV(1) = var(sky(s1).rawData(:,1)) + var(sky(s1).rawData(:,2));
                    %IQUV(2) = var(sky(s1).rawData(:,1)) - var(sky(s1).rawData(:,2));
                    %IQUV(3) = 2*mean(sky(s1).rawData(:,1).*conj(sky(s1).rawData(:,2)));
                    %disp(IQUV)
                    %keyboard
                end
            end
        end
    end
    
    %keyboard
    % Each FPGA gets data for 2 stations, first step through FPGA
    for fpgaIndex = 1:floor(length(stationList)/2)
        % Each FPGA gets data from up to 8 logical stations (2 stations x 4 substations)
        % Get a list of logical stations for this FPGA.
        disp(['Generating data for FPGA ' num2str(fpgaIndex)])
        logicalIDs = [];
        stationIDs = [];   % Each logicalID has a unique combination of stationID and a substationID
        substationIDs = [];
        channelCount = 0;   % count of the number of channels being sent to this FPGA. Can get up to 2*384.
        usedChannelCount = 0;  % count of the number of channels being sent to this FPGA that are non-zero.
        resampledPoints = modelConfig.runtime * 408 * 2048; % Number of sample points to generate
        
        if (((fpgaIndex-1)*2 + 2) <= length(stationList)) % two stations sending data to this FPGA
            stationsToThisFPGA = [stationList((fpgaIndex-1)*2 + 1) stationList((fpgaIndex-1)*2 + 2)];
        else
            stationsToThisFPGA = stationList((fpgaIndex-1)*2 + 1);  % one station sending data to this FPGA
        end
        for station = stationsToThisFPGA
            for s1 = 1:length(arrayConfigFull.stations)
                for s2 = 1:length(arrayConfigFull.stations(s1).stationIDs)
                    if (station == arrayConfigFull.stations(s1).stationIDs(s2))
                        logicalIDs = [logicalIDs arrayConfigFull.stations(s1).logicalIDs(s2)];
                        stationIDs = [stationIDs station];
                        substationIDs = [substationIDs arrayConfigFull.stations(s1).substationIDs(s2)];
                    end
                end
            end
        end
        
        % Find the number of channels coming to this FPGA
        totalChannels = 0;
        for logicalIndex = 1:length(logicalIDs)
            logicalID = logicalIDs(logicalIndex);
            subArrayIndex = -1;
            for sIndex = 1:length(arrayConfigFull.subArray)
                if any(logicalID == arrayConfigFull.subArray(sIndex).logicalIDs)
                    subArrayIndex = arrayConfigFull.subArray(sIndex).index;
                end
            end
            if (subArrayIndex ~= -1)
                % Find the stationBeams for this subArray, and add up the number of channels in each stationBeam
                for beamIndex = 1:length(arrayConfigFull.stationBeam)
                    if (arrayConfigFull.stationBeam(beamIndex).subArray == subArrayIndex)
                        totalChannels = totalChannels + length(arrayConfigFull.stationBeam(beamIndex).channels);
                    end
                end
            end
        end
        
        fpga(fpgaIndex).headers = zeros(114,totalChannels * modelConfig.runtime * 408,'uint8'); % Note runtime is in units of 0.9s integration periods; Each period is 408 packets worth. (408 * 2048 * 1080ns = 0.9024s)
        fpga(fpgaIndex).data = zeros(resampledPoints*2,1,'int8');  % x2 for dual polarisation.
        
        logicalChannel = 0;  % logical channel is used in the SPEAD header. It counts channels for a particular station.
        for logicalIndex = 1:length(logicalIDs)
            logicalID = logicalIDs(logicalIndex);
            stationID = stationIDs(logicalIndex);
            substationID = substationIDs(logicalIndex);
            if ((logicalIndex > 1) && (stationIDs(logicalIndex-1) ~= stationIDs(logicalIndex)))
                % New station, reset the logicalChannel count.
                logicalChannel = 0;
            end
            % Work out how many active channels there are for this logical station.
            % First find the subarray this logical station is in
            subArrayIndex = -1;
            for sIndex = 1:length(arrayConfigFull.subArray)
                if any(logicalID == arrayConfigFull.subArray(sIndex).logicalIDs)
                    subArrayIndex = arrayConfigFull.subArray(sIndex).index;
                end
            end
            if (subArrayIndex == -1)
                warning(['Logical station ' num2str(logicalID) ' is not assigned to a sub array']);
                channelList = [];
            else
                
                % Find the stationBeams for this subArray, and step through each of the channels for each stationBeam
                for beamIndex = 1:length(arrayConfigFull.stationBeam)
                    if (arrayConfigFull.stationBeam(beamIndex).subArray == subArrayIndex)
                        skyIndex = arrayConfigFull.stationBeam(beamIndex).skyIndex;
                        for channelI = 1:length(arrayConfigFull.stationBeam(beamIndex).channels)
                            channel = arrayConfigFull.stationBeam(beamIndex).channels(channelI);
                            % Count of the total number of channels coming to this FPGA
                            channelCount = channelCount + 1;
                            % Count of the channels coming from this station.
                            logicalChannel = logicalChannel + 1;
                            % Generate the first header for this logical station and channel in fpga(fpgaIndex).headers(1:114,channelCount)
                            %  channel is arrayConfigFull.stationBeam(beamIndex).channels(channelI) 
                            %  subArrayIndex is arrayConfigFull.subArray(XX).index
                            %  stationBeam is arrayConfigFull.stationBeam(beamIndex).index
                            % 
                            % Fields in each packet :
                            %  1:6    : Destination MAC 
                            %  7:12   : Source MAC
                            %  13:14  : EtherType = 0x0800 (indicates IPv4)
                            %  15:34  : IPv4 header
                            %  
                            %  35:42  : UDP header
                            %  43:114 : SPEAD header
                            hdr = zeros(114,1,'uint8');
                            hdr(1:6) = modelConfig.LFAAMAC(fpgaIndex,:);  % destination MAC address
                            hdr(7:12) = [0 0 0 0 0 0];    % source MAC address (low.CBF ignores this).
                            hdr(13:14) = [8 0];           % Ethertype
                            hdr(15) = hex2dec('45');      % IPv4, version + IHL
                            hdr(16) = 0;                  % IPv4, DSCP + ECN
                            hdr(17) = hex2dec('20');      % IPv4 length; Total length is (IPv4 Header) + (UDP header) + (SPEAD) + (data) = 20 + 8 + 72 + 8192 = 8292 bytes = 0x2064
                            hdr(18) = hex2dec('64');
                            hdr(19:20) = 0;               % IPv4, Identification
                            hdr(21:22) = 0;               % IPv4, Flags and fragment offset
                            hdr(23) = 8;                  % IPv4, time to live
                            hdr(24) = 17;                 % IPv4, protocol, 17 = UDP.
                            hdr(25:26) = 0;               % hdr(25:26) is the IPv4 checksum; fix this after the IP addresses are done.
                            hdr(27:30) = [192 168 255 1]; % Source IP address
                            hdr(31:34) = modelConfig.IP(fpgaIndex,:); % Destination IP address
                            % Calculate the IPv4 checksum, from hdr(15:34)
                            hdrd = double(hdr(15:34));
                            hdr16bit = 256 * hdrd(1:2:end) + hdrd(2:2:end);
                            hsum = sum(hdr16bit);
                            while (hsum > 65535)
                                hsum = mod(hsum,65536) + floor(hsum/65536);
                            end
                            hsum = 65535 - hsum;
                            hdr(25) = floor(hsum/256);
                            hdr(26) = mod(hsum,256);
                            % UDP header
                            hdr(35) = floor(modelConfig.LFAAUDPSrc/256);  % source UDP port
                            hdr(36) = mod(modelConfig.LFAAUDPSrc,256);
                            hdr(37) = floor(modelConfig.LFAAUDPDest/256); % destination UDP port
                            hdr(38) = mod(modelConfig.LFAAUDPDest,256);
                            hdr(39) = floor(8272/256);  % UDP Length = 8 (UDP header) + 72 (SPEAD header) + 8192 (data) = 8272
                            hdr(40) = mod(8272,256);
                            % calculate the UDP checksum
                            % temporarily put in 0 for the UDP checksum.
                            hdr(41) = 0;
                            hdr(42) = 0;
                            
                            % SPEAD - first 8 bytes
                            hdr(43) = 83;   % magic number, 0x53 = 83 decimal
                            hdr(44) = 4;    % version = 4
                            hdr(45) = 2;    % Item Pointer Width = 2
                            hdr(46) = 6;    % Heap Address Width = 6
                            hdr(47:48) = 0; % Reserved.
                            hdr(49) = 0;    % Number of Items
                            hdr(50) = 8;    
                            % hdr(51:58) = SPEAD heap_counter (ID = 0x0001)
                            hdr(51) = 128;  % msb = 1, indicates that heap_counter is immediate
                            hdr(52) = 1;    % heap_counter item ID
                            hdr(53) = floor((logicalChannel-1)/256);     % logical channel ID. Note it counts from 0.
                            hdr(54) = mod((logicalChannel-1),256);
                            hdr(55:58) = 0; % Packet counter. This counter is continuous for any given channel.
                            % hdr(59:66) = SPEAD pkt_len (ID = 0x0004)
                            hdr(59) = 128;  % msb = 1, indicates that pkt_len is immediate
                            hdr(60) = 4;    % item ID
                            hdr(61:66) = 0; % pkt_len is not used.
                            % hdr(67:74) = SPEAD sync_time (ID = 0x1027)
                            hdr(67) = 144;  % msb = 1, indicates sync_time is immediate, 0x80 + 0x10 = 0x90 = 144
                            hdr(68) = 39;   % 39 = 0x27
                            unixTime = posixtime(datetime(clock));  % Note this doesn't do timezones correctly.
                            hdr(69) = floor(unixTime/2^40);
                            hdr(70) = mod(floor(unixTime/2^32),256);
                            hdr(71) = mod(floor(unixTime/2^24),256);
                            hdr(72) = mod(floor(unixTime/2^16),256);
                            hdr(73) = mod(floor(unixTime/2^8),256);
                            hdr(74) = mod(unixTime,256);
                            % hdr(75:82) = SPEAD timestamp (ID = 0x1600). Note LFAA ICD is unclear - this either counts in ns or in units of 1.25ns
                            % Note 2048 samples per packet, with 1080 ns sampling period, so
                            %  * For 1 ns counter, this will increment by 2048 * 1080 = 2211840 between packets.
                            %  * For 1.25 ns counter, this will increment by 2048 * 1080/1.25 = 1769472 between packets.
                            %
                            hdr(75) = 150;  % msb = 1 for immediate value, and ID(15:8) = 0x16 so 0x96 = 150
                            hdr(76) = 0;
                            hdr(77) = floor(modelConfig.timestamp/2^40);
                            hdr(78) = mod(floor(modelConfig.timestamp/2^32),256);
                            hdr(79) = mod(floor(modelConfig.timestamp/2^24),256);
                            hdr(80) = mod(floor(modelConfig.timestamp/2^16),256);
                            hdr(81) = mod(floor(modelConfig.timestamp/2^8),256);
                            hdr(82) = mod(modelConfig.timestamp,256);
                            % hdr(83:90) = SPEAD center_freq (ID = 0x1011). Channel center frequency in Hz.  
                            hdr(83) = 144;  % 0x80 + 0x10 = 0x90 = 144;
                            hdr(84) = 17;   % 0x11
                            channelFrequency = channel * 781250; % Channel center frequency in Hz; Note that channel 65 (the first possible channel) has a center frequency of 50781250 Hz.
                            hdr(85) = floor(channelFrequency/2^40);
                            hdr(86) = mod(floor(channelFrequency/2^32),256);
                            hdr(87) = mod(floor(channelFrequency/2^24),256);
                            hdr(88) = mod(floor(channelFrequency/2^16),256);
                            hdr(89) = mod(floor(channelFrequency/2^8),256);
                            hdr(90) = mod(channelFrequency,256);
                            % hdr(91:98) = SPEAD csp_channel_info (ID = 0x3000)
                            hdr(91) = 176;   % msb = 1, then 0x30 from ID => 0x80 + 0x30 = 0xB0 = 176
                            hdr(92) = 0;     % low byte of ID = 0x00 
                            hdr(93) = 0;     % reserved = 0
                            hdr(94) = 0;     % reserved = 0
                            hdr(95) = 0;     % high byte of beam, always 0.
                            hdr(96) = arrayConfigFull.stationBeam(beamIndex).index;
                            hdr(97) = floor(channel/256);
                            hdr(98) = mod(channel,256);
                            % hdr(99:106) = SPEAD csp_antenna_info (ID = 0x3001)
                            hdr(99) = 176;   % msb = 1, then 0x30 from ID => 0x80 + 0x30 = 0xB0 = 176
                            hdr(100) = 1;    % low byte of ID = 0x01
                            hdr(101) = substationID;          % substation ID
                            hdr(102) = subArrayIndex;         % subarray ID
                            hdr(103) = floor(stationID/256);  % station ID (high byte)
                            hdr(104) = mod(stationID,256);    % station ID (low byte)
                            totalSubstations = 0;
                            for s1 = 1:length(arrayConfigFull.stations)
                                totalSubstations = totalSubstations + sum(arrayConfigFull.stations(s1).stationIDs == stationID);
                            end
                            nof_contributing_antenna = floor(256 / totalSubstations);  % 256 antennas divided equally among the substations
                            hdr(105) = floor(nof_contributing_antenna/256);            % nof_contributing_antenna (high byte)
                            hdr(106) = mod(nof_contributing_antenna,256);              % nof_contributing_antenna (low byte)
                            % hdr(107:114) = SPEAD sample_offset (ID = 0x3300) (not used, set to 0).
                            hdr(107) = 179;    % msb = 1, ID = 0x33 so 0xB3 = 179
                            hdr(108) = 0;
                            hdr(109:114) = 0;
                            fpga(fpgaIndex).headers(:,channelCount) = hdr;
                            
                            % Generate data for this logical channel and channel i.e. fpga(fpgaIndex).data(:,usedChannelCount)
                            % Need the sinusoid for the delay rather than the polynomial approximation.
                            skyIsNull = 1;
                            
                            for s1 = 1:length(sky)
                                hpol = zeros(resampledPoints,1);
                                vpol = zeros(resampledPoints,1);
                                % Find each object in the sky which is in the view of this beam, i.e. sky(s1).index == skyIndex
                                s2 = find(arrayConfigFull.stationsFull.logicalIDs == logicalID);
                                %thisStation.latitude = arrayConfigFull.stationsFull.latitude;
                                %thisStation.longitude = arrayConfigFull.stationsFull.longitude;
                                %thisStation.altitude = arrayConfigFull.stationsFull.altitude;
                                %thisStation.offsets = arrayConfigFull.stationsFull.offsets(s2,:);
                                
                                if ((sky(s1).index == skyIndex) && (~isempty(sky(s1).rawData)))   % When the sky comes from an image, the first entry in sky is for station beam pointing only and has no data associated with it.
                                    %delayFunction = getDelayFunctions(thisStation,sky(s1),modelConfig.time);
                                    % In fullDoppler mode, we don't ignore common motion, and we do put in the stationBeam doppler correction.
                                    % If not in fullDoppler mode, we ignore the common motion, so the only doppler inserted into the data is the rotation of the array relative to the sky.
                                    allDelayFunctions = getDelayFunctions(arrayConfigFull.stationsFull,sky(s1),modelConfig.time,~modelConfig.fullDoppler);
                                    
                                    % reformat delayFunction into the form required by resampleNU.
                                    % Note OAFPS = Offset, Amplitude, Frequency, Phase, Slope
                                    delayOAFPS(1) = allDelayFunctions.offset(s2) * 1e9;    % resampleNU requires nanseconds; getDelayFunctions returns seconds.
                                    delayOAFPS(2) = allDelayFunctions.amplitude(s2) * 1e9;
                                    delayOAFPS(3) = allDelayFunctions.rate(s2);
                                    delayOAFPS(4) = allDelayFunctions.phase(s2);
                                    if (modelConfig.fullDoppler)
                                        % Find the global Doppler correction. This corrects for common motion of all the stations in the array.
                                        % Slope is in ns/s, doppler from stationBeam is in m/s
                                        delayOAFPS(5) = 1e9 * arrayConfigFull.stationBeam(beamIndex).doppler / 2.99792458e8;
                                    else
                                        delayOAFPS(5) = 0;
                                    end
                                    % The vertical polarisation can have a different delay to account for different cable lengths in the receivers.
                                    delayOAFPSV = delayOAFPS;
                                    delayOAFPSV(1) = delayOAFPSV(1) + 1e9 * arrayConfigFull.stationsFull.polarisationOffsets(s2);
                                    % Get the data from the sky structure and resample it for this particular station.
                                    % channel - the LFAA channel this packet is supplying
                                    % sky(s1).LFAAChannels - list of coarse channels for this sky entry
                                    % sky(s1).rawData - raw data 
                                    if any(channel == sky(s1).LFAAChannels)
                                        thisData = find(channel == sky(s1).LFAAChannels);
                                        % Resample the data for this station, both polarisations
                                        disp(['resample data, station logical ID = ' num2str(logicalID) ' channel = ' num2str(channel) ' sky index = ' num2str(sky(s1).index) ' sky subIndex = ' num2str(sky(s1).subIndex)]);
                                        tic
                                        % oddly, double precision is faster than single precision. Unclear why.
                                        hpol = hpol + sqrt(sky(s1).power) * resampleNUfast(double(sky(s1).rawData(:,2*(thisData-1) + 1)),1080,channelFrequency,delayOAFPS,resampledPoints);
                                        vpol = vpol + sqrt(sky(s1).power) * resampleNUfast(double(sky(s1).rawData(:,2*(thisData-1) + 2)),1080,channelFrequency,delayOAFPSV,resampledPoints);
                                        skyIsNull = 0;
                                        toc
                                    end 
                                end
                            end
                            % Put the data in hpol and vpol into fpga(fpgaIndex).data(:,usedChannelCount)
                            % Note that the vpol and hpol samples alternate just as in the LFAA packet format.
                            % Note for fpga().dataPointers : 
                            %  .dataPointers = offset into .data field to take the data part of the packet from. Dimensions Mx2.
                            %                  At each row, first entry is the row offset, second entry is the column offset.
                            %                  column offset value of 1 or more selects a column from data.
                            %                  column offset value of 0 indicates filling the data part with zeros.
                            if (skyIsNull)
                                % just indicate that the data associated with this header is zeros.
                                fpga(fpgaIndex).dataPointers(channelCount,1) = 0;
                                fpga(fpgaIndex).dataPointers(channelCount,2) = 0;  % 0 for zeros.
                            else
                                % Scale hpol, vpol so that RMS is as specified
                                hpol = hpol * sky(s1).RMS / std(hpol);
                                vpol = vpol * sky(s1).RMS / std(vpol);
                                hpolReal = int8(real(hpol));
                                hpolImag = int8(imag(hpol));
                                vpolReal = int8(real(vpol));
                                vpolImag = int8(imag(vpol));
                                % Convert negative saturation to -127 instead of -128, as -128 is reserved for RFI.
                                hpolRealLow = find(hpolReal < -127);
                                hpolImagLow = find(hpolImag < -127);
                                vpolRealLow = find(vpolReal < -127);
                                vpolImagLow = find(vpolImag < -127);
                                
                                hpolReal(hpolRealLow) = -127; %#ok<FNDSB>
                                hpolImag(hpolImagLow) = -127; %#ok<FNDSB>
                                vpolReal(vpolRealLow) = -127; %#ok<FNDSB>
                                vpolImag(vpolImagLow) = -127; %#ok<FNDSB>
                                % put in the calculated data.
                                usedChannelCount = usedChannelCount + 1;
                                
                                fpga(fpgaIndex).data(1:2:end,usedChannelCount) = complex(vpolReal,vpolImag);  % resampledPoints*2,
                                fpga(fpgaIndex).data(2:2:end,usedChannelCount) = complex(hpolReal,hpolImag);
                                fpga(fpgaIndex).dataPointers(channelCount,1) = 1;   % This is for the first header, so start at the start of the data.
                                fpga(fpgaIndex).dataPointers(channelCount,2) = usedChannelCount;
                            end
                            
                        end
                    end
                end
            end
        end
        
        % At this point fpga(fpgaIndex).headers contains a list of headers for one instant in time.
        % Replicate and modify this list so that fpga(fpgaIndex).headers has all the headers for the 
        % Note : totalChannels is the number of headers for a particular time for a particular FPGA.
        %    fpga(fpgaIndex).headers = zeros(114,totalChannels * modelConfig.runtime * 408,'uint8');
        dp1 = fpga(fpgaIndex).dataPointers;
        fpga(fpgaIndex).dataPointers = zeros(totalChannels * modelConfig.runtime * 408,2);
        fpga(fpgaIndex).dataPointers(1:totalChannels,:) = dp1;
        
        for hdrCount = (totalChannels+1):(totalChannels * modelConfig.runtime * 408)
            % Copy the headers from the previous instance for this channel.
            fpga(fpgaIndex).headers(:,hdrCount) = fpga(fpgaIndex).headers(:,hdrCount - totalChannels);
            % Change things that are different from the previous header.
            % hdr(51:58) = SPEAD heap_counter (ID = 0x0001)
            % hdr(55:58) = 0; % Packet counter. This counter is continuous for any given channel.
            % add one to the existing packet counter
            curPacketCounter = double(fpga(fpgaIndex).headers(58,hdrCount)) + 256 * double(fpga(fpgaIndex).headers(57,hdrCount)) + 256*256*double(fpga(fpgaIndex).headers(56,hdrCount)) + 256*256*256*double(fpga(fpgaIndex).headers(55,hdrCount));
            curPacketCounter = curPacketCounter + 1;
            fpga(fpgaIndex).headers(55,hdrCount) = floor(curPacketCounter/(256*256*256));
            fpga(fpgaIndex).headers(56,hdrCount) = mod(floor(curPacketCounter/(256*256)),256);
            fpga(fpgaIndex).headers(57,hdrCount) = mod(floor(curPacketCounter/(256)),256);
            fpga(fpgaIndex).headers(58,hdrCount) = mod(curPacketCounter,256);

            % hdr(75:82) = SPEAD timestamp (ID = 0x1600). Note LFAA ICD is unclear - this either counts in ns or in units of 1.25ns
            %  Note 2048 samples per packet, with 1080 ns sampling period, so
            %   * For 1 ns counter, this will increment by 2048 * 1080 = 2211840 between packets.
            %   * For 1.25 ns counter, this will increment by 2048 * 1080/1.25 = 1769472 between packets.
            curTimeStamp = double(fpga(fpgaIndex).headers(82,hdrCount)) + (2^8)*double(fpga(fpgaIndex).headers(81,hdrCount)) + (2^16)*double(fpga(fpgaIndex).headers(80,hdrCount)) + ...
                           (2^24)*double(fpga(fpgaIndex).headers(79,hdrCount)) + (2^32)*double(fpga(fpgaIndex).headers(78,hdrCount)) + (2^40)*double(fpga(fpgaIndex).headers(77,hdrCount));
            curTimeStamp = curTimeStamp + 2048 * 1080/1.25;
            fpga(fpgaIndex).headers(77,hdrCount) = floor(curTimeStamp/(2^40));
            fpga(fpgaIndex).headers(78,hdrCount) = mod(floor(curTimeStamp/(2^32)),256);
            fpga(fpgaIndex).headers(79,hdrCount) = mod(floor(curTimeStamp/(2^24)),256);
            fpga(fpgaIndex).headers(80,hdrCount) = mod(floor(curTimeStamp/(2^16)),256);
            fpga(fpgaIndex).headers(81,hdrCount) = mod(floor(curTimeStamp/(2^8)),256);
            fpga(fpgaIndex).headers(82,hdrCount) = mod(curTimeStamp,256);
            % Put in the link to the data
            fpga(fpgaIndex).dataPointers(hdrCount,:) = fpga(fpgaIndex).dataPointers(hdrCount - totalChannels,:);
            fpga(fpgaIndex).dataPointers(hdrCount,1) = fpga(fpgaIndex).dataPointers(hdrCount,1) + 4096;  % 2048 samples per packet, data is interleaved between two polarisations, so plus 4096.
        end
        
        % Generate the transmit time for each header
        % Packets should be evenly spaced.
        % The number of packets that need to be sent in the time interval covered by the samples in the packet is "totalChannels"
        % The time interval covered by the samples is 2048 * 1080 = 2211840 ns
        packet_gap = 2048*1080 / totalChannels;  % Time gap between the start of successive packets in nanoseconds.
        fpga(fpgaIndex).txTimes = [0:(totalChannels * modelConfig.runtime * 408 - 1)] * packet_gap;
        
        % UDP checksum calculation for each header
        disp(['Calculating the UDP checksums for fpga ' num2str(fpgaIndex)]);
        tic
        for hdrCount = 1:(totalChannels * modelConfig.runtime * 408)
            pseudoHdr(1:8) = double(fpga(fpgaIndex).headers(27:34,hdrCount)); % src and destination IP address
            pseudoHdr(9) = 0;            % 0
            pseudoHdr(10) = 17;          % protocol
            pseudoHdr(11) = double(fpga(fpgaIndex).headers(39,hdrCount));     % UDP length field, same as length field in the UDP header.
            pseudoHdr(12) = double(fpga(fpgaIndex).headers(40,hdrCount));  
            pseudoHdr(13:18) = double(fpga(fpgaIndex).headers(35:40,hdrCount)); % source UDP port, destination UDP port, UDP length.
            pseudoHdr16bit = 256 * pseudoHdr(1:2:18) + pseudoHdr(2:2:18);
            % UDP checksum is header bytes 41, 42. The rest of the header (SPEAD) is from 43 to 114.
            restOfHdr16bit = 256 * double(fpga(fpgaIndex).headers(43:2:114,hdrCount)) + double(fpga(fpgaIndex).headers(44:2:114,hdrCount));
            % Data part
            if (fpga(fpgaIndex).dataPointers(hdrCount,2) == 0)
                % No data in .data, data is zeros.
                sumData16bit = 0;
            else
                data16bit = 256 * real(double(fpga(fpgaIndex).data(fpga(fpgaIndex).dataPointers(hdrCount,1),fpga(fpgaIndex).dataPointers(hdrCount,2)))) + ...
                                  imag(double(fpga(fpgaIndex).data(fpga(fpgaIndex).dataPointers(hdrCount,1),fpga(fpgaIndex).dataPointers(hdrCount,2))));
                sumData16bit = sum(data16bit);
            end
            hsum = sum(pseudoHdr16bit) + sum(restOfHdr16bit) + sumData16bit;
            while (hsum > 65535)
                hsum = mod(hsum,65536) + floor(hsum/65536);
            end
            hsum = 65535 - hsum;
            fpga(fpgaIndex).headers(41,hdrCount) = floor(hsum/256);  % Put the UDP checksum in the packet.
            fpga(fpgaIndex).headers(42,hdrCount) = mod(hsum,256);
            %keyboard
        end
        toc
        %disp('Done calculating UDP checksums');
        %keyboard
    end
    %keyboard
    savename = [rundir '/LFAA'];
    save(savename,'fpga');
else
    disp('LFAA data is already up to date in LFAA.mat');
end

%% 
%keyboard

%%
% Generate the register settings in the firmware from the LMC configuration
if (regOutOfDate)
    disp('Generating register settings as they are out of date.')
    [LRUreg, globalreg] = getRegisterSettings(arrayConfigFull,modelConfig);
    reg.LRU = LRUreg;
    reg.global = globalreg;
    regdata = jsonencode(reg);
    fid = fopen([rundir '/registerSettings.txt'],'w');
    fprintf(fid,regdata);
    fclose(fid);
else
    disp('Register settings in registerSettings.txt are already up to date.');
end

