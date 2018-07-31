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
%% 

%% Load station locations

fid = fopen([rundir '/stations.txt']);
station_locations_json = fread(fid,inf);
fclose(fid);
station_locations_json = char(station_locations_json');
stations = jsondecode(station_locations_json);


%% Load sky definition
fid = fopen([rundir '/sky.txt']);
sky_json = fread(fid,inf);
fclose(fid);
sky_json = char(sky_json');
sky1 = jsondecode(sky_json);
% expandSky adds the field "subIndex", and expands entries in sky that contain an image into a set of objects at each non-zero pixel in the image.
sky = expandSky(sky1.sky,rundir);


%% Load array setup
fid = fopen([rundir '/arrayConfig.txt']);
array_json = fread(fid,inf);
fclose(fid);
array_json = char(array_json');
arrayConfig = jsondecode(array_json);

% All about arrays...
% From LFAA we get (in each packet)
%  heap_counter 
%     -> logical_channel_id(47:32)
%     -> packet_counter(31:0)
%  csp_channel_info 
%     -> beam_id(31:16)
%     -> frequency_id(15:0)
%  csp_antenna_info 
%     -> substation_id(47:40)
%     -> subarray_id(39:32)
%     -> station_id(31:16)
%     -> nof_contributing_antenna(15:0)

%% Load Model setup
fid = fopen([rundir '/modelConfig.txt']);
model_json = fread(fid,inf);
fclose(fid);
model_json = char(model_json');
modelConfig = jsondecode(model_json);

% Put in LFAAMAC field if it is not already there.
if ~isfield(modelConfig,'LFAAMAC')
    modelConfig.LFAAMAC = zeros(300,6,'uint8');
    modelConfig.LFAAMAC(:,1) = 1;
    modelConfig.LFAAMAC(:,2) = 2;
    modelConfig.LFAAMAC(:,3) = 3;
    modelConfig.LFAAMAC(:,6) = 0;
    for n = 1:300
        modelConfig.LFAAMAC(n,4) = floor((n-1)/256);
        modelConfig.LFAAMAC(n,5) = mod((n-1),256);
    end
end

% Put in the IP field if it is not already there
if ~isfield(modelConfig,'IP')
    modelConfig.IP = zeros(300,4,'uint8');
    modelConfig.IP(:,1) = 192;
    modelConfig.IP(:,2) = 168;
    for n = 1:300
        modelConfig.IP(n,3) = floor((n-1)/256);
        modelConfig.IP(n,4) = mod((n-1),256);
    end
end

% Put in the LFAAUDPSrc and LFAAUDPDest fields (source and destination UDP ports for SPEAD packets from LFAA
if ~isfield(modelConfig,'LFAAUDPSrc')
    modelConfig.LFAAUDPSrc = 2000;
end
if ~isfield(modelConfig,'LFAAUDPDest')
    modelConfig.LFAAUDPDest = 2001;
end

% Put in timestamp. This is the offset of the first sample in the first packet from sync_time in the LFAA SPEAD header.
if ~isfield(modelConfig,'timestamp')
    modelConfig.timestamp = 0;
end


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

%%
if (lmcOutOfDate)
     disp('Generating lmcConfig as it is out of date.')
    %% 
    % Generate the LMC configuration
    %
    % Go through arrayConfig, do error checking, and insert any missing optional fields using the default values.
    % Also replace the skyIndex and skySubindex fields with delay polynomials.
    %
    % Top level fields in arrayConfig.txt are:
    %  "global" : These parameters are kept as is.
    %  "stations" : Combined with the information in stations.txt to generate a revised version ("stationsFull"),
    %               which describes the location, station id and substation id for each logical station.
    %  "subArray" : Specifies sub arrays, copied as is with error checking.
    %  "stationBeam" : combined with sky information, used to generate a revised version with delay polynomials for each stationBeam.
    %  "correlator" : copied as is, with error checking.
    %  "PSSBeam" : combined with sky information, used to generate a revised version with delay polynomials.
    %  "PSTBeam" : Similar to PSSBeam.
    %  "VLBIBeam" : Similar to PSSBeam.

    %% -- "global" -------------------------------------------------------------
    % copy global directly
    g1 = arrayConfig.global(1);
    if iscell(g1)
        g1Struct = g1{1};
    else
        g1Struct = g1;
    end
    arrayConfigFull.global = g1Struct;

    %% -- "stations" -------------------------------------------------------------
    % sets logical IDs, locations and station types for all the stations.
    % Two structures are generated :
    %  stationsFull : Used for generating delay polynomials later in this function. Also stored in arrayConfigFull for future reference.
    %  arrayConfigFull.stations : Part of the LMC configuration (written back to the file "lmcConfig.txt").
    %
    % "StationsFull" is a revised version of the stations structure, with fields:
    %  .latitude, .longitude, .altitude : Copied directly from stations structure
    %  .logicalIDs    : Vector of all logical IDs that are defined
    %  .stationIDs    : Vector of the same length as .logicalIDs. Holds the station ID for this station (1 to 512).
    %  .substationIDs : Vector of the same length as .logicalIDs. Holds the substation ID for this station (1 to 4).
    %  .stationType   : Vector of the same length as logicalIDs, 
    %                   where each entry is 1,2,4 to indicate 300MHz, 75MHz or 18.75MHz respectively 
    %                   (1,2,4 chosen since there are e.g. 2x the number of logical stations as stations for 75MHz type stations).
    %  .offsets       : Nx2 array, where N is the number of logicalIDs. First column is the offset East from the telescope center,
    %                   second column is the offset north.
    %
    % arrayConfigFull.stations has fields:
    %  "stationType" : one of "300MHz", "75MHz", or "18.75MHz", copy of the input setting.
    %  "logicalIDs" : same as stationsFull.logicalIDs. List of all logical IDs with this stationType.
    %  "stationIDs" : as per stationsFull.stationIDs.
    %  "substationIDs" : as per stationsFull.substationIDs.
    %
    stationsFull.latitude = stations.latitude;
    stationsFull.longitude = stations.longitude;
    stationsFull.altitude = stations.altitude;
    stationsFull.logicalIDs = [];
    stationsFull.stationIDs = [];
    stationsFull.substationIDs = [];
    stationsFull.stationTypes = [];
    stationsFull.offsets = [];
    totalLogicalIDs = 0;

    if (length(arrayConfig.stations) > 3)
        error('stations field in arrayConfig.txt can have at most 3 elements (one each for 300, 75 and 18.75 MHz station types)');
    end

    for s1 = 1:length(arrayConfig.stations)
        ss1 = arrayConfig.stations(s1);
        stationType_string = ss1.stationType{1};
        if strcmp(stationType_string,'300MHz')
            stationType = 1;
        elseif strcmp(stationType_string,'75MHz')
            stationType = 2;
        elseif strcmp(stationType_string,'18.75MHz')
            stationType = 4;
        else
            error('stations.stationType in arrayConfig.txt must be one of "300MHz","75MHz","18.75MHz"');
        end
        newLogicalIDs = (ss1.logicalIDs:(ss1.logicalIDs + stationType * length(ss1.location) - 1));
        stationsFull.logicalIDs((totalLogicalIDs+1):(totalLogicalIDs + stationType * length(ss1.location))) = newLogicalIDs;
        stationsFull.stationTypes((totalLogicalIDs+1):(totalLogicalIDs + stationType * length(ss1.location))) = stationType;
        if (stationType == 1)
            stationsFull.substationIDs((totalLogicalIDs+1):(totalLogicalIDs + stationType * length(ss1.location))) = 1;

            newStationIDs = zeros(1,length(ss1.location));
            newStationIDs(1:(stationType * length(ss1.location))) = ss1.location(1:(stationType * length(ss1.location)));
            stationsFull.stationIDs((totalLogicalIDs+1):(totalLogicalIDs + stationType * length(ss1.location))) = newStationIDs;
            stationsFull.offsets((totalLogicalIDs+1):(totalLogicalIDs + stationType * length(ss1.location)),1:2) = stations.offsets(ss1.location,:);
        elseif (stationType == 2)
            stationsFull.substationIDs((totalLogicalIDs+1):2:(totalLogicalIDs + stationType * length(ss1.location))) = 1;
            stationsFull.substationIDs((totalLogicalIDs+1+1):2:(totalLogicalIDs + stationType * length(ss1.location))) = 2;

            newStationIDs = zeros(1,2*length(ss1.location));
            newStationIDs(1:2:2*length(ss1.location)) = ss1.location(1:(length(ss1.location)));
            newStationIDs(2:2:2*length(ss1.location)) = ss1.location(1:(length(ss1.location)));

            stationsFull.stationIDs((totalLogicalIDs+1):(totalLogicalIDs + stationType * length(ss1.location))) = newStationIDs;

            stationsFull.offsets((totalLogicalIDs+1):2:(totalLogicalIDs + stationType * length(ss1.location)),1:2) = stations.offsets(ss1.location,:);
            stationsFull.offsets((totalLogicalIDs+1+1):2:(totalLogicalIDs + stationType * length(ss1.location)),1:2) = stations.offsets(ss1.location,:);
        else
            stationsFull.substationIDs((totalLogicalIDs+1):4:(totalLogicalIDs + stationType * length(ss1.location))) = 1;
            stationsFull.substationIDs((totalLogicalIDs+1+1):4:(totalLogicalIDs + stationType * length(ss1.location))) = 2;
            stationsFull.substationIDs((totalLogicalIDs+1+2):4:(totalLogicalIDs + stationType * length(ss1.location))) = 3;
            stationsFull.substationIDs((totalLogicalIDs+1+3):4:(totalLogicalIDs + stationType * length(ss1.location))) = 4;

            newStationIDs = zeros(1,4*length(ss1.location));
            newStationIDs(1:4:4*length(ss1.location)) = ss1.location(1:(length(ss1.location)));
            newStationIDs(2:4:4*length(ss1.location)) = ss1.location(1:(length(ss1.location)));
            newStationIDs(3:4:4*length(ss1.location)) = ss1.location(1:(length(ss1.location)));
            newStationIDs(4:4:4*length(ss1.location)) = ss1.location(1:(length(ss1.location)));
            stationsFull.stationIDs((totalLogicalIDs+1):(totalLogicalIDs + stationType * length(ss1.location))) = newStationIDs;
            stationsFull.offsets((totalLogicalIDs+1):4:(totalLogicalIDs + stationType * length(ss1.location)),1:2) = stations.offsets(ss1.location,:);
            stationsFull.offsets((totalLogicalIDs+1+1):4:(totalLogicalIDs + stationType * length(ss1.location)),1:2) = stations.offsets(ss1.location,:);
            stationsFull.offsets((totalLogicalIDs+1+2):4:(totalLogicalIDs + stationType * length(ss1.location)),1:2) = stations.offsets(ss1.location,:);
            stationsFull.offsets((totalLogicalIDs+1+3):4:(totalLogicalIDs + stationType * length(ss1.location)),1:2) = stations.offsets(ss1.location,:);
        end

        arrayConfigFull.stations(s1).stationType = stationType_string;
        arrayConfigFull.stations(s1).logicalIDs = newLogicalIDs;
        arrayConfigFull.stations(s1).stationIDs = newStationIDs;
        arrayConfigFull.stations(s1).substationIDs = stationsFull.substationIDs((totalLogicalIDs+1):(totalLogicalIDs + stationType * length(ss1.location)));

        totalLogicalIDs = totalLogicalIDs + stationType * length(ss1.location);

    end

    % Change the offsets for substations by 10 meters.
    for logicalID = 1:totalLogicalIDs
        if (stationsFull.stationTypes(logicalID) == 2)
            if (stationsFull.substationIDs(logicalID) == 1)
                stationsFull.offsets(logicalID,1) = stationsFull.offsets(logicalID,1) - 10;
            else
                stationsFull.offsets(logicalID,1) = stationsFull.offsets(logicalID,1) + 10;
            end
        elseif (stationsFull.stationTypes(logicalID) == 4)
            if (stationsFull.substationIDs(logicalID) == 1)
                stationsFull.offsets(logicalID,1) = stationsFull.offsets(logicalID,1) - 10;
                stationsFull.offsets(logicalID,2) = stationsFull.offsets(logicalID,2) - 10;
            elseif (stationsFull.substationIDs(logicalID) == 2)
                stationsFull.offsets(logicalID,1) = stationsFull.offsets(logicalID,1) + 10;
                stationsFull.offsets(logicalID,2) = stationsFull.offsets(logicalID,2) - 10;
            elseif (stationsFull.substationIDs(logicalID) == 3)
                stationsFull.offsets(logicalID,1) = stationsFull.offsets(logicalID,1) - 10;
                stationsFull.offsets(logicalID,2) = stationsFull.offsets(logicalID,2) + 10;
            else
                stationsFull.offsets(logicalID,1) = stationsFull.offsets(logicalID,1) + 10;
                stationsFull.offsets(logicalID,2) = stationsFull.offsets(logicalID,2) + 10;            
            end
        end
    end
    arrayConfigFull.stationsFull = stationsFull;

    %% -- "subArray" ------------------------------------------------------------ 
    % Get each subarray listed, error check, put in all fields in a consistent format, and insert into arrayConfigFull
    for s1 = 1:length(arrayConfig.subArray)
        sa1 = arrayConfig.subArray(s1);
        % defaults for all fields in the subarray structure
        sa2.index = 0;
        sa2.logicalIDs = 0;
        sa2.stationType = 384;
        % Error check and copy in the data.
        if iscell(sa1) % Matlab may or may not make it a cell array depending on whether all the entries in the array have the same set of fields.
            sa1 = sa1{1};
        end
        if isfield(sa1,'index')
            sa2.index = sa1.index;
        else
            error('Missing field : Each entry in subArray in arrayConfig.txt must have a field "index"');
        end
        % Station type - must be one of 300MHz, 75MHz, 18.75MHz
        if isfield(sa1,'stationType')
            stationType_string = sa1.stationType;
            if strcmp(stationType_string,'300MHz')
                stationType = 1;
            elseif strcmp(stationType_string,'75MHz')
                stationType = 2;
            elseif strcmp(stationType_string,'18.75MHz')
                stationType = 4;
            else
                error('subArray.stationType in arrayConfig.txt must be one of "300MHz","75MHz","18.75MHz"');
            end
            sa2.stationType = sa1.stationType;
        else
            error('Missing field : Each entry in subArray in arrayConfig.txt must have a field "stationType"');
        end

        if isfield(sa1,'logicalIDs')
            if isnumeric(sa1.logicalIDs)
                sa2.logicalIDs = sa1.logicalIDs;
            elseif iscell(sa1.logicalIDs)
                % Only other option should be a cell with a string "all"
                % Get all the logical IDs with this stationType
                for slist = 1:length(arrayConfigFull.stations)
                    if strcmp(arrayConfigFull.stations(slist).stationType,stationType_string)
                        sa2.logicalIDs = arrayConfigFull.stations(slist).logicalIDs;
                    end
                end
                if ~strcmp(sa1.logicalIDs{1},'all')
                    warning(['expected "all" for logicalIDs field in subarray in arrayConfig.txt, but got ' sa1.logicalIDs{1} '. Assuming all stations'])
                end
            else
                % Not recognised.
                error('logicalIDs field in subArray in arrayConfig.txt is neither a vector or a string');
            end
        else
            error('Missing field : Each entry in subArray in arrayConfig.txt must have a field "logicalIDs"');
        end
        sa2.logicalIDs = sort(sa2.logicalIDs);

        % force the index to the subArray to match subArray.index
        arrayConfigFull.subArray(sa2.index) = sa2;
    end

    % check - each logical ID (arrayConfigFull.subarray(x).logicalIDs) can only be assigned to a single subArray.
    IDUsedInArray = zeros(2048,1);
    IDUsedInArrayStationType = zeros(2048,1);
    for s1 = 1:length(arrayConfig.subArray)
        if strcmp(arrayConfigFull.subArray(s1).stationType,'300MHz')
            stype = 1;
        elseif strcmp(arrayConfigFull.subArray(s1).stationType,'75MHz')
            stype = 2;
        else
            stype = 4;
        end
        for id_index = 1:length(arrayConfigFull.subArray(s1).logicalIDs)
            thisID = arrayConfigFull.subArray(s1).logicalIDs(id_index);
            if (thisID > 2048)
                error(['Found logical ID of ' num2str(thisID) '. Logical IDs must be 2048 or less']);
            end
            if ~IDUsedInArray(thisID)
                IDUsedInArray(thisID) = 1;
                IDUsedInArrayStationType(thisID) = stype;
            else
                error(['Logical ID ' num2str(thisID) ' used in more than one subarray']);
            end
        end
    end

    % Check - each logical ID used in an array must be defined in one and only one place in arrayConfigFull.stations(x).logicalIDs
    for s1 = 1:2048
        found = 0;
        if IDUsedInArray(s1)
            % search arrayConfigFull.stations(x).logicalIDs for s1
            for slist = 1:length(arrayConfigFull.stations)
                if any(arrayConfigFull.stations(slist).logicalIDs == s1)
                    %keyboard
                    if ((strcmp(arrayConfigFull.stations(slist).stationType,'300MHz') && (IDUsedInArrayStationType(s1) == 1)) || ...
                        (strcmp(arrayConfigFull.stations(slist).stationType,'75MHz') && (IDUsedInArrayStationType(s1) == 2)) || ...
                        (strcmp(arrayConfigFull.stations(slist).stationType,'18.75MHz') && (IDUsedInArrayStationType(s1) == 4)))
                        found = 1;
                    else
                        error(['Logical station ' num2str(s1) ' has different station type in stations and subArray']);
                    end
                end
            end
            if (~found)
                error(['Logical station ' num2str(s1) ' is used in a subarray but is not defined in stations']);
            end
        end
    end

    %% -- "stationBeam" ------------------------------------------------------------ 
    % Get each station beam listed, error check, put in all fields in a consistent format, and insert into arrayConfigFull
    % 
    for s1 = 1:length(arrayConfig.stationBeam)
        sb1 = arrayConfig.stationBeam(s1);
        % Defaults for all fields in the station Beam structure
        sb2.index = 0;
        sb2.subArray = 0;
        sb2.delayPolynomial = [0 0 0 0];  % Note there is a 4-tap delay polynomial for every station in the subarray.
        sb2.channels = 0;
        sb2.doppler = 0;
        sb2.skyIndex = 0;  % Not part of the LMC configuration, but needed for generating the LFAA data for this station beam.
        % Error check input data and generate full LMC data.
        if iscell(sb1) % Matlab may or may not make it a cell array depending on whether all the entries in the array have the same set of fields.
            sb1 = sb1{1}; % if it was a cell, convert to a structure.
        end

        if isfield(sb1,'index')
            % index must be in the range 1 to 8
            if ((sb1.index < 1) || (sb1.index > 8))
                error('Station Beam index must be in the range 1 to 8');
            else
                sb2.index = sb1.index;
            end
        else
            error('Missing field : Each entry in stationBeam in arrayConfig.txt must have a field "index"');
        end

        if isfield(sb1,'subArray')
            % subarray must be in the range 1 to 16
            if ((sb1.subArray < 1) || (sb1.subArray > 16))
                error('Station Beam subArray must be in the range 1 to 16');
            else
                sb2.subArray = sb1.subArray;
            end
        else
            error('Missing field : Each entry in stationBeam in arrayConfig.txt must have a field "subArray"');
        end

        if isfield(sb1,'skyIndex')

            sb2.skyIndex = sb1.skyIndex;

            % Look up skyIndex in sky structure
            skyIndex = 0;
            for sk1 = 1:length(sky)
                if ((sky(sk1).index == sb1.skyIndex) && (sky(sk1).subIndex == 0))
                    skyIndex = sk1;
                end
            end
            if (skyIndex == 0)
                error('stationBeam.skyIndex in arrayConfig.txt references a non-existent entry in sky definition (sky.txt)');
            end

            delayPolys = getDelayFunctions(stationsFull,sky(skyIndex),modelConfig.time);
            % Only include the polynomials which relate to stations used in this subarray
            polyIndexList = [];
            for id = 1:length(arrayConfigFull.subArray(sb2.subArray).logicalIDs)
                % logicalID steps through the logical IDs of stations in this subarray
                logicalID = arrayConfigFull.subArray(sb2.subArray).logicalIDs(id);
                % find where this is in stationsFull.logicalIDs, as this will be the location in delayPolys
                polyIndexList = [polyIndexList find(stationsFull.logicalIDs == logicalID)];
            end
            sb2.delayPolynomial = delayPolys.poly(polyIndexList,:);
        else
            error('Missing field : Each entry in stationBeam in arrayConfig.txt must have a field "skyIndex"');
        end

        if isfield(sb1,'channels')
            % Channels field for the stationBeam. There are two special cases-
            %  0 - Full bandwidth. For 75 and 18.75 MHz stations, a second value is required to select the starting channel
            %  1 - Block of bandwidth. Two further numbers are required, to specify the starting channel and the size of the block.
            %  2 - Select the minimal set of coarse channels required to include the data from the sky that this beam is looking at.
            % Otherwise, this is a list of channels in the range 65 to 448, which must be in groups of 8 (except if this involves going past channel 448).
            if (sb1.channels(1) == 0)
                % Full bandwidth. Check station type for this subarray
                if strcmp(arrayConfigFull.subArray(sb1.subArray).stationType{1},'300MHz')
                    % Generate all 384 channels
                    sb2.channels = (65:448);
                elseif strcmp(arrayConfigFull.subArray(sb1.subArray).stationType{1},'75MHz')
                    if (length(sb1.channels) ~= 2)
                        error('When 0 is used for channels field in stationBeam for a 75MHz subarray in arrayConfig.txt, a second number is required to specify the starting channel');
                    else
                        if ((sb1.channels(2) < 65) || (sb1.channels(2) > 448))
                            error('When 0 is used for the first value in channels field in stationBeam, the second value must be in the range 65 to 448 inclusive');
                        else
                            if (sb1.channels(2) <= (448 - 96))
                                sb2.channels = (sb1.channels(2):(sb1.channels(2) + 95));
                            else
                                sb2.channels = (sb1.channels(2):448);
                            end
                        end
                    end
                else
                    % 18.75MHz subArray.
                    if (length(sb1.channels) ~= 2)
                        error('When 0 is used for channels field in stationBeam for an 18.75MHz subarray in arrayConfig.txt, a second number is required to specify the starting channel');
                    else
                        if ((sb1.channels(2) < 65) || (sb1.channels(2) > 448))
                            error('When 0 is used for the first value in channels field in stationBeam, the second value must be in the range 65 to 448 inclusive');
                        else
                            if (sb1.channels(2) <= (448 - 24))
                                sb2.channels = (sb1.channels(2):(sb1.channels(2) + 23));
                            else
                                sb2.channels = (sb1.channels(2):448);
                            end
                        end
                    end
                end
            elseif (sb1.channels(1) == 1)
                % Single contiguous block of channels
                if (length(sb1.channels) ~= 3)
                    error('When 1 is used for the channels field in stationBeam, exactly two more values must be provided to specify the starting channel and the number of channels');
                else
                    sb2.channels = (sb1.channels(2):(sb1.channels(2) + sb1.channels(3) - 1));
                    if ((max(sb2.channels) < 448) && (mod(length(sb2.channels),8) ~= 0))
                        error(['Set of parameters ' num2str(sb2.channels) ' for channels in stationBeam in arrayConfig.txt is not a multiple of 8 channels']);
                    end
                    if ((min(sb2.channels) < 65) || (max(sb2.channels > 448)))
                        error(['Set of parameters ' num2str(sb2.channels) ' for channels in stationBeam in arrayConfig.txt give channels outside the range 65 to 448']);
                    end
                end
            elseif (sb1.channels(1) == 2)
                % Get the minimal set of coarse channels needed to match the sky.
                % Look up skyIndex in sky structure
                skyIndex = 0;
                for sk1 = 1:length(sky)
                    if ((sky(sk1).index == sb1.skyIndex) && (sky(sk1).subIndex == 0))
                        skyIndex = sk1;
                    end
                end
                if (skyIndex == 0)
                    error('stationBeam.channels in arrayConfig.txt references a non-existent entry in sky definition (sky.txt)');
                end
                % min and max frequency in Hz
                skyMinFrequency = sky(skyIndex).coarse * 781250 + sky(skyIndex).fine * (32/27)*781250/4096 - (sky(skyIndex).BW/2) * (32/27)*781250/4096;
                skyMaxFrequency = sky(skyIndex).coarse * 781250 + sky(skyIndex).fine * (32/27)*781250/4096 + (sky(skyIndex).BW/2) * (32/27)*781250/4096;
                % min and max frequency by coarse channel
                skyMinCoarse = round(skyMinFrequency/781250);
                skyMaxCoarse = round(skyMaxFrequency/781250);
                totalCoarse = skyMaxCoarse - skyMinCoarse + 1;
                totalCoarse = ceil(totalCoarse/8) * 8;
                maxCoarse = skyMinCoarse + totalCoarse - 1;
                if (maxCoarse > 448)
                    maxCoarse = 448;
                end
                sb2.channels = (skyMinCoarse:maxCoarse);
            else
                % Check the values are in contiguous blocks of 8, are in the range 65:448 and do not exceed the total bandwidth allowed.
                sb2.channels = sb1.channels;
                for k1 = 1:length(sb2.channels)
                    if (mod(k1-1,8) ~= 0) 
                        % This is not the start of a block of 8 channels, so this channel must be one more than the previous channel
                        if ((sb2.channels(k1-1) + 1) ~= sb2.channels(k1))
                            error('stationBeam.channels in arrayConfig.txt must be blocks of 8 continuous coarse channels');
                        end
                    end
                    if ((sb2.channels(k1) < 65) || (sb2.channels(k1) > 448))
                        error('stationBeam.channels in arrayConfig.txt must be in the range 65 to 448 inclusive');
                    end
                end
                if strcmp(arrayConfigFull.subArray(sb1.subArray).stationType{1},'75MHz')
                    if (length(sb2.channels) > 96)
                        error('Too many channels specified in stationBeam.channels in arrayConfig.txt (maximum is 96 for 75MHz station type)');
                    end
                end
                if strcmp(arrayConfigFull.subArray(sb1.subArray).stationType{1},'18.75MHz')
                    if (length(sb2.channels) > 24)
                        error('Too many channels specified in stationBeam.channels in arrayConfig.txt (maximum is 24 for 18.75MHz station type)');
                    end
                end
                if (sb2.channels(end) ~= 448)
                    % then the number of channels must be a multiple of 8
                    if (mod(length(sb2.channels),8) ~= 0)
                        error('Number of channels specified in stationBeam.channels in arrayConfig.txt must be a multiple of 8 when the last channel is not 448');
                    end
                end
            end

        else
            error('Missing field : Each entry in stationBeam in arrayConfig.txt must have a field "channels"');
        end

        if isfield(sb1,'doppler')
            sb2.doppler = sb1.doppler;
            if (length(sb2.doppler) ~= length(sb2.channels))
                error('stationBeam.doppler must be an array with one element for each channel specified in stationBeam.channels');
            end
        end

        arrayConfigFull.stationBeam(s1) = sb2;
    end

    %% -- "correlator" ------------------------------------------------------------ 
    % Parse the setup for the correlator, and insert into arrayConfigFull
    % 
    if isfield(arrayConfig,'correlator')
        for c1 = 1:length(arrayConfig.correlator)
            cor1 = arrayConfig.correlator(c1);
            if iscell(cor1)
                cor1 = cor1{1};
            end
            % define default values for all the fields in this correlator setup
            cor2.subarray = 0;
            cor2.mode = 'SPECTRAL_LINE';
            cor2.stationBeam = 0;
            cor2.staticRFI = 0; % 0 = false = static RFI detection disabled.
            cor2.dynamicRFI = 0; % 0 = false = dynamic RFI detection disabled.
            cor2.RFIexcision = [];
            cor2.zoom.BW = 384;
            cor2.zoom.startChannel = 65;

            % copy from cor1 to cor2
            if (isfield(cor1,'subArray') && isfield(cor1,'stationBeam'))
                cor2.subArray = cor1.subArray;
                cor2.stationBeam = cor1.stationBeam;
                % check this subArray/stationBeam combination exists
                found = 0;
                for b1 = 1:length(arrayConfigFull.stationBeam)
                    if ((arrayConfigFull.stationBeam(b1).index == cor2.stationBeam) && (arrayConfigFull.stationBeam(b1).subArray == cor2.subArray))
                        found = 1;
                    end
                end
                if (~found)
                    error('Did not find a station beam matching the values specified in correlator.subArray and correlator.stationBeam in arrayConfig.txt');
                end
            else
                error('Each correlator entry in arrayConfig.txt must have fields correlator.subArray and correlator.stationBeam');
            end

            if isfield(cor1,'mode')
                if (strcmp(cor1.mode,'SPECTRAL_LINE') || strcmp(cor1.mode,'ZOOM'))
                    cor2.mode = cor1.mode;
                else
                    error('correlator.mode in arrayConfig.txt must be either "SPECTRAL_LINE" or "ZOOM"');
                end
            else
                error('Each correlator entry in arrayConfig.txt must have a mode field (to set either "SPECTRAL_LINE" or "ZOOM")');
            end

            if isfield(cor1,'staticRFI')
                if strcmpi(cor1.staticRFI,'TRUE')
                    cor2.staticRFI = 1;
                else
                    cor2.staticRFI = 0;
                end
            end

            if isfield(cor1,'dynamicRFI')
                if strcmpi(cor1.dynamicRFI,'TRUE')
                    cor2.staticRFI = 1;
                else
                    cor2.staticRFI = 0;
                end
            end

            if isfield(cor1,'RFIexcision')
                cor2.RFIexcision = cor1.RFIexcision;
            end

            if isfield(cor1,'zoom')
                cor2.zoom.BW = cor1.zoom.BW;
                cor2.zoom.startChannel = cor1.zoom.startChannel;
            end

            % Insert into arrayConfigFull
            arrayConfigFull.correlator(c1) = cor2;
        end
    end

    %% -- "PSSBeam" ------------------------------------------------------------
    % 

    beamsTotal = 0;
    if isfield(arrayConfig,'PSSBeam')
        for p1 = 1:length(arrayConfig.PSSBeam)
            pss1 = arrayConfig.PSSBeam(p1);
            if iscell(pss1)
                pss1 = pss1{1};
            end

            % Check for unrecognized field names
            allnames = fieldnames(pss1);
            for a1 = 1:length(allnames)
                f1 = allnames{a1};
                if (~strcmp(f1,'stationBeam') && ~strcmp(f1,'fill') && ~strcmp(f1,'skySubIndex') && ~strcmp(f1,'subArray') && ~strcmp(f1,'frequency') && ~strcmp(f1,'stationBeam') && ...
                    ~strcmp(f1,'weights') && ~strcmp(f1,'jones') && ~strcmp(f1,'IP') && ~strcmp(f1,'MAC') && ~strcmp(f1,'port') && ~strcmp(f1,'enable') && ~strcmp(f1,'staticRFI') && ...
                    ~strcmp(f1,'dynamicRFI') && ~strcmp(f1,'RFIExcision') && ~strcmp(f1,'index'))
                    error(['Unrecognised field in PSSBeam in arrayConfig.txt, ' f1]);
                end
            end

            % Find the station beam this PSSBeam is using, and then find the sky index for that station beam.
            if ~isfield(pss1,'stationBeam')
                error('Each PSSBeam in arrayConfig.txt must have a field stationBeam');
            end
            stationBeamIndex = 0;
            for b1 = 1:length(arrayConfigFull.stationBeam)
                if ((arrayConfigFull.stationBeam(b1).index == pss1.stationBeam) && (arrayConfigFull.stationBeam(b1).subArray == pss1.subArray))
                    stationBeamIndex = b1;
                    skyIndex = arrayConfigFull.stationBeam(stationBeamIndex).skyIndex;                
                end
            end

            if (stationBeamIndex == 0)
                error('PSSBeam.stationBeam is not defined in stationBeam in arrayConfig.txt');
            end

            % Copy in values from pss1. First determine if we are generating multiple beams.
            if isfield(pss1,'fill')
               if strcmpi(pss1.fill,'TRUE')
                   % Find the number of sky sub indexes for this skyIndex
                   beams = 0;
                   for s1 = 1:length(sky)
                       if (sky(s1).index == skyIndex)
                           beams = beams + 1;
                       end
                   end
               else
                   beams = 1;
               end
            else
                beams = 1;
            end

            for b1 = 1:beams
                % Define defaults for this PSSBeam
                pss2.index = pss1.index + b1 - 1;
                pss2.skySubIndex = 0;  % Record kept but not a part of LMC configuration
                pss2.subArray = 0;
                pss2.delayPolynomial = [0, 0, 0, 0];  % There is a delay polynomial for each logical station which contributes to the stationBeam being used for this PSS beam. It is an offset from the stationBeam polynomial.
                pss2.frequency = 0;
                pss2.stationBeam = 0;
                pss2.weights = 0;
                pss2.jones = 0;
                pss2.IP = [0,0,0,0];
                pss2.MAC = [0,0,0,0,0,0];
                pss2.port = 2000;
                pss2.enable = 1;
                pss2.staticRFI = 0;  % 0 = false;
                pss2.dynamicRFI = 0;   % 0 = false;
                pss2.RFIExcision = [];           
                % Go through and replace defaults with data.

                % .skySubIndex
                if (beams == 1)
                    if isfield(pss1,'skySubIndex')
                        pss2.skySubIndex = pss1.skySubIndex;
                    else
                        error('Each PSSBeam in arrayConfig.txt must have a field skySubIndex');
                    end
                else
                    pss2.skySubIndex = b1;
                end

                % .subArray
                if isfield(pss1,'subArray')
                    pss2.subArray = pss1.subArray;
                else
                    error('Each PSSBeam in arrayConfig.txt must have a field subArray');
                end

                % .delayPolynomial 
                % Get the delay polynomial for this subindex
                if (pss2.skySubIndex == -1)
                    % boresight, so delay polynomial is just 0.
                    pss2.delayPolynomial = zeros(length(arrayConfigFull.subArray(pss2.subArray).logicalIDs),4);
                else
                    % Get the sky coordinates for this subIndex
                    skyLoc = 0;
                    skyLocBoresight = 0;
                    for s1 = 1:length(sky)
                        if ((sky(s1).index == skyIndex) && (sky(s1).subIndex == pss1.skySubIndex))
                            skyLoc = s1;
                        end
                        if ((sky(s1).index == skyIndex) && (sky(s1).subIndex == 0))
                            skyLocBoresight = s1;
                        end
                    end
                    if (skyLoc == 0)
                        error('Failed to find an entry in sky with required value for sky.subIndex in PSS beam configuration (while processing arrayConfig.txt)');
                    end
                    % Actually want the difference from arrayConfigFull.stationBeam(stationBeamIndex).delayPolynomial, which is derived by pointing at skyLocBoresight
                    delayPolysBoresight = getDelayFunctions(stationsFull,sky(skyLocBoresight),modelConfig.time);
                    delayPolys = getDelayFunctions(stationsFull,sky(skyLoc),modelConfig.time);
                    polyDiff = getDelayPolyDiff(delayPolysBoresight,delayPolys,modelConfig.time,10);

                    % Only include the polynomials which relate to stations used in this subarray
                    polyIndexList = [];
                    for id = 1:length(arrayConfigFull.subArray(pss2.subArray).logicalIDs)
                        % logicalID steps through the logical IDs of stations in this subarray
                        logicalID = arrayConfigFull.subArray(pss2.subArray).logicalIDs(id);
                        % find where this is in stationsFull.logicalIDs, as this will be the location in delayPolys
                        polyIndexList = [polyIndexList find(stationsFull.logicalIDs == logicalID)];
                    end

                    pss2.delayPolynomial = polyDiff(polyIndexList,:);
                end

                %.frequency
                if isfield(pss1,'frequency')
                    pss2.frequency = pss1.frequency;
                else
                    error('Each PSSBeam in arrayConfig.txt must have a field frequency');
                end

                %.stationBeam 
                if isfield(pss1,'stationBeam')
                    pss2.stationBeam = pss1.stationBeam;  % Note we have already checked that this exists.
                else
                    error('Each PSSBeam in arrayConfig.txt must have a field stationBeam');
                end

                %.weights
                if isfield(pss1,'weights')
                    if (length(arrayConfigFull.subArray(pss2.subArray).logicalIDs) ~= length(pss1.weights))
                        error('The length of the weights array in PSSBeam in arrayConfig.txt must match the number of logical stations in the subarray');
                    end
                    if ((max(pss1.weights) > 1) || (min(pss1.weights) < 0))
                        error('Each weight in a PSSBeam in arrayConfig.txt must be in the range 0 to 1');
                    end
                    pss2.weights = pss1.weights;
                else
                    pss2.weights = ones(1,length(arrayConfigFull.subArray(pss2.subArray).logicalIDs));
                end

                %.jones
                % must be an array with dimensions (stations)x(coarse channels)x2x2
                if isfield(pss1,'jones')
                    % Summary versions can be in the following formats
                    %  -  2x2 : single matrix replicated for all stations and coarse channels
                    %  -  (stations)x2x2 : different matrices for each station, replicated across all channels.
                    %  -  (stations)x(coarse channels)x2x2 : Full array.
                    %  -  string 'sky' : Inverse of matrices specified for the object in the sky.
                    warning('Jones matrices not supported');
                else

                end

                %.IP = [0,0,0,0];
                if isfield(pss1,'IP')
                    pss2.IP = pss1.IP;
                else
                    pss2.IP = [0,0,0,floor(pss2.index/3)];
                end

                %.MAC = [0,0,0,0,0,0];
                if isfield(pss1,'MAC')
                    pss2.MAC = pss1.MAC;
                else
                    pss2.MAC = [0,0,0,0,0,floor(pss2.index/3)];
                end

                %.port = 2000;
                if isfield(pss1,'port')
                    pss2.port = pss1.port;
                else
                    pss2.port = 2000 + 2*pss2.index;
                end

                %.enable = 1;
                if isfield(pss1,'enable')
                    if isstring(pss1.enable)
                        if strcmpi(pss1.enable,'true')
                            pss2.enable = 1;
                        else
                            pss2.enable = 0;
                        end
                    else
                        pss2.enable = pss1.enable;
                    end
                else
                    pss2.enable = 1;
                end

                %.staticRFI = 0;  % 0 = false;
                if isfield(pss1,'staticRFI')
                    if isstring(pss1.staticRFI)
                        if strcmpi(pss1.staticRFI,'true')
                            pss2.staticRFI = 1;
                        else
                            pss2.staticRFI = 0;
                        end
                    else
                        pss2.staticRFI = pss1.staticRFI;
                    end
                else
                    pss2.staticRFI = 1;
                end

                %.dynamicRFI = 0;   % 0 = false;
                if isfield(pss1,'dynamicRFI')
                    if isstring(pss1.dynamicRFI)
                        if strcmpi(pss1.dynamicRFI,'true')
                            pss2.dynamicRFI = 1;
                        else
                            pss2.dynamicRFI = 0;
                        end
                    else
                        pss2.dynamicRFI = pss1.dynamicRFI;
                    end
                else
                    pss2.dynamicRFI = 1;
                end            


                %.RFIexcision = []; 
                if isfield(pss1,'RFIExcision')
                    warning('RFIExcision is not supported yet.');
                end

                beamsTotal = beamsTotal + 1;
                arrayConfigFull.PSSBeam(beamsTotal) = pss2;
            end

        end
    end

    if (beamsTotal > 500)
        error(['Found ' num2str(beamsTotal) ' PSS beams. Maximum allowed is 16']);
    end

    %% -- "PSTBeam" ------------------------------------------------------------
    % 

    beamsTotal = 0;
    if isfield(arrayConfig,'PSTBeam')
        for p1 = 1:length(arrayConfig.PSTBeam)
            pst1 = arrayConfig.PSTBeam(p1);
            if iscell(pst1)
                pst1 = pst1{1};
            end

            % Check for unrecognized field names
            allnames = fieldnames(pst1);
            for a1 = 1:length(allnames)
                f1 = allnames{a1};
                if (~strcmp(f1,'index') && ~strcmp(f1,'skySubIndex') && ~strcmp(f1,'subArray') && ~strcmp(f1,'stationBeam') && ~strcmp(f1,'weights') && ...
                    ~strcmp(f1,'jones') && ~strcmp(f1,'enable') &&  ~strcmp(f1,'IP') && ~strcmp(f1,'MAC') && ~strcmp(f1,'port')  && ~strcmp(f1,'staticRFI') && ...
                    ~strcmp(f1,'dynamicRFI') && ~strcmp(f1,'RFIExcision'))
                    error(['Unrecognised field in PSTBeam in arrayConfig.txt, ' f1]);
                end
            end

            % Find the station beam this PSTBeam is using, and then find the sky index for that station beam.
            matchingBeams = 0;
            if ~isfield(pst1,'stationBeam')
                % Not defined, so pick the station beam based on the subarray (there should only be one stationBeam for a PST beam).
                stationBeamIndex = 0;
                for b1 = 1:length(arrayConfigFull.stationBeam)
                    if (arrayConfigFull.stationBeam(b1).subArray == pst1.subArray)
                        stationBeamIndex = b1;
                        skyIndex = arrayConfigFull.stationBeam(stationBeamIndex).skyIndex;
                        pst1.stationBeam = arrayConfigFull.stationBeam(stationBeamIndex).index;
                        matchingBeams = matchingBeams + 1;
                    end
                end
            else
                stationBeamIndex = 0;
                for b1 = 1:length(arrayConfigFull.stationBeam)
                    if ((arrayConfigFull.stationBeam(b1).index == pst1.stationBeam) && (arrayConfigFull.stationBeam(b1).subArray == pst1.subArray))
                        stationBeamIndex = b1;
                        skyIndex = arrayConfigFull.stationBeam(stationBeamIndex).skyIndex;
                        matchingBeams = matchingBeams + 1;
                    end
                end

                if (stationBeamIndex == 0)
                    error('PSTBeam.stationBeam is not defined in stationBeam in arrayConfig.txt');
                end
            end
            if (matchingBeams ~= 1)
                error(['Problem with PSTBeam ' num2str(pst1.index) '. Each PSTBeam must be assigned to a subarray with a single stationBeam']);
            end

            % Copy in values from pst1.

            % Define defaults for this PSTBeam
            pst2.index = pst1.index + b1 - 1;
            pst2.skySubIndex = 0;  % Record kept but not a part of LMC configuration
            pst2.subArray = 0;
            pst2.delayPolynomial = [0, 0, 0, 0];  % There is a delay polynomial for each logical station which contributes to the stationBeam being used for this PSS beam. It is an offset from the stationBeam polynomial.
            pst2.stationBeam = 0;
            pst2.weights = 0;
            pst2.jones = 0;
            pst2.IP = [0,0,0,0];
            pst2.MAC = [0,0,0,0,0,0];
            pst2.port = 2000;
            pst2.enable = 1;
            pst2.staticRFI = 0;  % 0 = false;
            pst2.dynamicRFI = 0;   % 0 = false;
            pst2.RFIExcision = [];           
            % Go through and replace defaults with data.

            % .skySubIndex
            if isfield(pst1,'skySubIndex')
                pst2.skySubIndex = pst1.skySubIndex;
            else
                error('Each PSTBeam in arrayConfig.txt must have a field skySubIndex');
            end

            % .subArray
            if isfield(pst1,'subArray')
                pst2.subArray = pst1.subArray;
            else
                error('Each PSTBeam in arrayConfig.txt must have a field subArray');
            end

            % .delayPolynomial 
            % Get the delay polynomial for this subindex
            if (pst2.skySubIndex == -1)
                % boresight, so delay polynomial is just 0.
                pst2.delayPolynomial = zeros(length(arrayConfigFull.subArray(pst2.subArray).logicalIDs),4);
            else
                % Get the sky coordinates for this subIndex
                skyLoc = 0;
                skyLocBoresight = 0;
                for s1 = 1:length(sky)
                    if ((sky(s1).index == skyIndex) && (sky(s1).subIndex == pst1.skySubIndex))
                        skyLoc = s1;
                    end
                    if ((sky(s1).index == skyIndex) && (sky(s1).subIndex == 0))
                        skyLocBoresight = s1;
                    end
                end
                if (skyLoc == 0)
                    error('Failed to find an entry in sky with required value for sky.subIndex in PST beam configuration (while processing arrayConfig.txt)');
                end
                % Actually want the difference from arrayConfigFull.stationBeam(stationBeamIndex).delayPolynomial, which is derived by pointing at skyLocBoresight
                delayPolysBoresight = getDelayFunctions(stationsFull,sky(skyLocBoresight),modelConfig.time);
                delayPolys = getDelayFunctions(stationsFull,sky(skyLoc),modelConfig.time);
                polyDiff = getDelayPolyDiff(delayPolysBoresight,delayPolys,modelConfig.time,10);

                % Only include the polynomials which relate to stations used in this subarray
                polyIndexList = [];
                for id = 1:length(arrayConfigFull.subArray(pst2.subArray).logicalIDs)
                    % logicalID steps through the logical IDs of stations in this subarray
                    logicalID = arrayConfigFull.subArray(pst2.subArray).logicalIDs(id);
                    % find where this is in stationsFull.logicalIDs, as this will be the location in delayPolys
                    polyIndexList = [polyIndexList find(stationsFull.logicalIDs == logicalID)];
                end
                pst2.delayPolynomial = polyDiff(polyIndexList,:);

            end

            %.stationBeam
            if isfield(pst1,'stationBeam')
                pst2.stationBeam = pst1.stationBeam;  % Note we have already checked that this exists.
            else
                error('Each PSTBeam in arrayConfig.txt must have a field stationBeam');
            end

            %.weights
            if isfield(pst1,'weights')
                if (length(arrayConfigFull.subArray(pst2.subArray).logicalIDs) ~= length(pst1.weights))
                    error('The length of the weights array in PSTBeam in arrayConfig.txt must match the number of logical stations in the subarray');
                end
                if ((max(pst1.weights) > 1) || (min(pst1.weights) < 0))
                    error('Each weight in a PSTBeam in arrayConfig.txt must be in the range 0 to 1');
                end
                pst2.weights = pst1.weights;
            else
                pst2.weights = ones(1,length(arrayConfigFull.subArray(pst2.subArray).logicalIDs));
            end

            %.jones
            % must be an array with dimensions (stations)x(coarse channels)x2x2
            if isfield(pst1,'jones')
                % Summary versions can be in the following formats
                %  -  2x2 : single matrix replicated for all stations and coarse channels
                %  -  (stations)x2x2 : different matrices for each station, replicated across all channels.
                %  -  (stations)x(coarse channels)x2x2 : Full array.
                %  -  string 'sky' : Inverse of matrices specified for the object in the sky.
                warning('Jones matrices not supported for PST at present');
            end

            %.IP = [0,0,0,0];
            if isfield(pst1,'IP')
                pst2.IP = pst1.IP;
            else
                pst2.IP = [0,0,0,pst2.index];
            end

            %.MAC = [0,0,0,0,0,0];
            if isfield(pst1,'MAC')
                pst2.MAC = pst1.MAC;
            else
                pst2.MAC = [0,0,0,0,0,pst2.index];
            end

            %.port = 2000;
            if isfield(pst1,'port')
                pst2.port = pst1.port;
            else
                pst2.port = 2000 + 2*pst2.index;
            end

            %.enable = 1;
            if isfield(pst1,'enable')
                if isstring(pst1.enable)
                    if strcmpi(pst1.enable,'true')
                        pst2.enable = 1;
                    else
                        pst2.enable = 0;
                    end
                else
                    pst2.enable = pst1.enable;
                end
            else
                pst2.enable = 1;
            end

            %.staticRFI = 0;  % 0 = false;
            if isfield(pst1,'staticRFI')
                if isstring(pst1.staticRFI)
                    if strcmpi(pst1.staticRFI,'true')
                        pst2.staticRFI = 1;
                    else
                        pst2.staticRFI = 0;
                    end
                else
                    pst2.staticRFI = pst1.staticRFI;
                end
            else
                pst2.staticRFI = 1;
            end

            %.dynamicRFI = 0;   % 0 = false;
            if isfield(pst1,'dynamicRFI')
                if isstring(pst1.dynamicRFI)
                    if strcmpi(pst1.dynamicRFI,'true')
                        pst2.dynamicRFI = 1;
                    else
                        pst2.dynamicRFI = 0;
                    end
                else
                    pst2.dynamicRFI = pst1.dynamicRFI;
                end
            else
                pst2.dynamicRFI = 1;
            end


            %.RFIexcision = []; 
            if isfield(pst1,'RFIExcision')
                warning('RFIExcision is not supported yet.');
            end

            beamsTotal = beamsTotal + 1;
            arrayConfigFull.PSTBeam(beamsTotal) = pst2;

        end
    end

    if (beamsTotal > 16)
        error(['Found ' num2str(beamsTotal) ' PST beams. Maximum allowed is 16']);
    end

    %% -- "VLBIBeam" ------------------------------------------------------------
    % 
    beamsTotal = 0;
    if isfield(arrayConfig,'VLBIBeam')

        for p1 = 1:length(arrayConfig.VLBIBeam)
            vlbi1 = arrayConfig.VLBIBeam(p1);
            if iscell(vlbi1)
                vlbi1 = vlbi1{1};
            end

            % Check for unrecognized field names
            allnames = fieldnames(vlbi1);
            for a1 = 1:length(allnames)
                f1 = allnames{a1};
                if (~strcmp(f1,'index') && ~strcmp(f1,'skySubIndex') && ~strcmp(f1,'subArray') && ~strcmp(f1,'stationBeam') && ~strcmp(f1,'weights') && ...
                    ~strcmp(f1,'BW') && ~strcmp(f1,'overSampling') &&  ~strcmp(f1,'subbandEnable') && ~strcmp(f1,'frequency') && ~strcmp(f1,'address')  && ~strcmp(f1,'staticRFI') && ...
                    ~strcmp(f1,'dynamicRFI') && ~strcmp(f1,'bits') && ~strcmp(f1,'jones') && ~strcmp(f1,'RFIExcision'))
                    error(['Unrecognised field in VLBIBeam in arrayConfig.txt, ' f1]);
                end
            end

            % Find the station beam this VLBIBeam is using, and then find the sky index for that station beam.
            matchingBeams = 0;
            if ~isfield(vlbi1,'stationBeam')
                % Not defined, so pick the station beam based on the subarray (there should only be one stationBeam for a PST beam).
                stationBeamIndex = 0;
                for b1 = 1:length(arrayConfigFull.stationBeam)
                    if (arrayConfigFull.stationBeam(b1).subArray == vlbi1.subArray)
                        stationBeamIndex = b1;
                        skyIndex = arrayConfigFull.stationBeam(stationBeamIndex).skyIndex;
                        vlbi1.stationBeam = arrayConfigFull.stationBeam(stationBeamIndex).index;
                        matchingBeams = matchingBeams + 1;
                    end
                end
            else
                stationBeamIndex = 0;
                for b1 = 1:length(arrayConfigFull.stationBeam)
                    if ((arrayConfigFull.stationBeam(b1).index == vlbi1.stationBeam) && (arrayConfigFull.stationBeam(b1).subArray == vlbi1.subArray))
                        stationBeamIndex = b1;
                        skyIndex = arrayConfigFull.stationBeam(stationBeamIndex).skyIndex;
                        matchingBeams = matchingBeams + 1;
                    end
                end

                if (stationBeamIndex == 0)
                    error('PSTBeam.stationBeam is not defined in stationBeam in arrayConfig.txt');
                end
            end
            if (matchingBeams ~= 1)
                error(['Problem with VLBIBeam ' num2str(vlbi1.index) '. Each VLBIBeam must be assigned to a subarray with a single stationBeam']);
            end

            % Copy in values from vlbi1.

            % Define defaults for this vlbiBeam
            vlbi2.index = vlbi1.index + b1 - 1;
            vlbi2.skySubIndex = 0;  % Record kept but not a part of LMC configuration
            vlbi2.subArray = 0;
            vlbi2.delayPolynomial = [0, 0, 0, 0];  % There is a delay polynomial for each logical station which contributes to the stationBeam being used for this PSS beam. It is an offset from the stationBeam polynomial.
            vlbi2.stationBeam = 0;
            vlbi2.weights = 0;
            vlbi2.BW = 1;
            vlbi2.oversampling = 0;
            vlbi2.subbandEnable = 1;
            vlbi2.frequency = 0;
            vlbi2.address = [0,0,0,0];
            vlbi2.staticRFI = 0;  % 0 = false;
            vlbi2.dynamicRFI = 0;   % 0 = false;
            vlbi2.bits = 8;
            vlbi2.jones = 0;
            vlbi2.enable = 1;
            vlbi2.RFIExcision = [];           
            % Go through and replace defaults with data.

            % .skySubIndex
            if isfield(vlbi1,'skySubIndex')
                vlbi2.skySubIndex = vlbi1.skySubIndex;
            else
                error('Each VLBIBeam in arrayConfig.txt must have a field skySubIndex');
            end

            % .subArray
            if isfield(vlbi1,'subArray')
                vlbi2.subArray = vlbi1.subArray;
            else
                error('Each VLBIBeam in arrayConfig.txt must have a field subArray');
            end

            % .delayPolynomial 
            % Get the delay polynomial for this subindex
            if (vlbi2.skySubIndex == -1)
                % boresight, so delay polynomial is just 0.
                vlbi2.delayPolynomial = zeros(length(arrayConfigFull.subArray(vlbi2.subArray).logicalIDs),4);
            else
                % Get the sky coordinates for this subIndex
                skyLoc = 0;
                skyLocBoresight = 0;
                for s1 = 1:length(sky)
                    if ((sky(s1).index == skyIndex) && (sky(s1).subIndex == vlbi1.skySubIndex))
                        skyLoc = s1;
                    end
                    if ((sky(s1).index == skyIndex) && (sky(s1).subIndex == 0))
                        skyLocBoresight = s1;
                    end
                end
                if (skyLoc == 0)
                    error('Failed to find an entry in sky with required value for sky.subIndex in PST beam configuration (while processing arrayConfig.txt)');
                end
                % Actually want the difference from arrayConfigFull.stationBeam(stationBeamIndex).delayPolynomial, which is derived by pointing at skyLocBoresight
                delayPolysBoresight = getDelayFunctions(stationsFull,sky(skyLocBoresight),modelConfig.time);
                delayPolys = getDelayFunctions(stationsFull,sky(skyLoc),modelConfig.time);
                polyDiff = getDelayPolyDiff(delayPolysBoresight,delayPolys,modelConfig.time,10);

                % Only include the polynomials which relate to stations used in this subarray
                polyIndexList = [];
                for id = 1:length(arrayConfigFull.subArray(vlbi2.subArray).logicalIDs)
                    % logicalID steps through the logical IDs of stations in this subarray
                    logicalID = arrayConfigFull.subArray(vlbi2.subArray).logicalIDs(id);
                    % find where this is in stationsFull.logicalIDs, as this will be the location in delayPolys
                    polyIndexList = [polyIndexList find(stationsFull.logicalIDs == logicalID)];
                end
                vlbi2.delayPolynomial = polyDiff(polyIndexList,:);

            end

            %.stationBeam
            vlbi2.stationBeam = vlbi1.stationBeam;  % Note we have already checked that this exists.

            %.weights
            if isfield(vlbi1,'weights')
                if (length(arrayConfigFull.subArray(vlbi2.subArray).logicalIDs) ~= length(vlbi1.weights))
                    error('The length of the weights array in VLBIBeam in arrayConfig.txt must match the number of logical stations in the subarray');
                end
                if ((max(vlbi1.weights) > 1) || (min(vlbi1.weights) < 0))
                    error('Each weight in a VLBIBeam in arrayConfig.txt must be in the range 0 to 1');
                end
                vlbi2.weights = vlbi1.weights;
            else
                vlbi2.weights = ones(1,length(arrayConfigFull.subArray(vlbi2.subArray).logicalIDs));
            end

            %.BW = 1 to 256 MHz;
            if isfield(vlbi1,'BW')
                vlbi2.BW = vlbi1.BW;
            end

            %.oversampling = 16 to 32;
            if isfield(vlbi1,'oversampling')
                vlbi2.oversampling = vlbi1.oversampling;
            end

            %.subbandEnable = list of integers in range 1 to 16;
            if isfield(vlbi1,'subbandEnable')
                vlbi2.subbandEnable = vlbi1.subbandEnable;
            end

            %.frequency = 0;
            if isfield(vlbi1,'frequency')
                vlbi2.frequency = vlbi1.frequency;
            end

            %.address = TBD.
            if isfield(vlbi1,'address')
                vlbi2.address = vlbi1.address;
            end

            %.staticRFI = 0;  % 0 = false;
            if isfield(vlbi1,'staticRFI')
                if isstring(vlbi1.staticRFI)
                    if strcmpi(vlbi1.staticRFI,'true')
                        vlbi2.staticRFI = 1;
                    else
                        vlbi2.staticRFI = 0;
                    end
                else
                    vlbi2.staticRFI = vlbi1.staticRFI;
                end
            else
                vlbi2.staticRFI = 0;
            end

            %.dynamicRFI = 0;   % 0 = false;
            if isfield(vlbi1,'dynamicRFI')
                if isstring(vlbi1.dynamicRFI)
                    if strcmpi(vlbi1.dynamicRFI,'true')
                        vlbi2.dynamicRFI = 1;
                    else
                        vlbi2.dynamicRFI = 0;
                    end
                else
                    vlbi2.dynamicRFI = vlbi1.dynamicRFI;
                end
            else
                vlbi2.dynamicRFI = 0;
            end

            %.bits = 2,4 or 8.
            if isfield(vlbi1,'bits')
                vlbi2.bits = vlbi1.bits;
            end
            if ((vlbi2.bits ~= 2) && (vlbi2.bits ~= 4) && (vlbi2.bits ~= 8))
                error('VLBI.bits must be one of 2, 4 or 8');
            end

            %.jones
            % must be an array with dimensions (stations)x(coarse channels)x2x2
            if isfield(vlbi1,'jones')
                % Summary versions can be in the following formats
                %  -  2x2 : single matrix replicated for all stations and coarse channels
                %  -  (stations)x2x2 : different matrices for each station, replicated across all channels.
                %  -  (stations)x(coarse channels)x2x2 : Full array.
                %  -  string 'sky' : Inverse of matrices specified for the object in the sky.
                warning('Jones matrices not supported for VLBI at present');
            end

            %.RFIexcision = []; 
            if isfield(vlbi1,'RFIExcision')
                warning('RFIExcision is not supported yet.');
            end

            beamsTotal = beamsTotal + 1;
            arrayConfigFull.vlbiBeam(beamsTotal) = vlbi2;

        end

    end

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
    %   .txTimes : Time to send each packet in seconds. Dimensions M x 1. Doubles.
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
    keyboard
    for s1 = 1:length(sky)
        % add to the sky structure a list of coarse channels that are non-zero, and a raw data structure
        % which is MxN, where M is the number of samples generated and N is the number of non-zero coarse channels.
        
    end
    
    %keyboard
    % Each FPGA gets data for 2 stations, first step through FPGA
    for fpgaIndex = 1:floor(length(stationList)/2)
        % Each FPGA gets data from up to 8 logical stations (2 stations x 4 substations)
        % Get a list of logical stations for this FPGA.
        
        logicalIDs = [];
        stationIDs = [];   % Each logicalID has a unique combination of stationID and a substationID
        substationIDs = [];
        channelCount = 0;   % count of the number of channels being sent to this FPGA. Can get up to 2*384.
        usedChannelCount = 0;  % count of the number of channels being sent to this FPGA that are non-zero.
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
        
        %keyboard
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
                            % Count of the channels coming from this station
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
                            
                            % UDP checksum calculation; do this later when data is available.
%                             pseudoHdr(1:8) = hdr(27:34); % src and destination IP address
%                             pseudoHdr(9) = 0;            % 0
%                             pseudoHdr(10) = 17;          % protocol
%                             pseudoHdr(11) = hdr(39);     % UDP length field, same as length field in the UDP header.
%                             pseudoHdr(12) = hdr(40);  
%                             pseudoHdr(13:18) = hdr(35:40); % source UDP port, destination UDP port, UDP length.
%                             pseudoHdr16bit = 256 * pseudoHdr(1:2:18) + pseudoHdr(2:2:18); 
%                             hsum = sum(pseudoHdr16bit);
%                             while (hsum > 65535)
%                                 hsum = mod(hsum,65536) + floor(hsum/65536);
%                             end
%                             hsum = 65535 - hsum;
%                             hdr(41) = floor(hsum/256);  % UDP checksum
%                             hdr(42) = mod(hsum,256);
                            
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
                            
                            % Generate data for this logical channel and channel i.e. fpga(fpgaIndex).data(:,usedChannelCount)
                            % Need the sinusoid for the delay rather than the polynomial approximation.
                            for s1 = 1:length(sky)
                                % Find each object in the sky which is in the view of this beam, i.e. sky(s1).index == skyIndex
                                s2 = find(arrayConfigFull.stationsFull.logicalIDs == logicalID);
                                thisStation.latitude = arrayConfigFull.stationsFull.latitude;
                                thisStation.longitude = arrayConfigFull.stationsFull.longitude;
                                thisStation.altitude = arrayConfigFull.stationsFull.altitude;
                                thisStation.offsets = arrayConfigFull.stationsFull.offsets(s2,:);
                                if (sky(s1).index == skyIndex)
                                    delayFunction = getDelayFunctions(thisStation,sky(s1),modelConfig.time);
                                    % reformat delayFunction into the form required by resampleNU.
                                    delayOAFP(1) = delayFunction.offset;
                                    delayOAFP(2) = delayFunction.amplitude;
                                    delayOAFP(3) = delayFunction.rate;
                                    delayOAFP(4) = delayFunction.phase;
                                    resampledPoints = modelConfig.runtime * 408 * 2048; % Number of sample points to generate
                                    % Get the data.
                                    
                                    
                                    % Resample the data for this station
                                    resampleNU(din,1080,channelFrequency,delayOAFP,resampledPoints)
                                    keyboard
                                end
                                
                            end
                        end
                    end
                end
            
%            for channel = 1:something
%                for beamIndex = 1:length(beamList)
                
                
                
                
%            end
            end
        end
        % At this point fpga(fpgaIndex).headers contains a list of headers for one instant in time.
        % Replicate and modify this list so that fpga(fpgaIndex).headers has all the headers for the 
        % full simulation time.
%        for sampleTime = 1:something
%        end
        
        
    end
    
    
    
else
    disp('LFAA data is already up to date in LFAA.mat');
    
end

%% 


%%
% Generate the register settings in the firmware from the LMC configuration


