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
% 
% 1. Station processing:
%     For each FPGA
%        load input data from LFAA.mat
%        do processing (LFAA SPEAD ingest, doppler, coarse splitter)
%        save results to stationX.mat
%  
% 2. Array Processing:
%     For each FPGA
%        load input data from stationX.mat
%        do processing (Coarse corner turn, filterbanks, fine delays)
%        save results to array.mat
%
% 3. Output Processing:
%     For each FPGA
%        load input data from array.mat
%        do processing (packetising for SDP, PSS and PST)
%        save results to output.mat
%        
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
% Station Processing

% Get the LFAA data
for lru = 1:modelConfig.LRUs
    
end

%% -----------------------------------------------------------------------
% Array Processing


%% -----------------------------------------------------------------------
% Output Processing
