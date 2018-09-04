function [modelConfig] = parseModelConfig(rundir)
% Parse 'modelConfig.txt', a json file describing the configuration of the model.
% Checks that the parameters are valid and adds in defaults for any missing parameters.
% Returns the full configuration in the modelConfig structure.
%
%
%

%% Load Model setup
fid = fopen([rundir '/modelConfig.txt']);
model_json = fread(fid,inf);
fclose(fid);
model_json = char(model_json');
modelConfig = jsondecode(model_json);

%% Check compulsory fields
if ~isfield(modelConfig,'SNR')
    error('SNR field is missing from modelConfig.txt');
end

if ~isfield(modelConfig,'time')
    error('time field is missing from modelConfig.txt');
end

if ~isfield(modelConfig,'runtime')
    error('runtime field is missing from modelConfig.txt');
end

if ~isfield(modelConfig,'configuration')
    error('configuration field is missing from modelConfig.txt');
end

if ~isfield(modelConfig,'runCorrelator')
    error('runCorrelator field is missing from modelConfig.txt');
end

if ~isfield(modelConfig,'runPSS')
    error('runPSS field is missing from modelConfig.txt');
end

if ~isfield(modelConfig,'runPST')
    error('runPST field is missing from modelConfig.txt');
end

% Get the number of FPGAs for this configuration
if (modelConfig.configuration == 0) 
    modelConfig.LRUs = 3;
elseif (modelConfig.configuration == 1)
    modelConfig.LRUs = 12;
elseif (modelConfig.configuration == 2)
    modelConfig.LRUs = 36;
elseif (modelConfig.configuration == 3) 
    modelConfig.LRUs = 48;
elseif (modelConfig.configuration == 4)
    modelConfig.LRUs = 144;
elseif (modelConfig.configuration == 5)
    modelConfig.LRUs = 288;
else
    error('configuration field in modelConfig.txt must be in the range 0 to 5');
end

%% Check or insert Optional Fields

if ~isfield(modelConfig,'stationMap')
    modelConfig.stationMap = 1:(2*modelConfig.LRUs);
end

if ~isfield(modelConfig,'keepLFAASPEADDecoder')
    modelConfig.keepLFAASPEADDecoder = 0;
else
    if (modelConfig.keepLFAASPEADDecoder(1) == -1)
        modelConfig.keepLFAASPEADDecoder = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepLFAASPEADDecoder) < -1) || (max(modelConfig.keepLFAASPEADDecoder) > modelConfig.LRUs))
        error('keepLFAASPEADDecoder field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepLFAASPEADDecoder = unique(modelConfig.keepLFAASPEADDecoder);
    end
end

if ~isfield(modelConfig,'keepLocalDoppler')
    modelConfig.keepLocalDoppler = 0;
else
    if (modelConfig.keepLocalDoppler(1) == -1)
        modelConfig.keepLocalDoppler = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepLocalDoppler) < -1) || (max(modelConfig.keepLocalDoppler) > modelConfig.LRUs))
        error('keepLocalDoppler field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepLocalDoppler = unique(modelConfig.keepLocalDoppler);
    end
end

if ~isfield(modelConfig,'keepCoarseSplitter')
    modelConfig.keepCoarseSplitter = 0;
else
    if (modelConfig.keepCoarseSplitter(1) == -1)
        modelConfig.keepCoarseSplitter = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepCoarseSplitter) < -1) || (max(modelConfig.keepCoarseSplitter) > modelConfig.LRUs))
        error('keepCoarseSplitter field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepCoarseSplitter = unique(modelConfig.keepCoarseSplitter);
    end
end

if ~isfield(modelConfig,'keepCoarseCornerTurner')
    modelConfig.keepCoarseCornerTurner = 0;
else
    if (modelConfig.keepCoarseCornerTurner(1) == -1)
        modelConfig.keepCoarseCornerTurner = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepCoarseCornerTurner) < -1) || (max(modelConfig.keepCoarseCornerTurner) > modelConfig.LRUs))
        error('keepCoarseCornerTurner field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepCoarseCornerTurner = unique(modelConfig.keepCoarseCornerTurner);
    end
end

if ~isfield(modelConfig,'keepCorrelatorFilterbank')
    modelConfig.keepCorrelatorFilterbank = 0;
else
    if (modelConfig.keepCorrelatorFilterbank(1) == -1)
        modelConfig.keepCorrelatorFilterbank = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepCorrelatorFilterbank) < -1) || (max(modelConfig.keepCorrelatorFilterbank) > modelConfig.LRUs))
        error('keepCorrelatorFilterbank field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepCorrelatorFilterbank = unique(modelConfig.keepCorrelatorFilterbank);
    end
end

if ~isfield(modelConfig,'keepPSSFilterbank')
    modelConfig.keepPSSFilterbank = 0;
else
    if (modelConfig.keepPSSFilterbank(1) == -1)
        modelConfig.keepPSSFilterbank = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSSFilterbank) < -1) || (max(modelConfig.keepPSSFilterbank) > modelConfig.LRUs))
        error('keepPSSFilterbank field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSSFilterbank = unique(modelConfig.keepPSSFilterbank);
    end
end

if ~isfield(modelConfig,'keepPSTFilterbank')
    modelConfig.keepPSTFilterbank = 0;
else
    if (modelConfig.keepPSTFilterbank(1) == -1)
        modelConfig.keepPSTFilterbank = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSTFilterbank) < -1) || (max(modelConfig.keepPSTFilterbank) > modelConfig.LRUs))
        error('keepPSTFilterbank field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSTFilterbank = unique(modelConfig.keepPSTFilterbank);
    end
end

if ~isfield(modelConfig,'keepCorrelatorFineDelay')
    modelConfig.keepCorrelatorFineDelay = 0;
else
    if (modelConfig.keepCorrelatorFineDelay(1) == -1)
        modelConfig.keepCorrelatorFineDelay = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepCorrelatorFineDelay) < -1) || (max(modelConfig.keepCorrelatorFineDelay) > modelConfig.LRUs))
        error('keepCorrelatorFineDelay field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepCorrelatorFineDelay = unique(modelConfig.keepCorrelatorFineDelay);
    end
end

if ~isfield(modelConfig,'keepPSSFineDelay')
    modelConfig.keepPSSFineDelay = 0;
else
    if (modelConfig.keepPSSFineDelay(1) == -1)
        modelConfig.keepPSSFineDelay = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSSFineDelay) < -1) || (max(modelConfig.keepPSSFineDelay) > modelConfig.LRUs))
        error('keepPSSFineDelay field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSSFineDelay = unique(modelConfig.keepPSSFineDelay);
    end
end

if ~isfield(modelConfig,'keepPSTFineDelay')
    modelConfig.keepPSTFineDelay = 0;
else
    if (modelConfig.keepPSTFineDelay(1) == -1)
        modelConfig.keepPSTFineDelay = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSTFineDelay) < -1) || (max(modelConfig.keepPSTFineDelay) > modelConfig.LRUs))
        error('keepPSTFineDelay field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSTFineDelay = unique(modelConfig.keepPSTFineDelay);
    end
end

if ~isfield(modelConfig,'keepCorrelatorFineCornerTurner')
    modelConfig.keepCorrelatorFineCornerTurner = 0;
else
    if (modelConfig.keepCorrelatorFineCornerTurner(1) == -1)
        modelConfig.keepCorrelatorFineCornerTurner = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepCorrelatorFineCornerTurner) < -1) || (max(modelConfig.keepCorrelatorFineCornerTurner) > modelConfig.LRUs))
        error('keepCorrelatorFineCornerTurner field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepCorrelatorFineCornerTurner = unique(modelConfig.keepCorrelatorFineCornerTurner);
    end
end

if ~isfield(modelConfig,'keepPSSFineCornerTurner')
    modelConfig.keepPSSFineCornerTurner = 0;
else
    if (modelConfig.keepPSSFineCornerTurner(1) == -1)
        modelConfig.keepPSSFineCornerTurner = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSSFineCornerTurner) < -1) || (max(modelConfig.keepPSSFineCornerTurner) > modelConfig.LRUs))
        error('keepPSSFineCornerTurner field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSSFineCornerTurner = unique(modelConfig.keepPSSFineCornerTurner);
    end
end

if ~isfield(modelConfig,'keepPSTFineCornerTurner')
    modelConfig.keepPSTFineCornerTurner = 0;
else
    if (modelConfig.keepPSTFineCornerTurner(1) == -1)
        modelConfig.keepPSTFineCornerTurner = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSTFineCornerTurner) < -1) || (max(modelConfig.keepPSTFineCornerTurner) > modelConfig.LRUs))
        error('keepPSTFineCornerTurner field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSTFineCornerTurner = unique(modelConfig.keepPSTFineCornerTurner);
    end
end

if ~isfield(modelConfig,'keepCorrelator')
    modelConfig.keepCorrelator = 0;
else
    if (modelConfig.keepCorrelator(1) == -1)
        modelConfig.keepCorrelator = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepCorrelator) < -1) || (max(modelConfig.keepCorrelator) > modelConfig.LRUs))
        error('keepCorrelator field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepCorrelator = unique(modelConfig.keepCorrelator);
    end
end

if ~isfield(modelConfig,'keepPSSBeamformer')
    modelConfig.keepPSSBeamformer = 0;
else
    if (modelConfig.keepPSSBeamformer(1) == -1)
        modelConfig.keepPSSBeamformer = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSSBeamformer) < -1) || (max(modelConfig.keepPSSBeamformer) > modelConfig.LRUs))
        error('keepPSSBeamformer field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSSBeamformer = unique(modelConfig.keepPSSBeamformer);
    end
end

if ~isfield(modelConfig,'keepPSTBeamformer')
    modelConfig.keepPSTBeamformer = 0;
else
    if (modelConfig.keepPSTBeamformer(1) == -1)
        modelConfig.keepPSTBeamformer = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSTBeamformer) < -1) || (max(modelConfig.keepPSTBeamformer) > modelConfig.LRUs))
        error('keepPSTBeamformer field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSTBeamformer = unique(modelConfig.keepPSTBeamformer);
    end
end

if ~isfield(modelConfig,'keepSDPOutput')
    modelConfig.keepSDPOutput = 1:modelConfig.LRUs;
else
    if (modelConfig.keepSDPOutput(1) == -1)
        modelConfig.keepSDPOutput = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepSDPOutput) < -1) || (max(modelConfig.keepSDPOutput) > modelConfig.LRUs))
        error('keepSDPOutput field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepSDPOutput = unique(modelConfig.keepSDPOutput);
    end
end

if ~isfield(modelConfig,'keepPSSOutput')
    modelConfig.keepPSSOutput = 1:modelConfig.LRUs;
else
    if (modelConfig.keepPSSOutput(1) == -1)
        modelConfig.keepPSSOutput = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSSOutput) < -1) || (max(modelConfig.keepPSSOutput) > modelConfig.LRUs))
        error('keepPSSOutput field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSSOutput = unique(modelConfig.keepPSSOutput);
    end
end

if ~isfield(modelConfig,'keepPSTOutput')
    modelConfig.keepPSTOutput = 1:modelConfig.LRUs;
else
    if (modelConfig.keepPSTOutput(1) == -1)
        modelConfig.keepPSTOutput = 1:modelConfig.LRUs;
    elseif ((min(modelConfig.keepPSTOutput) < -1) || (max(modelConfig.keepPSTOutput) > modelConfig.LRUs))
        error('keepPSTOutput field in moduleConfig.txt contains values that are out of range');
    else
        modelConfig.keepPSTOutput = unique(modelConfig.keepPSTOutput);
    end
end


% Put in fullDoppler option
if ~isfield(modelConfig,'fullDoppler')
    modelConfig.fullDoppler = 0;
end

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