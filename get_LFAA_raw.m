function get_LFAA_raw(rundir)
% Loads LFAA.mat data file and converts to the raw format used to send data to the FPGAs.
% In the raw format, there are two data files per FPGA.
%  - "LFAAData_fpgaX.raw"
%    This is a direct transcription of fpga(X).data
%  - "LFAAHeaders_fpgaX.raw"
%    List of headers to send, each header is 148 bytes, containing
%    1:8    : Byte offset into LFAAData_fpgaX.raw for the data part of the packet
%    9:12   : Number of bytes to send in the data part of the packet (usually 8192)
%    13:20  : Time at which to send the packet
%    21:32  : Unused (0)
%    33:146 : SPEAD header
%    147:148: Unused (0).
%

%% Load the data file

fname = [rundir '/LFAA'];

if exist(fullfile(cd,rundir,'LFAA.mat'),'file')
    load(fname)
else
    error(['Cannot find ' fname '. Run create_config.m to generate it']);
end


% Convert to raw format.
for f1 = 1:length(fpga)
   %% First the data file, LFAAData_fpgaX.raw
   fname = [rundir '/LFAAData_fpga' num2str(f1) '.raw'];
   fid = fopen(fname,'w');
   s1 = size(fpga(f1).data);
   bytesPerChannel = 2*s1(1);
   thisChannel = zeros(bytesPerChannel,1);
   for b1 = 1:s1(2)
       thisChannel(1:2:end) = real(fpga(f1).data(:,b1));  % Note fpga(XX).data is complex int8 data.
       thisChannel(2:2:end) = imag(fpga(f1).data(:,b1));
       fwrite(fid,thisChannel,'int8');
   end
   fclose(fid);
   
   % read it back to check
   %fid = fopen(fname,'r');
   %check = fread(fid);
   %fclose(fid);
   
   %% The header data file
   fname = [rundir '/LFAAHeaders_fpga' num2str(f1) '.raw'];
   fid = fopen(fname,'w');
   s1 = size(fpga(f1).headers);
   hdr = zeros(148,1,'uint8');
   for h1 = 1:s1(2)
       % ptr(1) is the offset within the channel of the start of the data. Note this is samples (with 1 based indexing), so we need 2*(ptr(1)-1) for the offset in bytes.
       % ptr(2) is the channel. If ptr(2) = 0, the data is zeros, and we mark the byte_offset field in LFAAHeaders_fpgaX.raw as 0xffffffffffffffff
       ptr = fpga(f1).dataPointers(h1,:);
       if (ptr(2) == 0)
           hdr(1:8) = 255;
       else
           fullOffset = (ptr(2) - 1) * bytesPerChannel + (ptr(1)-1)*2;
           hdr(1:2) = 0;  % Cannot need more than 6 bytes (=262144 Gbytes)
           hdr(3) = uint8(floor(fullOffset/2^40));
           hdr(4) = uint8(mod(floor(fullOffset/2^32),256));
           hdr(5) = uint8(mod(floor(fullOffset/2^24),256));
           hdr(6) = uint8(mod(floor(fullOffset/2^16),256));
           hdr(7) = uint8(mod(floor(fullOffset/2^8),256));
           hdr(8) = uint8(mod(fullOffset,256));
       end
       % Number of bytes in the data part of the packet = 8192 = 0x2000
       hdr(9) = 0;
       hdr(10) = 0;
       hdr(11) = 32;
       hdr(12) = 0;
       % Time to send the packet in nanoseconds
       hdr(13) = mod(floor(fpga(f1).txTimes(h1)/2^56),256);
       hdr(14) = mod(floor(fpga(f1).txTimes(h1)/2^48),256);
       hdr(15) = mod(floor(fpga(f1).txTimes(h1)/2^40),256);
       hdr(16) = mod(floor(fpga(f1).txTimes(h1)/2^32),256);
       hdr(17) = mod(floor(fpga(f1).txTimes(h1)/2^24),256);
       hdr(18) = mod(floor(fpga(f1).txTimes(h1)/2^16),256);
       hdr(19) = mod(floor(fpga(f1).txTimes(h1)/2^8),256);
       hdr(20) = mod(fpga(f1).txTimes(h1),256);
       % unused
       hdr(21:32) = 0;
       hdr(33:146) = fpga(f1).headers(:,h1);
       hdr(147:148) = 0;
       fwrite(fid,hdr,'uint8');
   end
   fclose(fid);
   %keyboard
end


