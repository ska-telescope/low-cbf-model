function [skyout] = expandSky(skyin,rundir)
% Adds the field "subIndex", and expands entries in sky that contain an image into a set of objects at each non-zero pixel in the image.
%
%

skyoutCount = 0;
for skyobject = 1:length(skyin)
    skyimgName = skyin(skyobject).skyImage{1};
    if (~isempty(skyimgName) > 0)
        % Filename provided for an image of the sky, expand it into a set of objects for each non-zero pixel
        % Each object created has a different subindex field. 
        % subindex = 0 is not an object, it is just used for pointing to the center of the field.
        
        % Create subindex 0
        skyoutCount = skyoutCount + 1;
        skyout(skyoutCount) = skyin(skyobject);
        subindex = 0;
        skyout_subindex(skyoutCount) = subindex;  % subindex goes into the struct at the end, otherwise we cannot assign from skyin to skyout as the structures become different.
        
        % Create entries for each non-zero point in the image
        skyimg = imread([rundir '/' skyimgName]);
        % Convert to greyscale, with range 0 to 1
        if (length(size(skyimg)) == 3)
            skyimg = sum(skyimg,3); % just add up all the colours.
        end
        skyimg = skyimg/max(max(skyimg)); % scale to range 0 to 1
        skysize = size(skyimg);
        for r = 1:skysize(1)
            for c = 1:skysize(2)
                if (skyimg(r,c) ~= 0)
                    skyoutCount = skyoutCount + 1;
                    subindex = subindex + 1;
                    skyout(skyoutCount) = skyin(skyobject);
                    skyout_power(skyoutCount) = 1;
                    skyout(skyoutCount).ascension = skyout(skyoutCount).ascension - skyout(skyoutCount).skyImageDimension * (((r-1) - (skysize(1)-1)/2)/(skysize(1)-1));
                    skyout(skyoutCount).declination = skyout(skyoutCount).declination - skyout(skyoutCount).skyImageDimension * (((c-1) - (skysize(2)-1)/2)/(skysize(2)-1));
                    skyout_power(skyoutCount) = skyout_power(skyoutCount) * skyimg(r,c);
                    skyout_subindex(skyoutCount) = subindex;
                end
            end
        end
    else
        % no filename, just copy the entry in skyin, and assign subindex
        skyoutCount = skyoutCount + 1;
        skyout(skyoutCount) = skyin(skyobject);
        skyout_power(skyoutCount) = 1;  % Only used for images.
        skyout_subindex(skyoutCount) = 0;
    end
end

for skyobject = 1:skyoutCount
    skyout(skyobject).subIndex = skyout_subindex(skyobject);
    skyout(skyobject).power = skyout_power(skyobject);
end


