function [im,info] = stackread(fullpath, maxFrame)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This program reads in a multiplaned TIFF file
    % into a single variable
    %
    % 11/12/09 Patrick Oakes
    % poakes@gmail.com
    %
    % Updated
    % 5/26/16 - Ian Linsmeier
    % - Changed function name
    % - Compatible with 8,16, & 32-bit images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin == 1
        maxFrame = 0;
    end

    if ~exist('fullpath') || isempty(fullpath)
        [file,path] = uigetfile({'*.tif';'*.png';'*.jpg';'*.pdf';'*.gif'});
        if isequal(path, 0) || isequal(file, 0), return; end
        fullpath = fullfile(path,file);
    end

    info = imfinfo(fullpath);
    dim = size(info);
    if dim(1) > dim(2)
        num_planes = dim(1);
    else
        num_planes = dim(2);
    end
    im = zeros(info(1).Height,info(1).Width,num_planes);

    if double(info(1).BitDepth) == 24
        im = uint8(zeros(info(1).Height,info(1).Width,3,num_planes));
    end

    if maxFrame ~= 0
        num_planes = maxFrame;
    end

    if num_planes == 1
        im = imread(fullpath);
    else
        for ii = 1:num_planes
            disp(['frame ', num2str(ii)])
            if double(info(1).BitDepth) == 24
                im(:,:,:,ii) = imread(fullpath,ii);
            else
                im(:,:,ii) = imread(fullpath,ii);
            end
        end
    end

    if double(info(1).BitDepth) == 8
        im=uint8(im);
    elseif double(info(1).BitDepth) == 16
        im=uint16(im);
    elseif double(info(1).BitDepth) == 32
        im=uint32(im);
    end
end
