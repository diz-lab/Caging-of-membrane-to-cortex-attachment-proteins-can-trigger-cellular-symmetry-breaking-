% function alignment(dataDir, imgName, varargin)
%
% Wrapper for FFTAlignment. Calculates the alignment as for each frame of a tiff stack
% named imgName in the full path dataDir and saves the output to dataDir/alignment.
%
% Options are passed to FFTAlignment by name/value pairs, e.g. call
% alignment(dataDir, imgName, 'checkpoint', 100, 'winsize', 33)
%
%
% INPUT
% =====
% dataDir     = full directory to data to be analyzed.
% imgName     = name of image inside of dataDir to analyze.
% winsize     = Size of the sub windows. Must be odd. (Optional, default is 33)
% overlap     = How much the windows overlap. Value should be between 0
%               (complete overlap) and 1 (no overlap). (Optional, default is
%               0.5)
% st          =   order parameter spacing. Considers a window of size
%               2*st+1 to compare vectors to. Depending on the
%               overlap size, you'll want to change this. (Optional, default is
%               2)
% checkpoint  = threshold value above which to do the caculation. The
%               sum of the image intensity in the window needs to be
%               above this value for the routine to calculate a vector
%               for the given window (optional, default is 0 - i.e. every
%               window)
% mask_method = Determines the method to mask the FFT when
%               determining the moment calculation. Global (val = 1) uses
%               a circle for every window. Local (val = 2) uses a
%               local threshold for each subwindow. (Optional, default = 1)
% figures     = Plot figures. 0 = No; 1 = Yes; (Optional, default = 0)
%
% OUTPUT
% ======
% saves all data to dataDir/alignment
%
% Written by: Shiladitya Banerjee
function alignment(dataDir, imgName, varargin)

    % Default values and acceptable options for varargin
    options = struct('winsize', 33, 'overlap', 0.5, 'st', 2, 'checkpoint', 0, 'mask_method', 1, 'figures', 0);
    optionNames = fieldnames(options);

    if nargin < 2
        error('alignment needs at least 2 inputs. See documentation with `help alignment`')
    end

    nArgs = length(varargin);
    if round(nArgs/2) ~= nArgs/2
        error('ALIGNMENT needs name/value pairs')
    end

    for pair = reshape(varargin,2,[])
        inputName = lower(pair{1}); % make case insensitive
        if any(strcmp(inputName, optionNames))
            % overwrite options
            options.(inputName) = pair{2};
        else
            error('%s is not a recognized parameter name', inputName)
        end
    end

    tmax = numel(imfinfo(fullfile(dataDir, imgName)));
    fname = fullfile(dataDir, imgName);
    outdir = fullfile(dataDir, 'alignment');

    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end

    for t=1:tmax
        filename1 = ['order_' num2str(t) '.mat'];
        filename2 = ['align_' num2str(t) '.mat'];

        matrixfile1 = fullfile(outdir, filename1);
        matrixfile2 = fullfile(outdir, filename2);

        im_actin = imread(fname,t);
        Pt = FFTAlignment(im_actin, options.winsize, options.overlap, options.st,...
                          options.checkpoint, options.mask_method, options.figures);

        pos_list = Pt.pos;
        vec_data = Pt.vec;
        pos_vec = [pos_list, vec_data];
        order = Pt.ordermat;

        save(matrixfile1, 'order');
        save(matrixfile2, 'pos_vec');
        sprintf([num2str(t), '\b']);
    end
        
    %added easy visualization APT 1802
    %figure(3), imshow(mat2gray(im_actin))
    figure(3), imshow(im_actin,[0, 1]); 
    set(gca,'YDir','normal')
   %figure(3), pcolor(mat2gray(im_actin))
    hold on
    p=pos_vec;
    quiver(p(:,1),p(:,2),p(:,3),p(:,4),'ShowArrowHead','off','color',[1 0 0]) %This is the correct orientation
    %quiver(p(:,1),p(:,2),p(:,3),p(:,4),'ShowArrowHead','off','color',[1 0 0]) %This is the correct orientation
    %set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse','Visible','off','LineWidth',1.5);
    hold off
    saveas(gcf,fullfile(outdir,'nematicAlignment.pdf'),'pdf')
    %saveas(gcf,fullfile(outdir,'nematicAlignment.png'),'png')
    %saveas(gcf,fullfile(outdir,'nematicAlignment.fig'),'fig')

    save(fullfile(dataDir, 'alignmentParams.mat'), 'options')  
    
    
end
