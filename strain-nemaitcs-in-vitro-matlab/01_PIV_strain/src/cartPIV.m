function [x, y, dxt, dyt, vxt, vyt, expt, pivspec, filtspec] = cartPIV(stack, expt, pivspec, filtspec, datapath)

    if exist('stack')~=1 || isempty(stack)
        [filename,pathname]=uigetfile({'*.tif';'*.png';'*.jpg';'*.pdf';'*.gif'});
        stack = stackread(fullfile(pathname,filename));
    end

    mkexp =0;
    if exist('expt')~=1 || isempty(expt), mkexp = 1; end
    if mkexp == 1 || myIsField(expt,'dscl')~=1, expt.dscl = 1; end
    if mkexp == 1 || myIsField(expt,'tscl')~=1, expt.tscl = 1; end

    mkpivspec = 0;
    if exist('pivspec')~=1 || isempty(pivspec), mkpivspec = 1; end
    if mkpivspec == 1 || myIsField(pivspec,'imdist')~=1, pivspec.imdist = 1; end
    if mkpivspec == 1 || myIsField(pivspec,'firstim')~=1, pivspec.firstim = 1; end
    if mkpivspec == 1 || myIsField(pivspec,'lastim')~=1, pivspec.lastim = 0; end
    if mkpivspec == 1 || myIsField(pivspec,'winsize')~=1, pivspec.winsize=64; end
    if mkpivspec == 1 || myIsField(pivspec,'overlap')~=1, pivspec.overlap=0.5; end
    if mkpivspec == 1 || myIsField(pivspec,'method')~=1, pivspec.method='single'; end

    mkfiltspec = 0;
    if exist('filtspec')~=1 || isempty(filtspec), mkfiltspec = 1; end
    if mkfiltspec == 1 || myIsField(filtspec,'snr')~=1, filtspec.snr = 1.001; end
    if mkfiltspec == 1 || myIsField(filtspec,'globfiltstd')~=1, filtspec.globfiltstd = 3; end %standard deviation used as threshold for globalfilt function
    if mkfiltspec == 1 || myIsField(filtspec,'locfiltstd')~=1, filtspec.locfiltstd = 2; end %standard deviation used as threshold for localfilt function
    if mkfiltspec == 1 || myIsField(filtspec,'kernal')~=1, filtspec.kernal = 3; end %kernal size used in localfilt function
    if mkfiltspec == 1 || myIsField(filtspec,'locfiltmeth')~=1, filtspec.locfiltmeth = 'median'; end %method used in localfilt function

    stackdim=size(stack);
    if pivspec.lastim == 0
       pivspec.lastim = stackdim(3);
    end

    for i=1:(pivspec.lastim-pivspec.imdist)
        ii = i + pivspec.imdist;
        disp(['Analyzing Image Pair: ',num2str(i)]);
        %Read in the image and display it
        img1=stack(:,:,i);
        img2=stack(:,:,ii);

        %Calculate vector field between the two chosen images
        [xpix,ypix,upivpix,vpivpix,snrout]=matpiv(img1,img2,pivspec.winsize,pivspec.imdist,pivspec.overlap,pivspec.method);

        if myIsField(filtspec,'filt') == 0 || filtspec.filt == 1
            [upix,vpix] = PIVfilt(xpix,ypix,upivpix,vpivpix,snrout,filtspec);
        else
            upix = upivpix;
            vpix = vpivpix;
        end

        x = xpix * expt.dscl;
        y = ypix * expt.dscl;

        if i < pivspec.firstim
            dxt = NaN([size(xpix), pivspec.lastim]);
            dyt = NaN([size(xpix), pivspec.lastim]);
            vxt = NaN([size(xpix), pivspec.lastim]);
            vyt = NaN([size(xpix), pivspec.lastim]);
        end

        dxt(:, :, ii) = upix * (expt.dscl);
        dyt(:, :, ii) = vpix * (expt.dscl);
        vxt(:, :, ii) = upix * (expt.dscl ./ expt.tscl);
        vyt(:, :, ii) = vpix * (expt.dscl ./ expt.tscl);
    end

    clear unan vnan su sv gu gv lu lv upivpix vpivpix upix vpix img1 img2 i

    % save([pwd,'\ExpParams.mat'],'exp','pivspec','filtspec');
    % clear exp

    if exist(fullfile(datapath, 'WS'), 'dir') == 0, mkdir(fullfile(datapath, 'WS')); end
    save(fullfile(datapath, 'WS', 'cartPIV_WS.mat'), '-regexp', '^(?!(stack)$).');

end
