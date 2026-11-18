function runCartPIV(datapath)

    %% Check to see if the stack has been preloaded & load the stack if true
    preloadstack = 0;
    if exist(fullfile(datapath, 'ExpParams.mat'))==2
        load(fullfile(datapath, 'ExpParams.mat'));
        if myIsField(expt, 'path') == 1 && myIsField(expt, 'filename') && ~isempty(expt.path) && ~isempty(expt.filename)
            stack = stackread(fullfile(expt.path, expt.filename), expt.lastim);
            stackdim = size(stack);
            if length(stackdim) == 2, stackdim(3) = 1; end
            preloadstack = 1;
        end
    end
    
    datapath = expt.path;

    clearvars expt

    %% Select image stack to be analyzed
    if preloadstack == 0
        [FileName, PathName] = uigetfile('*.tif');
        stack = stackread(fullfile(PathName,  FileName), expt.lastim);
        stackdim = size(stack);
        if length(stackdim) == 2, stackdim(3) = 1; end
    end

    %% Load Predefined Experimental Parameters
    load(fullfile(datapath, 'ExpParams.mat'));

    %% Run PIV
    [x, y, dxt, dyt, vxt, vyt, expt, pivspec, filtspec] = cartPIV(stack, expt, pivspec, filtspec,datapath);

    %% Calculate PIV statistics
    [dx,dy,vx,vy,vmag,ddiv,vdiv,strain,strainr,straina,curlz,Shear,Shear_diag,Shear_offdiag] = cartPIVstats(x, y, dxt, dyt, vxt, vyt);

    %% Plot the divergence of the velocity over time

    % [vdiv.stdplot] = errorbar_lims(expt.time,vdiv.std);
    % ydn = min((vdiv.mean-vdiv.stdplot));
    % yup = max((vdiv.mean+vdiv.stdplot));
    % [yllim,yulim] = plot_lims(ydn,yup);
    %
    % figure, hold on
    % h = errorbar(expt.time,vdiv.mean,vdiv.std,'ok','MarkerFaceColor','k','MarkerEdgeColor','k');
    % axis square, set(gca,'box','on')
    % xlim([min(expt.time) max(expt.time)])
    % if ~isequaln(yllim,yulim), ylim([yllim yulim]); end
    % xlabel('$time\;(sec)$','interpreter','latex')
    % ylabel(['$\nabla\cdot v\;(s^{-1})$'],'interpreter','latex')
    % set(gca,'fontsize',16)
    %
    %
    % % [savetag] = tcorrsavetag([],[],[],[],[],[],[],1);
    % saveloc({'Plots'},['DivergenceOfVelocityVsTime'],1,h)


    %% Plot PIV overlay
    if myIsField(pivspec, 'pivoverlay')==1 && pivspec.pivoverlay == 1
        savename = 'VelocityVectorFieldOverlay';
        plotPIVoverlay(stack, expt, x, y, vxt, vyt, 0, 0, 1, 1, savename, 0, 1.25, 2, datapath);
    end
end
