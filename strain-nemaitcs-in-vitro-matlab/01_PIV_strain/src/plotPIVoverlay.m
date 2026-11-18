function plotPIVoverlay(stack, expt, x, y, u, v, range, scale, plottype, bsave, savename, bdisplay, lw, mhs, datapath)
% Plots PIV vector field overlay onto an image stack
% plottype - specifies vector components to overlay onto the images
% 1: Resultant vectors
% 2: Radial vector components
% 3: Tangential vector components
% 4: Radial & Tangential vector components
% 5: Radial & Resultant vectors
% 6: Tangential & Resultant vectors
%
% Input Parameters
% expt - structure containing experimental parameters
% x - x-coordinates (um)
% y - y-coordinates (um)
% u - x-component of vector field (um/s)
% v - y-component of vector field (um/s)
% range - range (pixels) over which to plot the vector field. range = 0
%         will plot all vectors in the field
% scale - vector scaling factor. scale = 0 use the default vector scaling
% plottype - see above
% bsave - save boolean specifying whether to save the vector overlays
%       bsave = 0 - don't save the PIV overlays
%       bsave = 1 - just save tiffs
%       bsave = 2 - save tiffs and figs
%       bsave = 3 - save tiffs, figs, and epsc
% bdisplay - display boolean specifiying whether to show the images as the
%            functon runs
% lw - Line Width of the overaid vectors. Default = 1.75
% mhs - Max head size of overlaid vectors. Default = 2

if bsave == 1
    if exist(fullfile(datapath,'PIVOverlay'), 'dir') == 0; mkdir(fullfile(datapath,'PIVOverlay'));end
end

if exist('savename') == 0 || isempty(savename)
    switch plottype
        case 1
            savename = 'ResultantVectorOverlay';
%             if exist(fullfile(datapath,'PIVOverlay', 'ResultantVectorOverlay'), 'dir') == 0, mkdir(fullfile(datapath,'PIVOverlay', 'ResultantVectorOverlay'));end
        case 2
            savename = 'RadialVectorOverlay';
%             if exist(fullfile(datapath,'PIVOverlay', 'RadialVectorOverlay'), 'dir') == 0,mkdir(fullfile(datapath,'PIVOverlay', 'RadialVectorOverlay'));end
        case 3
            savename = 'TangentialVectorOverlay';
%             if exist(fullfile(datapath,'PIVOverlay', 'TangentialVectorOverlay'), 'dir') == 0,mkdir(fullfile(datapath,'PIVOverlay', 'TangentialVectorOverlay'));end
        case 4
            savename = 'Radial&TangentialVectorOverlay';
%             if exist(fullfile(datapath,'PIVOverlay', 'Radial&TangentialVectorOverlay'), 'dir') == 0,mkdir(fullfile(datapath,'PIVOverlay', 'Radial&TangentialVectorOverlay'));end
        case 5
            savename = 'Radial&ResultantVectorOverlay';
%             if exist(fullfile(datapath,'PIVOverlay', 'Radial&ResultantVectorOverlay'), 'dir') == 0,mkdir(fullfile(datapath,'PIVOverlay', 'Radial&ResultantVectorOverlay'));end
        case 6
            savename = 'Tangential&ResultantVectorOverlay';
%             if exist(fullfile(datapath,'PIVOverlay', 'Tangential&ResultantVectorOverlay'), 'dir') == 0,mkdir(fullfile(datapath,'PIVOverlay', 'Tangential&ResultantVectorOverlay'));end
        otherwise
           plottype = 1; bsave = 0; bdisplay = 1;
    end
    if exist(fullfile(datapath,'PIVOverlay',savename), 'dir') == 0,mkdir(fullfile(datapath,'PIVOverlay',savename));end
else
    if exist(fullfile(datapath,'PIVOverlay',savename), 'dir') == 0,mkdir(fullfile(datapath,'PIVOverlay',savename));end
end

if exist('lw')~=1, lw = 1.75; end
if exist('mhs')~=1, mhs=2; end

xpix = x./expt.dscl;
ypix = y./expt.dscl;
[xpol,ypol,th,r] = cart2polarcoord(expt.cp,xpix,ypix);
xspace = xpix(1,2)-xpix(1,1);
yspace = ypix(2,1)-ypix(1,1);

if range == 0, range = [0,inf];end
bbil = find(r(:)<(range(1)*expt.dscl));
bbig = find(r(:)>(range(2)*expt.dscl));
bbi = [bbil;bbig];

for i=1:size(u,3)
    vxtemp = u(:,:,i);
    vxtemp(bbi) = nan;
    vxrange(:,:,i) = vxtemp;
    vytemp = v(:,:,i);
    vytemp(bbi) = nan;
    vyrange(:,:,i) = vytemp;
end

[vr,vt,~,~] = cart2polvectfield(xpol,ypol,vxrange,vyrange);

Uscale=xspace/max(abs(vxrange(:)));
Vscale=yspace/max(abs(vyrange(:)));
defaultscale = max([Uscale,Vscale]);
disp(['Default Scale = ',num2str(defaultscale)]);
if scale == 0, scale = defaultscale; end

% stack8 = im2uint8(stack);
for f=1:size(stack,3)
    sadj(:,:,f) = imadjust(stack(:,:,f));
end
stack8 = im2uint8(sadj);

for i=1:size(u,3)
    tic
%     img=stack(:,:,i);
    img=stack8(:,:,i);

    vxplot = vxrange(:,:,i);
    vyplot = vyrange(:,:,i);

    vr_xplot = cos(th).*vr(:,:,i);
    vr_yplot = sin(th).*vr(:,:,i);

    vt_xplot = -1*sin(th).*vt(:,:,i);
    vt_yplot = cos(th).*vt(:,:,i);

    %% Plot Vector Overlays
    switch plottype
        case 1
            if bdisplay == 0, h = figure('visible','off'); end
            imagesc(img);
            colormap(gray);
            hold on

            h1=quiver(xpix,ypix,vxplot,vyplot,0,'r','LineWidth',lw,'MaxHeadSize',mhs);
            hU=get(h1,'UData');
            hV=get(h1,'VData');
            set(h1,'UData',scale*hU,'VData',scale*hV);
            set(gca,'LineWidth',2);
            axis equal
            axis tight
            axis off;
            width = size(img,2)./100; % inches
            height = size(img,1)./100; % inches
            tightPlot(h,width,height);

            if bsave >= 1
                saveas(gca, fullfile(datapath,'PIVOverlay',savename, ['TimePoint_',num2str(i),'.tif'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 2
                saveas(gca, fullfile(datapath,'PIVOverlay',savename, ['TimePoint_',num2str(i),'.fig'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 3
                saveas(gca, fullfile(datapath,'PIVOverlay',savename, ['TimePoint_',num2str(i),'.epsc'])); %Save the image overlayed with the quiver plot
            end
            close all
        case 2
            if bdisplay == 0, figure('visible','off'); end
            %Plot quiver output as image overlays for radial component of velocity vector
            imagesc(img); %loads the second image of the image pair, scales the image data to the full range of the current colormap and displays the image
            colormap(gray);%sets the colormap to gray
            hold on %retains the displayed image

            h_r=quiver(xpix,ypix,vr_xplot,vr_yplot,0,'c','LineWidth',lw,'MaxHeadSize',mhs);%Generates radial velocity vector field overlay plots and assigns it to a handle
            hU_r=get(h_r,'UData');%retrieves x-component of the radial velocity vector field plot data from quiver handle
            hV_r=get(h_r,'VData');%retrieves y-component of the radial velocity vector field plot data from quiver handle
            set(h_r,'UData',scale*hU_r,'VData',scale*hV_r);%plots the two velocity vector components as single vector field and scales them accordingly
            set(gca,'LineWidth',2);%plots velocity vector field overlay on current axes
            axis equal;%squares the axes
            axis off;%removes the axis labels

            if bsave >= 1
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 2
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 3
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.epsc'])); %Save the image overlayed with the quiver plot
            end
%             if bsave == 1
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif']));%Save the image overlayed with the radial velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig']));%Save the image overlayed with the radial velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.eps']));%Save the image overlayed with the radial velocity quiver plot
%             end
            close all
        case 3
            if bdisplay == 0, figure('visible','off'); end
            %Plot quiver output as image overlays for tangential component of velocity vector
            imagesc(img); %loads the second image of the image pair, scales the image data to the full range of the current colormap and displays the image
            colormap(gray);%sets the colormap to gray
            hold on %retains the displayed image

            h_th=quiver(xpix,ypix,vt_xplot,vt_yplot,0,'m','LineWidth',lw,'MaxHeadSize',mhs);%Generates tangential velocity vector field overlay plots and assigns it to a handle
            hU_th=get(h_th,'UData');%retrieves x-component of the tangential velocity vector field plot data from quiver handle
            hV_th=get(h_th,'VData');%retrieves y-component of the tangential velocity vector field plot data from quiver handle
            set(h_th,'UData',scale*hU_th,'VData',scale*hV_th);%plots the two velocity vector components as single vector field and scales them accordingly
            set(gca,'LineWidth',2);%plots velocity vector field overlay on current axes
            axis equal;%squares the axes
            axis off;%removes the axis labels

            if bsave >= 1
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 2
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 3
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.epsc'])); %Save the image overlayed with the quiver plot
            end
%             if bsave == 1
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif']));%Save the image overlayed with the tangential velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig']));%Save the image overlayed with the tangential velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.eps']));%Save the image overlayed with the tangential velocity quiver plot
%             end
            close all
        case 4
            if bdisplay == 0, figure('visible','off'); end
            %Plot quiver output as image overlays for tangential component of velocity vector
            imagesc(img); %loads the second image of the image pair, scales the image data to the full range of the current colormap and displays the image
            colormap(gray);%sets the colormap to gray
            hold on %retains the displayed image

            h_r=quiver(xpix,ypix,vr_xplot,vr_yplot,0,'c','LineWidth',lw,'MaxHeadSize',mhs);%Generates radial velocity vector field overlay plots and assigns it to a handle
            hU_r=get(h_r,'UData');%retrieves x-component of the radial velocity vector field plot data from quiver handle
            hV_r=get(h_r,'VData');%retrieves y-component of the radial velocity vector field plot data from quiver handle
            set(h_r,'UData',scale*hU_r,'VData',scale*hV_r);%plots the two velocity vector components as single vector field and scales them accordingly
            set(gca,'LineWidth',2);%plots velocity vector field overlay on current axes
            hold on
            h_th=quiver(xpix,ypix,vt_xplot,vt_yplot,0,'m','LineWidth',lw,'MaxHeadSize',mhs);%Generates tangential velocity vector field overlay plots and assigns it to a handle
            hU_th=get(h_th,'UData');%retrieves x-component of the tangential velocity vector field plot data from quiver handle
            hV_th=get(h_th,'VData');%retrieves y-component of the tangential velocity vector field plot data from quiver handle
            set(h_th,'UData',scale*hU_th,'VData',scale*hV_th);%plots the two velocity vector components as single vector field and scales them accordingly
            set(gca,'LineWidth',2);%plots velocity vector field overlay on current axes
            axis equal;%squares the axes
            axis off;%removes the axis labels

            if bsave >= 1
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 2
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 3
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.epsc'])); %Save the image overlayed with the quiver plot
            end
%             if bsave == 1
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif']));%Save the image overlayed with the tangential velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig']));%Save the image overlayed with the tangential velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.eps']));%Save the image overlayed with the tangential velocity quiver plot
%             end
            close all
        case 5
            if bdisplay == 0, figure('visible','off'); end
            %Plot quiver output as image overlays for velocity vector
            imagesc(img); %loads the second image of the image pair, scales the image data to the full range of the current colormap and displays the image
            colormap(gray);%sets the colormap to gray
            hold on %retains the displayed image

            h1=quiver(xpix,ypix,vxplot,vyplot,0,'r','LineWidth',lw,'MaxHeadSize',mhs);
            hU=get(h1,'UData');
            hV=get(h1,'VData');
            set(h1,'UData',scale*hU,'VData',scale*hV);
            set(gca,'LineWidth',2);
            hold on
            h_r=quiver(xpix,ypix,vr_xplot,vr_yplot,0,'c','LineWidth',lw,'MaxHeadSize',mhs);%Generates radial velocity vector field overlay plots and assigns it to a handle
            hU_r=get(h_r,'UData');%retrieves x-component of the radial velocity vector field plot data from quiver handle
            hV_r=get(h_r,'VData');%retrieves y-component of the radial velocity vector field plot data from quiver handle
            set(h_r,'UData',scale*hU_r,'VData',scale*hV_r);%plots the two velocity vector components as single vector field and scales them accordingly
            set(gca,'LineWidth',2);%plots velocity vector field overlay on current axes
            axis equal;%squares the axes
            axis off;%removes the axis labels

            if bsave >= 1
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 2
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 3
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.epsc'])); %Save the image overlayed with the quiver plot
            end
%             if bsave == 1
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif']));%Save the image overlayed with the radial velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig']));%Save the image overlayed with the radial velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.eps']));%Save the image overlayed with the radial velocity quiver plot
%             end
            close all
        case 6
            if bdisplay == 0, figure('visible','off'); end
            %Plot quiver output as image overlays for velocity vector
            imagesc(img); %loads the second image of the image pair, scales the image data to the full range of the current colormap and displays the image
            colormap(gray);%sets the colormap to gray
            hold on %retains the displayed image

            h1=quiver(xpix,ypix,vxplot,vyplot,0,'r','LineWidth',lw,'MaxHeadSize',mhs);
            hU=get(h1,'UData');
            hV=get(h1,'VData');
            set(h1,'UData',scale*hU,'VData',scale*hV);
            set(gca,'LineWidth',2);
            hold on
            h_r=quiver(xpix,ypix,vr_xplot,vr_yplot,0,'c','LineWidth',lw,'MaxHeadSize',mhs);%Generates radial velocity vector field overlay plots and assigns it to a handle
            hU_r=get(h_r,'UData');%retrieves x-component of the radial velocity vector field plot data from quiver handle
            hV_r=get(h_r,'VData');%retrieves y-component of the radial velocity vector field plot data from quiver handle
            set(h_r,'UData',scale*hU_r,'VData',scale*hV_r);%plots the two velocity vector components as single vector field and scales them accordingly
            set(gca,'LineWidth',2);%plots velocity vector field overlay on current axes
            axis equal;%squares the axes
            axis off;%removes the axis labels

            if bsave >= 1
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 2
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig'])); %Save the image overlayed with the quiver plot
            end
            if bsave >= 3
                saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.epsc'])); %Save the image overlayed with the quiver plot
            end
%             if bsave == 1
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.tif']));%Save the image overlayed with the radial velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.fig']));%Save the image overlayed with the radial velocity quiver plot
%                 saveas(gca, fullfile(datapath,'PIVOverlay',savename,['TimePoint_',num2str(i),'.eps']));%Save the image overlayed with the radial velocity quiver plot
%             end
            close all
    end
    disp(['Finished plottinng PIV overlay on frame ', num2str(i), ' out of ', num2str(size(u,3))])
    toc
end
end


