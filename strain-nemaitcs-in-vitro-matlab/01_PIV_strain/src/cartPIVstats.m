function [dx,dy,vx,vy,vmag,ddiv,vdiv,strain,strainr,straina,curlz,Shear,Shear_diag,Shear_offdiag] = cartPIVstats(x,y,dxt,dyt,vxt,vyt)
dim = size(dxt,3);

dx.t = dxt;
dx.mean = nanmean(reshape(dxt,[size(dxt,1)*size(dxt,2), size(dxt,3)]));
dx.std = nanstd(reshape(dxt,[size(dxt,1)*size(dxt,2), size(dxt,3)]));

dy.t = dyt;
dy.mean = nanmean(reshape(dyt,[size(dyt,1)*size(dyt,2), size(dyt,3)]));
dy.std = nanstd(reshape(dyt,[size(dyt,1)*size(dyt,2), size(dyt,3)]));

vx.t = vxt;
vx.mean = nanmean(reshape(vxt,[size(vxt,1)*size(vxt,2), size(vxt,3)]));
vx.std = nanstd(reshape(vxt,[size(vxt,1)*size(vxt,2), size(vxt,3)]));

vy.t = vyt;
vy.mean = nanmean(reshape(vyt,[size(vyt,1)*size(vyt,2), size(vyt,3)]));
vy.std = nanstd(reshape(vyt,[size(vyt,1)*size(vyt,2), size(vyt,3)]));
   
vmagt = sqrt(vxt.^2+vyt.^2);
vmag.t = vmagt;
vmag.mean = nanmean(reshape(vmagt,[size(vmagt,1)*size(vmagt,2), size(vmagt,3)]));
vmag.std = nanstd(reshape(vmagt,[size(vmagt,1)*size(vmagt,2), size(vmagt,3)]));

ddivt = NaN(size(dxt));
vdivt = NaN(size(vxt));

for i=1:dim
    ddivt(:,:,i) = divergence(x,y,dxt(:,:,i),dyt(:,:,i));
    vdivt(:,:,i) = divergence(x,y,vxt(:,:,i),vyt(:,:,i));
end

ddiv.t = ddivt;
ddiv.mean = nanmean(reshape(ddivt,[size(ddivt,1)*size(ddivt,2), size(ddivt,3)]));
ddiv.std = nanstd(reshape(ddivt,[size(ddivt,1)*size(ddivt,2), size(ddivt,3)]));

vdiv.t = vdivt;
vdiv.mean = nanmean(reshape(vdivt,[size(vdivt,1)*size(vdivt,2), size(vdivt,3)]));
vdiv.std = nanstd(reshape(vdivt,[size(vdivt,1)*size(vdivt,2), size(vdivt,3)]));

%Calculate strain, strain rate, and acceleration
strain= PIVstrain(x,y,dx.t,dy.t);
strainrate= PIVstrainrate(x,y,vx.t,vy.t);
strainaccel = PIVstrainaccel(x,y,vxt,vyt);

straina.t=strainaccel;
straina.mean = nanmean(reshape(strainaccel,[size(strainaccel,1)*size(strainaccel,2), size(strainaccel,3)]));
straina.std= nanstd(reshape(strainaccel,[size(strainaccel,1)*size(strainaccel,2), size(strainaccel,3)]));

strainr.t=strainrate;
strainr.mean = nanmean(reshape(strainrate,[size(strainrate,1)*size(strainrate,2), size(strainrate,3)]));
strainr.std = nanstd(reshape(strainrate,[size(strainrate,1)*size(strainrate,2), size(strainrate,3)]));

s=strain;
clear strain
strain.t=s;
strain.mean = nanmean(reshape(s,[size(s,1)*size(s,2), size(s,3)]));
strain.std = nanstd(reshape(s,[size(s,1)*size(s,2), size(s,3)]));

%Calculate curl of velocity vectors
for i=1:size(vxt,3)
    [kurl,~]=curl(x,y,vxt(:,:,i),vyt(:,:,i));
    curlz(:,:,i)=kurl;
end

%Calculate shear of velocity vectors

for j = 1: size(vxt,3)
    [Sheartot,Ushear1, Ushear2] = PIVShear(x,y,vxt(:,:,j),vyt(:,:,j));
    Shear(:,:,j) = Sheartot;
    Shear_diag(:,:,j) = Ushear1;
    Shear_offdiag(:,:,j) = Ushear2;
end

clearvars -except dx dy vx vy vmag ddiv vdiv strain strainr straina curlz Shear Shear_diag Shear_offdiag

if exist([pwd,'/WS'], 'dir') == 0, mkdir([pwd,'/WS']); end
save([pwd,'/WS/cartPIVstats_WS.mat']);
end     

