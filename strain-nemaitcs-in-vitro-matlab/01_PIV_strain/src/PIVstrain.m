function [strain] = PIVstrain(x,y,u,v)
% Calculates the cumulative strain as the divergence of the displacement
% field (u,v) at PIV vector locations (x,y).

%Check for timepoint with all nans 
%APT added possibilities of multiple nanframes 1902
%APT instead of removing nanframes and readding, just set content of nanframes =0

nf=0;
for i=1:size(u,3)
    if sum(sum(isnan(u(:,:,i))))==size(u,1)*size(u,2)
        if nf==0
            nanframe=[]; % do this for exist() below
            nf=1;
        end
        nanframe=[nanframe,i];
    end
end

% if exist('nanframe','var')
%     u(:,:,nanframe)=[];
%     v(:,:,nanframe)=[];
% end

if exist('nanframe','var')
    for jj=1:length(nanframe)
        u(:,:,nanframe(jj))=zeros(size(u,1),size(u,2));
        v(:,:,nanframe(jj))=zeros(size(u,1),size(u,2));
    end
end

%keyboard

usum = cumsum(u,3);
vsum = cumsum(v,3);

for i=1:size(usum,3)
    strain(:,:,i) = divergence(x(:,:,1),y(:,:,1),usum(:,:,i),vsum(:,:,i) );
end

% %Replace the nan frame in strain, if it existed
% if exist('nanframe','var')
%     %filler=NaN(size(u,1),size(u,2));
%     filler=zeros(size(u,1),size(u,2));
%     for jj=1:length(nanframe)
%     nf=nanframe(jj);
%     
%     if nf==1
%         strain=cat(3,filler,strain);
%     elseif nf==(size(u,3)+1)
%         strain=cat(3,strain,filler);
%     else
%         strain1=cat(3,strain(:,:,(1:(nf-1))),filler);
%         strain=cat(3,strain1,strain(:,:,(nf+1):end));
%         disp('You deleted a slice from the middle of the image');
%     end
%     end
% end
%keyboard

end

