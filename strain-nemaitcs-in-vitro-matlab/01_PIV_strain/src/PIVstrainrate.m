function [sr] = PIVstrainrate(x,y,u,v)
% Calculates the instantaneous strain rate as the divergence of the
% velocity field (u,v) at PIV vector locations (x,y).

for i=1:size(u,3)
    sr(:,:,i) = divergence(x(:,:,1),y(:,:,1),u(:,:,i),v(:,:,i));
end
end

