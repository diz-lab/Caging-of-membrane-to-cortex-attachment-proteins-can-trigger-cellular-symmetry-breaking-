function [accel] = PIVstrainaccel(x,y,u,v)
% Calculates the instantaneous strain acceleration as the divergence of the
% acceleration field (u,v) at PIV vector locations (x,y).
%APT
%Added smoothing to test 180416
%u=smoothdata(u,3);
%v=smoothdata(v,3);
for i=1:size(u,3)
    accel(:,:,i) = divergence(x(:,:,1),y(:,:,1),u(:,:,i),v(:,:,i));
end

end

