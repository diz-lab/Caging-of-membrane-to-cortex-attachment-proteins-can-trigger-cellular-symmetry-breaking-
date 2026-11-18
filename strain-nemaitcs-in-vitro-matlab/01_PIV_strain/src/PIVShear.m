function [Usheartotal,Ushear1, Ushear2] = PIVShear(x,y,u,v)
%written by Zachary Gao Sun July 29th, 2020.
tau1 = [1 0;0 -1];
tau2 = [0 1;1 0];
Ushear1 = zeros(size(u));
Ushear2 = zeros(size(Ushear1));
Usheartotal = zeros(size(Ushear2));
    for ii = 1: size(u,3)
        [Ux(:,:,ii),Uy(:,:,ii)] = gradient(u(:,:,ii),x(1,2)-x(1,1),y(2,1)-y(1,1));
        [Vx(:,:,ii),Vy(:,:,ii)] = gradient(v(:,:,ii),x(1,2)-x(1,1),y(2,1)-y(1,1));
        

    for jj = 1:size(u,1)
        for kk = 1:size(u,2)
            Ushear1(jj,kk,ii) = tau1(1,1)*Ux(jj,kk,ii)+tau1(1,2)*Uy(jj,kk,ii)+tau1(2,1)*Vx(jj,kk,ii)+tau1(2,2)*Vy(jj,kk,ii);
            Ushear1(jj,kk,ii) = abs(Ushear1(jj,kk,ii));
            Ushear2(jj,kk,ii) = tau2(1,1)*Ux(jj,kk,ii)+tau2(1,2)*Uy(jj,kk,ii)+tau2(2,1)*Vx(jj,kk,ii)+tau2(2,2)*Vy(jj,kk,ii);
            Ushear2(jj,kk,ii) = abs(Ushear2(jj,kk,ii));
            Usheartotal(jj,kk,ii) = sqrt((Ushear1(jj,kk,ii))^2+(Ushear2(jj,kk,ii))^2);
        end
    end
    end
end