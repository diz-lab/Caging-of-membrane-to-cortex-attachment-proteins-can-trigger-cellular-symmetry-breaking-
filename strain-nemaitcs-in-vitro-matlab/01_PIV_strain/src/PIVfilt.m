function [upix,vpix] = PIVfilt(xpix,ypix,upivpix,vpivpix,snrout,filtspec)

unan = find(isnan(upivpix)==1);
vnan = find(isnan(vpivpix)==1);

if ~isempty(filtspec.snr)
    [su,sv] = snrfilt(xpix,ypix,upivpix,vpivpix,snrout,filtspec.snr);
else
    su = upivpix;
    sv = vpivpix;
end
    
if ~isempty(filtspec.globfiltstd)
    [gu,gv] = globfilt(xpix,ypix,su,sv,filtspec.globfiltstd);
else
    gu = su;
    gv = sv;
end
    
    
if ~isempty(filtspec.locfiltstd)
    [lu,lv] = localfilt(xpix,ypix,gu,gv,filtspec.locfiltstd,filtspec.locfiltmeth,filtspec.kernal);
else
    lu = gu;
    lv = gv;
end


if myIsField(filtspec,'naninterp')==1 && filtspec.naninterp == 0 %interpolate by default
    upix = lu;
    vpix = lv;
else
    [upix,vpix] = naninterp(lu,lv,xpix,ypix);
end



upix(unan) = nan;
vpix(vnan) = nan;

end

