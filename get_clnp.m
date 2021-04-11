function [ cl_np ] = get_clnp( fwhm, noise, lmax, t_cmb1)
fwhm_rad=fwhm/60*pi/180;
noise1=noise/(180/pi*60);
l=0:1:lmax;
if (noise~=0)
cl_np=(noise1/t_cmb1)^2*exp(l.*(l+1)*(fwhm_rad^2)/8/log(2));
else
cl_np=zeros(1,lmax+1);    
end

