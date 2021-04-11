function [qmap_lensb, umap_lensb] = gene_lensedmapnott(resol, number, angu_size, ...
    fwhm,dx_angle,dy_angle,qmap,umap)
% this code is to generate lensed T,Q,U maps according to the input
% unlensed power spectra, resolution, angular size and deflection angles.
% The flat approximation has been used which means we use Fourier transform 
% instead of spherical harmonic transform. And the boundary condition has
% been applied.
%--------------------------------------------------------------------------
% INPUT:
% .resol                 : resolution of lensed maps
% .resou                 : resolution of unlensed maps and deflection
%                          field. This value should be smaller than
%                          resol and resol/resou should be integer.
% .angu_size             : angular size of map
% .fwhm                  : FWHM (arcmin), if input zero no beam effect
%                          applys to  the lensed maps
% .dx_angle              : deflection angle in x direction
% .dy_angle              : deflection angle in y direction
% 
% OUTPUT:
% .tmap_lensed           : lensed T map
% .qmap_lensed           : lensed Q map
% .umap_lensed           : lensed U map
%--------------------------------------------------------------------------
% Wei-Hsiang Teng, NTU, 2010
%------------------------------parameter-----------------------------------
% the number of grid that 1 sigma (beam size) will take 
if (fwhm ~=0)
beamsigma=resol*fwhm/2.35482/60/angu_size;
reso_g=round(8*beamsigma);
if (floor(reso_g/2)~=reso_g/2)
    reso_g=reso_g-1;
end
[fx,fy]=meshgrid(-reso_g/2:reso_g/2-1);
norm=1/2/pi/beamsigma/beamsigma;
% the beam
g=norm*exp((-fx.^2-fy.^2)/2/beamsigma/beamsigma);
end

[qmap_lens,umap_lens]=lensing_nott(resol,number, angu_size,dx_angle,dy_angle, ...
           qmap,umap);
if (fwhm~=0)
% tmap_lensb = imfilter(tmap_lens,g,'circular','conv');
qmap_lensb = imfilter(qmap_lens,g,'circular','conv');
umap_lensb = imfilter(umap_lens,g,'circular','conv');
else
% tmap_lensb = tmap_lens;
qmap_lensb = qmap_lens;
umap_lensb = umap_lens;    
end
