function [defmap_mean]=get_meanfield(resol, angu_size, seed_n, lmax, lmin, lmin0, cl_ee, cl_bb, cl_np, cl_dd, al, number, dx_angler, dy_angler)
% this code is to generate the mean field in equation (66) of 
% "Reconstruction of lensing from the cosmic microwave background 
% polarization", arXiv:astro-ph/0306354.
% -------------------------------------------------------------------------
reso=size(dx_angler,1);
if (size(dy_angler,1)~=reso)
   disp('The size of input maps are not the same!'); return 
end
defmap_mean=zeros(reso,reso);
for seed=seed_n:seed_n+number-1
[ emap_n ] = gene_map( resol, angu_size, seed, lmax, 0, cl_np);
[ bmap_n ] = gene_map( resol, angu_size, seed+2000, lmax, 0, cl_np);
[qmap_noise, umap_noise] = eb2qu(angu_size, emap_n, bmap_n);
[qmap_noise1]=lensing1(resol,angu_size,dx_angler,dy_angler,qmap_noise);
[umap_noise1]=lensing1(resol,angu_size,dx_angler,dy_angler,umap_noise);
% [tmap_noise1,qmap_noise1,umap_noise1]=lensing(resol,angu_size, ...
%      dx_angler,dy_angler,zeros(resol,resol),qmap_noise,umap_noise);
% [qmap_noise1] = lensing_harmonic1(angu_size, lmax, qmap_noise, defmap);
% [umap_noise1] = lensing_harmonic1(angu_size, lmax, umap_noise, defmap);
[defmap_mean1]=ebee_estimatorclnn(angu_size,0,lmax,lmin,lmin0,al,qmap_noise1,umap_noise1,0, ...
    cl_ee,cl_ee,cl_bb,cl_bb,cl_dd,cl_np,cl_np,'',al);
defmap_mean=defmap_mean+defmap_mean1;
end
defmap_mean=defmap_mean/number;
end

