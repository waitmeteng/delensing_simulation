% This code is to simulate the method of "moving notch filtering" in the
% case there is primordial B modes at high-l region (in our case is vector
% B mode from cosmic string). The estimated primordial B power will be
% given.
%--------------------------------------------------------------------------

tic
%---------------parameters-------------------------------------------------
frac=1;
reso=600*frac; % resolution of unlensed maps
lmax=5400;
resol=600*frac; % resolution of lensed maps
angu_size=10*frac;
fwhm=4;        % FWHM in unit of arcmin
begin=3;    
final=7;    % final-begin+1 is the number of map
noise_p=1.41;  % noise level in unit of microK arcmin
tcmb=2.725e6;  % microK
lreso=360/angu_size;   % the spacing of multipole number
%--------------binning-----------------------------------------------------
bin=[36,180,288:108:1044];
bin_n1=size(bin,2)-1;        % number of bins
%--------------parameters for iterative estimator--------------------------
%----for details see Hirata & Seljak's paper-------------------------------
con_factor=0.12;        % convergence factor
mean_number=3;          % the number of maps for mean field
iter_number=64;         % to do 'iter_number' times of iteration
iter_s=num2str(iter_number);
do_iter=1;              % if perform the iterative estimator (1:yes, otherwise:no)
%--------------------------------------------------------------------------
r=0;                    % if r is not zero, include the tensor fields    
gmu=1.0e-7;             % if gmu is not zero, include the vector fields
esti_iter=0;
esti_num=1;
save_file=1;            % if save the results
%-----------------initialize the file names--------------------------------
reso_s=num2str(reso);
angu_s=num2str(angu_size);
noise_s=num2str(noise_p);
begin_s=num2str(begin);
r_s=num2str(r);
gmu_s=num2str(gmu);
fwhm_s=num2str(fwhm);
num=final-begin+1;
num_s=num2str(num);
file_nameq=strcat('./result/estimated_prim_bpower_quad_r',r_s,'_gmu',gmu_s,'_',num_s ...
    ,'maps_reso',reso_s,'_angu',angu_s,'_noise',noise_s,'_fwhm',fwhm_s,'_start',begin_s,'.txt');
file_namei=strcat('./result/estimated_prim_bpower_iter_r',r_s,'_gmu',gmu_s,'_',num_s ...
    ,'maps_reso',reso_s,'_angu',angu_s,'_noise',noise_s,'_fwhm',fwhm_s,'_start',begin_s,'.txt');
%------------------read the power spetra-----------------------------------

lnumber=9999;
[l,cl_tt, cl_ee, cl_te, cl_phi]= ...
    read_cambpower('./data/camb_r005_scalcls.dat', 'scalar',lnumber );

cl_tt=[0,0,cl_tt./l./(l+1)*2*pi];
cl_tt=cl_tt(1:lmax+1);
cl_ee=[0,0,cl_ee./l./(l+1)*2*pi];
cl_ee=cl_ee(1:lmax+1);
cl_te=[0,0,cl_te./l./(l+1)*2*pi];
cl_te=cl_te(1:lmax+1);
cl_dd=(l+1)./(l.^3).*cl_phi;
cl_phi=cl_dd./l./(l+1);
l=[0,1,l];
cl_dd=[0,0,cl_dd(1:lmax-1)];
cl_phi=[0,0,cl_phi(1:lmax-1)];
cl_bb=zeros(1,lmax+1);

[ll,cl_ttl, cl_eel, cl_bbl, cl_tel]= ...
    read_cambpower('./data/camb_r005_lensedcls.dat', 'lensed',9283);
cl_eel=[0,0,cl_eel./ll./(ll+1)*2*pi];
cl_eel=cl_eel(1:lmax+1);
cl_bbl=[0,0,cl_bbl./ll./(ll+1)*2*pi];
cl_bbl=cl_bbl(1:lmax+1);

if (r~=0)
[lll,cl_ttt, cl_eet, cl_bbt, cl_tet]= ...
    read_cambpower('./data/camb_r005_tenscls.dat', 'tensor',lnumber);
cl_ttt=[0,0,cl_ttt./lll./(lll+1)*2*pi]/0.05*r;
cl_eet=[0,0,cl_eet./lll./(lll+1)*2*pi]/0.05*r;
cl_bbt=[0,0,cl_bbt./lll./(lll+1)*2*pi]/0.05*r;
cl_tet=[0,0,cl_tet./lll./(lll+1)*2*pi]/0.05*r;
lll=[0,1,lll];
end

if (gmu~=0)
[l4,cl_tts, cl_ees, cl_bbs, cl_tes]= ...
    read_cambpower('./data/cmbact1_gmu_eminus7.dat', 'tensor',2999);
cl_tts=[0,0,cl_tts./l4./(l4+1)*2*pi]/(1.0e-7)*gmu;
cl_ees=[0,0,cl_ees./l4./(l4+1)*2*pi]/(1.0e-7)*gmu;
cl_bbs=[0,0,cl_bbs./l4./(l4+1)*2*pi]/(1.0e-7)*gmu;
cl_tes=[0,0,cl_tes./l4./(l4+1)*2*pi]/(1.0e-7)*gmu;
l4=[0,1,l4];
end

cl_np=get_clnp(fwhm,noise_p,lmax,tcmb);

%--------------------------------------------------------------------------
% produce the renormalization factor F_l
%---------------------for quadratic estimator -----------------------------
[al_eb] = gene_al_eb(resol,angu_size,lmax,cl_ee,cl_eel,cl_bb,cl_bbl,cl_np,cl_np);
al=zeros(1,lmax+1);
linput=0:lreso:lreso*(resol/2-1);
lout=lreso:lmax;
al_ebt=al_eb.*linput.*(linput+1)/2/pi;
al(lout+1)=interp1(linput,al_ebt,lout,'cubic');
al(lout+1)=al(lout+1)./lout./(lout+1)*2*pi;

[al_ee] = gene_al_ee(resol,angu_size,lmax,cl_ee,cl_eel,cl_np);
al1=zeros(1,lmax+1);
linput=0:lreso:lreso*(resol/2-1);
lout=lreso:lmax;
al_eet=al_ee.*linput.*(linput+1)/2/pi;
al1(lout+1)=interp1(linput,al_eet,lout,'cubic');
al1(lout+1)=al1(lout+1)./lout./(lout+1)*2*pi;

al2=1./(1./al1+1./al);
al2(1:lreso)=0;
%--------------------------------------------------------------------------
% for iterative estimator
[al_eb] = gene_al_eb(resol,angu_size,lmax,cl_ee,cl_ee,cl_bb,cl_bb,cl_np,cl_np);
al=zeros(1,lmax+1);
linput=0:lreso:lreso*(resol/2-1);
lout=lreso:lmax;
al_ebt=al_eb.*linput.*(linput+1)/2/pi;
al(lout+1)=interp1(linput,al_ebt,lout,'cubic');
al(lout+1)=al(lout+1)./lout./(lout+1)*2*pi;

[al_ee] = gene_al_ee(resol,angu_size,lmax,cl_ee,cl_ee,cl_np);
al1=zeros(1,lmax+1);
linput=0:lreso:lreso*(resol/2-1);
lout=lreso:lmax;
al_eet=al_ee.*linput.*(linput+1)/2/pi;
al1(lout+1)=interp1(linput,al_eet,lout,'cubic');
al1(lout+1)=al1(lout+1)./lout./(lout+1)*2*pi;

al2i=1./(1./al1+1./al);
al2i(1:lreso)=0;
% -----------------calculate the needed Wigner-d---------------------------
% LMAX=3000;
% [weight,z,N]=gausslegendrecof(3*LMAX,'jacobi',[-1 1]);
% nanw=isnan(weight);
% infw=isinf(weight);
% ini=1;
% jni=1;
% for k=1:N
% if (nanw(k)==1 || infw(k)==1)
%     xi(ini)=k;
%     ini=ini+1;
% else
%     x(jni)=k;
%     y(jni)=weight(k);
%     jni=jni+1;
% end
% end
% % Because there is problem to calculate weights for large l
% % I use interpolation to generate some weights.
% weight(xi)=interp1(x,y,xi,'spline');
% thetav=acos(z)*180/pi;
% %-------------------------------------------------------------------------- 
% % calculate the needed Wigner d function
% % d^l_33, d^l_3-3 (l>=3)
% % d^l_31, d^l_3-1 (l>=3)
% % d^l_11, d^l_1-1 (l>=1)
% % d^l_22, d^l_2-2 (l>=2)
% % --------------------------------------------------------------------------
% d33=blanco_m(LMAX,3,3,thetav')';
% d3m3=blanco_m(LMAX,3,-3,thetav')';
% d31=blanco_m(LMAX,3,1,thetav')';
% d3m1=blanco_m(LMAX,3,-1,thetav')';
% d11=blanco_m(LMAX,1,1,thetav')';
% d1m1=blanco_m(LMAX,1,-1,thetav')';
% d22=blanco_m(LMAX,2,2,thetav')';
% d2m2=blanco_m(LMAX,2,-2,thetav')';

%----------------let's start !!!!!!!!!!!!!!!!!!!---------------------------
bpower1=zeros(num,bin_n1);
if (do_iter==1)
  bpower2=zeros(num,bin_n1);  
end
if (esti_iter==1)
  bpower3=zeros(num,bin_n1);  
end
for seed_n=begin:final
% produce the deflection angle
seed_d=seed_n;    % random seed 
[defmap, dx_angle, dy_angle] = ...
    gene_defmap(reso,angu_size,seed_d,lmax,l(1),cl_dd);  
% produce the unlensed fields
seed_n1=seed_n+2^25;   % random seed    
[tmap_unlens,qmap_unlens,umap_unlens]=gene_tqu(reso,angu_size,seed_n1, ...
    lmax,l(1),cl_tt,cl_ee,cl_bb,cl_te);
% produce the tensor fields if r is not zero
if (r~=0)
    seed_nr=seed_n+2^26;   % random seed 
    [tmapr,qmapr,umapr]=gene_tqu(reso,angu_size,seed_nr, ...
    lmax,lll(1),cl_ttt,cl_eet,cl_bbt,cl_tet);
    tmap_unlens=tmap_unlens+tmapr;
    qmap_unlens=qmap_unlens+qmapr;
    umap_unlens=umap_unlens+umapr;
end

if (gmu~=0)
    seed_ns=seed_n+2^26+2^25;   % random seed 
    [tmaps,qmaps,umaps]=gene_tqu(reso,angu_size,seed_ns, ...
    3000,l4(1),cl_tts,cl_ees,cl_bbs,cl_tes);
    tmap_unlens=tmap_unlens+tmaps;
    qmap_unlens=qmap_unlens+qmaps;
    umap_unlens=umap_unlens+umaps;
end
% produce the lensed fields by remapping of unlensed fields
[qmap_lens, umap_lens] = gene_lensedmapnott(resol, 1, angu_size, ...
    0,dx_angle,dy_angle,qmap_unlens,umap_unlens);
clear tmap_unlens qmap_unlens umap_unlens 
% produce the noise if noise_p is not zero
if (noise_p~=0)      
seed_ne=seed_n+2^27;
seed_nb=seed_n+2^27+2^25;
[ emap_n ] = gene_map( resol, angu_size, seed_ne, lmax, 0, cl_np);
[ bmap_n ] = gene_map( resol, angu_size, seed_nb, lmax, 0, cl_np);
[qmap_noise, umap_noise] = eb2qu(angu_size, emap_n, bmap_n);
else 
   qmap_noise=zeros(resol,resol);
   umap_noise=zeros(resol,resol);
end
clear emap_n bmap_n
% the simulated observation maps
qmap_obs=qmap_lens+qmap_noise;
umap_obs=umap_lens+umap_noise;

%----------------start the "moving notch filter"---------------------------
bmap_delens1=zeros(reso,reso); % for quadratic
bmap_delens2=zeros(reso,reso); % for iterative
for k=1:bin_n1
%----------use quadratic estimator to get deflection field-----------------
[defmapr]=ebee_estimatorclnn(angu_size,0,lmax,bin(k+1),bin(k),al2,qmap_obs, ...
    umap_obs,0,cl_ee,cl_eel,cl_bb,cl_bbl,cl_dd,cl_np,cl_np,'wiener',al2);
defmapq=defmapr;
%-----------------------get delensed B map---------------------------------
[emap_obs, bmap_obs] = qu2eb(angu_size, lmax, qmap_obs, umap_obs, 0);
[ emap_w ] = wiener_filter( angu_size, lmax, emap_obs, cl_ee, cl_np);
[bmap_delens] = lensing_harmonic(angu_size, lmax, bmap_obs, emap_w, defmapq);
[bmap_delens] = get_annulus_map( angu_size, bin(k), bin(k+1), bmap_delens);
bmap_delens1=bmap_delens1+bmap_delens;
%--------------------------------------------------------------------------
%-------------------------start the iterative procedure--------------------
if (do_iter==1)
    
for kk=1:iter_number
seed_m=kk+2^24;
%--------------------------------------------------------------------------
% % produce the mean field due to instrumental noise
% % the order of a few realizations is sufficient.
[ dx_angler, dy_angler ] = def2ang(angu_size, defmapr);
[defmap_mean]=get_meanfield(resol, angu_size, seed_m, lmax, bin(k+1), bin(k), ...
    cl_ee, cl_bb, cl_np, cl_dd, al2i, mean_number, -dx_angler,-dy_angler);
[qmap_obs1]=lensing1(resol,angu_size,-dx_angler,-dy_angler,qmap_obs);
[umap_obs1]=lensing1(resol,angu_size,-dx_angler,-dy_angler,umap_obs);
%--------------------------------------------------------------------------
[defmapri]=ebee_estimatorclnn(angu_size,0,lmax,bin(k+1), bin(k), al2i,qmap_obs1, ...
    umap_obs1,0,cl_ee,cl_ee,cl_bb,cl_bb,cl_dd,cl_np,cl_np,'',al2i);
[defmap_wi] = wiener_filter( angu_size, lmax, defmapri-defmap_mean, ...
    cl_dd(1:lmax+1)*con_factor, al2i+cl_dd(1:lmax+1)*(1-con_factor));
[defmap_w] = wiener_filter( angu_size, lmax, defmapr, ...
    al2i*(1-con_factor)+cl_dd(1:lmax+1), al2i*con_factor);
[defmap_wi]=lensing1(resol,angu_size,dx_angler,dy_angler,defmap_wi);
diff_map=defmap_wi+defmap_w-defmapr;
defmapr=defmap_wi+defmap_w;

end
%----------------------to get delensed B map-------------------------------
[bmap_delens] = lensing_harmonic(angu_size, lmax, bmap_obs, emap_w, defmapr);
[bmap_delens] = get_annulus_map( angu_size, bin(k), bin(k+1), bmap_delens);
bmap_delens2=bmap_delens2+bmap_delens;
%--------------------------------------------------------------------------
end
end%-------------------end of "moving notch filter"------------------------
% calculate the delensing B power 
[l_arrayb1,bpower1(seed_n-begin+1,:)]=get_power(angu_size,0,bmap_delens1,bin_n1,bin);
if (do_iter==1)
[l_arrayb2,bpower2(seed_n-begin+1,:)]=get_power(angu_size,0,bmap_delens2,bin_n1,bin);    
end

% % start iteration to get unbiased estimation of primordial B power
%     cl_pri=zeros(1,lmax+1);
%     %cl_pri(1:3001)=cl_bbs;
%     cl_temp=bpower1(seed_n-begin+1,:);
%     %if (esti_iter~=1)
%         %cl_pri=zeros(1,lmax+1);
%     cl_nn=cal_recon_noise(LMAX,cl_ee',cl_bbl'+cl_pri',cl_np',cl_np' ...
%     ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
%     cl_nn1=[0;cl_nn];
%     [cl_res]=cal_res_bpower(LMAX,cl_phi',cl_ee',cl_eel',cl_np',cl_nn1, ...
%         d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
%     cl_res1=(cl_res(round(l_arrayb1)-1)'+cl_np(round(l_arrayb1)+1)).*l_arrayb1.*(l_arrayb1+1)/2/pi;
%     bpower1(seed_n-begin+1,:)=bpower1(seed_n-begin+1,:)-cl_res1;
%     else
%         for it=1:esti_num+1
%             cl_pri1=cl_pri;
%             
%     cl_nn=cal_recon_noise(LMAX,cl_ee',cl_bbl'+cl_pri1',cl_np',cl_np' ...
%     ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
%     cl_nn1=[0;cl_nn];
%     [cl_res]=cal_res_bpower(LMAX,cl_phi',cl_ee',cl_eel',cl_np',cl_nn1, ...
%         d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
%     cl_res1=(cl_res(round(l_arrayb1)-1)'+cl_np(round(l_arrayb1)+1)).*l_arrayb1.*(l_arrayb1+1)/2/pi;
%     cl_est=bpower1(seed_n-begin+1,:)-cl_res1;
%     ltemp=ceil(l_arrayb1(1)):floor(l_arrayb1(bin_n1));
%     cl_prig=interp1(l_arrayb1,cl_est,ltemp,'linear');
%     cl_pri=[zeros(1,ceil(l_arrayb1(1))),cl_prig,zeros(1,lmax+1-ceil(l_arrayb1(1))-size(cl_prig,2))];
%     cl_pri=cl_pri./lll./(lll+1)*2*pi;
%         end
%     bpower3(seed_n-begin+1,:)=cl_est;    
%     end

end%-------------------end of simulation-----------------------------------

%--------------------------saving the results------------------------------
if (save_file==1)
fid1=fopen(file_nameq,'w');
fprintf(fid1,'%4d %10.5e ',l_arrayb1,bpower1');
fclose(fid1);

if (do_iter==1)
fid1=fopen(file_namei,'w');
fprintf(fid1,'%4d %10.5e ',l_arrayb2,bpower2');
fclose(fid1);
end
end
toc