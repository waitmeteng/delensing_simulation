tic
% delensing_forecast_Smith.m
% -------------------------------------------------------------------------
% This script is to calculate the improvement alpha=sigma_0(r)/sigma(r)
% in the statistical error on r due to polarization delensing, for 
% varying noise level and beam.
% 
% The algorithm and formalism are defined in 
% Smith, Hanson, Loverde, Hirata & Zahn (2010)
% "arXiv:1010.0048, Delensing CMB Polarization with External Datasets".
% -------------------------------------------------------------------------
%
% Wei-Hsiang Teng, NTU, Dec 2010
% 
%--------------------------------------------------------------------------
% In this script you have to give the 
% . range of noise level                 : NOISE_LEVEL
%   for polarization (micronK-arcmin)
% . beam size (arcmin)                   : BEAM
% . whether to use iterative method      : ITERATE
%   0: not use    1: use
% . the maximum l (l<7000)               : LMAX
% . the percentage which stops           : STOP_RATIO
%   iteration
% . the maximum l for averaging          : LCUT
% . whether to use flat approximation    : FLAT
%   0: not use    1: use
%--------------------------------------------------------------------------
% NOISE_LEVEL=[0.1:0.1:0.4,0.5:0.25:3,3.5:0.5:6];
NOISE_LEVEL=1.41;
BEAM=5;
ITERATE=0;
LMAX=4000;
STOP_RATIO=1;
LCUT=180;
FLAT=0;
if (FLAT==1)
   lreso=20; 
   larray=lreso+1:LMAX;
   larray=larray';
end
function [cl_residual]=get_residual(noise_p, BEAM, LMAX, STOP_RATIO, cl_phi, ...
    cl_eeg, cl_ee, cl_bbg, cl_bb, d33, d3m3, d31, d3m1, d11, d1m1, d22, d2m2)

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% start forecasting !!!
%--------------------------------------------------------------------------
lmun=LCUT-1;
if (FLAT==0)   % using spherical harmonic to forecast

if (ITERATE==1)
cl_resav1=sum(cl_bb(3:LCUT+1))/lmun;
cl_resav0=0;
cl_res1=cl_bb;
while (abs(cl_resav1-cl_resav0)/cl_resav1*100>STOP_RATIO)
cl_nn=cal_recon_noise(LMAX,cl_ee,cl_res1,cl_np,cl_np ...
    ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
cl_nn1=[0;cl_nn];
[cl_out,cl_res]=cal_res_bpower(LMAX,cl_phi,cl_eeg,cl_ee,cl_np,cl_nn1 ...
    ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
cl_res1=[0;0;cl_res];
cl_resav0=cl_resav1;
cl_resav1=sum(cl_res1(3:LCUT+1))/lmun;
end
elseif (ITERATE==0)
cl_nn=cal_recon_noise(LMAX,cl_ee,cl_bb,cl_np,cl_np ...
    ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
cl_nn1=[0;cl_nn];
[cl_out,cl_res]=cal_res_bpower(LMAX,cl_phi,cl_eeg,cl_ee,cl_np,cl_nn1 ...
    ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
cl_res1=[0;0;cl_res];
end
alpha_l=(cl_bb(3:LMAX+1)+cl_np(3:LMAX+1))./(cl_res1(3:LMAX+1)+cl_np(3:LMAX+1));
alpha(kk)=sum(alpha_l(1:LCUT-1))/lmun;

% using Fourier (flat approximation) to forecast
elseif (FLAT==1)

if (ITERATE==1)
cl_resav1=sum(cl_bb(3:LCUT+1))/lmun;
cl_resav0=0;
cl_res1=cl_bb;
while (abs(cl_resav1-cl_resav0)/cl_resav1*100>STOP_RATIO)
[cl_nnt,lout1]=calreconnoise_flat(LMAX,lreso,cl_ee,cl_res1,cl_np,cl_np);
lout=lout1(2:size(lout1,2));
cl_nn=cl_nnt(2:size(lout1,2));
cl_nn1=interp1(log10(lout),log10(cl_nn.*lout.^3.*(lout+1)/2/pi),log10(larray),'cubic');
cl_nn2=[zeros(1,lreso+1)';(10.^cl_nn1)./(larray.^3)./(larray+1)*2*pi];
[cl_res,lout]=calresbpower_flat(LMAX,lreso,cl_eeg,cl_np,cl_phi,cl_nn2);
cl_rest=interp1(lout,cl_res,larray,'cubic');
cl_res1=[0;0;ones(lreso-1,1)*cl_rest(larray(1),1);cl_rest];
cl_resav0=cl_resav1;
cl_resav1=sum(cl_res1(3:LCUT+1))/lmun;

end
elseif (ITERATE==0)
[cl_nnt,lout1]=calreconnoise_flat(LMAX,lreso,cl_ee,cl_bb,cl_np,cl_np);
lout=lout1(2:size(lout1,2));
cl_nn=cl_nnt(2:size(lout1,2));
cl_nn1=interp1(log10(lout),log10(cl_nn.*lout.^3.*(lout+1)/2/pi),log10(larray),'cubic');
cl_nn2=[zeros(1,lreso+1)';(10.^cl_nn1)./(larray.^3)./(larray+1)*2*pi];
[cl_res,lout]=calresbpower_flat(LMAX,lreso,cl_eeg,cl_np,cl_phi,cl_nn2);
% cl_res1=interp1(lout,cl_res,larray,'cubic');
cl_rest=interp1(lout,cl_res,larray,'cubic');
cl_res1=[0;0;ones(lreso-1,1)*cl_rest(larray(1),1);cl_rest];
end
% top=find(lout>=LCUT);
% res_sum=sum(cl_res1(2:top(1)))/size(cl_res1(2:top(1)),2);
% lens_sum=sum(cl_bb(3:LCUT+1)+cl_np(3:LCUT+1))/size(cl_bb(3:LCUT+1),1);
% alpha(kk)=lens_sum/res_sum;
alpha_l=(cl_bb(3:LMAX+1)+cl_np(3:LMAX+1))./(cl_res1(3:LMAX+1)+cl_np(3:LMAX+1));
alpha(kk)=sum(alpha_l(1:LCUT-1))/lmun;
% end of forecasting !!!

end
end
toc