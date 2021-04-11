function [cl_res1]=get_residual(noise_p, tcmb, BEAM, LMAX, ITERATE, cl_phi, ...
    cl_eeg, cl_ee, cl_bbg, cl_bb, d33, d3m3, d31, d3m1, d11, d1m1, d22, d2m2, weight)

%--------------------------------------------------------------------------
LCUT=200;
lmun=LCUT-1;
STOP_RATIO=1;
cl_np=get_clnp(BEAM,noise_p,LMAX,tcmb);
cl_np=cl_np';
%--------------------------------------------------------------------------
% start calculating !!!
%--------------------------------------------------------------------------

if (ITERATE==1)
cl_resav1=sum(cl_bb(3:LCUT+1))/lmun;
cl_resav0=0;
cl_res1=cl_bb;
while (abs(cl_resav1-cl_resav0)/cl_resav1*100>STOP_RATIO)
cl_nn=cal_recon_noise(LMAX,cl_ee,cl_res1+cl_bbg,cl_np,cl_np ...
    ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
cl_nn1=[0;cl_nn];
[cl_res]=cal_res_bpower(LMAX,cl_phi,cl_eeg,cl_ee,cl_np,cl_nn1 ...
    ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
cl_res1=[0;0;cl_res];
cl_resav0=cl_resav1;
cl_resav1=sum(cl_res1(3:LCUT+1))/lmun;
end
elseif (ITERATE==0)
cl_nn=cal_recon_noise(LMAX,cl_ee,cl_bb+cl_bbg,cl_np,cl_np ...
    ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
cl_nn1=[0;cl_nn];
[cl_res]=cal_res_bpower(LMAX,cl_phi,cl_eeg,cl_ee,cl_np,cl_nn1 ...
    ,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight);
cl_res1=[0;0;cl_res];
end
cl_res1=cl_res1+cl_np;
end
