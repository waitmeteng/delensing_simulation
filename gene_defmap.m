function [defmap, defmap_x, defmap_y] = gene_defmap(reso,angu_size,seed_n,lmax_cut,l,cl_dd)
% this code is to produce deflection field and deflection angles using
% input deflection field power.
%--------------------------------------------------------------------------
% INPUT:
% .reso               : resolution of map
% .angu_size          : angular size of map (degree)
% .seed_n             : random seed of map
% .lmax_cut           : maximum l that not equals to zero
% .l                  : minimum l corresponds to input power
% 
% OUTPUT:
% .defmap             : deflection field
% .defmap_x           : deflection angle in x direction
% .defmap_y           : deflection angle in y direction
%--------------------------------------------------------------------------
% Wei-Hsiang Teng, NTU, 2010
%--------------------------------------------------------------------------

% [fid,message]=fopen('camb_r005_scalcls.dat','r');
% if fid==-1
%     disp(message);
% else
%     disp(fid);
% end
% power_spectrum=fscanf(fid,'%g',[6,9999]);
% cl_dd=(power_spectrum(1,:)+1)./(power_spectrum(1,:).^3).*power_spectrum(5,:);
% cl_dd=[0,0,cl_dd];
% fclose(fid);
% %--------------------parameter---------------------------------------------
% reso=300; % pixels number
% angu_size=10; % size of map
% seed_n=100; 
% file_choice=0; % 0:binary, 1:ASCII
% %--------------------------------------------------------------------------

angu_size_rad=angu_size*pi/180;
area = angu_size_rad*angu_size_rad;
% reso_s=num2str(reso);
% angu_s=num2str(angu_size);

delta_x=angu_size/reso*pi/180;
d_l=360/angu_size;

if (2*lmax_cut >= reso/2*d_l)
   disp('warning: the resolution is too small that may induce aliasing problem!');  
end

l_m=zeros(reso/2+1,reso);
lx_m=zeros(reso/2+1,reso);
ly_m=zeros(reso/2+1,reso);
def_l=zeros(reso,reso);
cl_ddm=zeros(reso/2+1,reso);
for i=1:reso
    for j=1:reso/2+1
        kx=i-(reso/2+1);
        ky=j-(reso/2+1);    
        ll=round(sqrt(kx^2+ky^2)*d_l); 
        if (ll<=lmax_cut)
           cl_ddm(j,i)=cl_dd(ll-l+1);
        end
        l_m(j,i)=ll;
        lx_m(j,i)=kx*d_l;
        ly_m(j,i)=ky*d_l;
    end
end
norm_l=(cl_ddm.^0.5)*sqrt(area)/sqrt(2);
%-------------------------generate deflection field------------------------

seed_s=num2str(seed_n);
randn('seed',seed_n);
def_l(1:reso/2+1,:)=(randn(reso/2+1,reso)+1i*randn(reso/2+1,reso)).*norm_l;
def_l(1,1)=real(def_l(1,1));
def_l(1,reso/2+1)=real(def_l(1,reso/2+1));
def_l(reso/2+1,1)=real(def_l(reso/2+1,1));
def_l(reso/2+2:reso,1)=conj(def_l(reso/2:-1:2,1));
def_l(1,reso/2+2:reso)=conj(def_l(1,reso/2:-1:2));
def_l(reso/2+1:reso,reso/2+1:reso)=conj(def_l(reso/2+1:-1:2,reso/2+1:-1:2));
def_l(reso/2+2:reso,2:reso/2)=conj(def_l(reso/2:-1:2,reso:-1:reso/2+2));
defmap=ifft2(ifftshift(def_l))./(delta_x^2);
def_lx=1i*lx_m./l_m.*def_l(1:reso/2+1,:);
def_ly=1i*ly_m./l_m.*def_l(1:reso/2+1,:);
def_lx(reso/2+1,reso/2+1)=0;
def_ly(reso/2+1,reso/2+1)=0;
def_lx(1,1)=real(def_lx(1,1));
def_lx(1,reso/2+1)=real(def_lx(1,reso/2+1));
def_lx(reso/2+1,1)=real(def_lx(reso/2+1,1));
def_lx(reso/2+2:reso,1)=conj(def_lx(reso/2:-1:2,1));
def_lx(1,reso/2+2:reso)=conj(def_lx(1,reso/2:-1:2));
def_lx(reso/2+1:reso,reso/2+1:reso)=conj(def_lx(reso/2+1:-1:2,reso/2+1:-1:2));
def_lx(reso/2+2:reso,2:reso/2)=conj(def_lx(reso/2:-1:2,reso:-1:reso/2+2));
def_ly(1,1)=real(def_ly(1,1));
def_ly(1,reso/2+1)=real(def_ly(1,reso/2+1));
def_ly(reso/2+1,1)=real(def_ly(reso/2+1,1));
def_ly(reso/2+2:reso,1)=conj(def_ly(reso/2:-1:2,1));
def_ly(1,reso/2+2:reso)=conj(def_ly(1,reso/2:-1:2));
def_ly(reso/2+1:reso,reso/2+1:reso)=conj(def_ly(reso/2+1:-1:2,reso/2+1:-1:2));
def_ly(reso/2+2:reso,2:reso/2)=conj(def_ly(reso/2:-1:2,reso:-1:reso/2+2));
defmap_x=ifft2(ifftshift(def_lx))./(delta_x^2);
defmap_y=ifft2(ifftshift(def_ly))./(delta_x^2);

%-----------------------save the data--------------------------------------
% file_name=strcat('./deflection_field/defmap_',reso_s,'_',angu_s);
% file_namex=strcat('./deflection_field/defmapx_',reso_s,'_',angu_s);
% file_namey=strcat('./deflection_field/defmapy_',reso_s,'_',angu_s);
% if (file_choice==0)
%     file_name0=strcat(file_name,'.bin',seed_s);
%     fid0=fopen(file_name0,'w');
%     count=fwrite(fid0,def_map,'double');
%     fclose(fid0);
%     if (count~=reso^2)
%        disp('something wrong happened when saving the data!!!'); 
%     end
%     file_namex0=strcat(file_namex,'.bin',seed_s);
%     fidx0=fopen(file_namex0,'w');
%     count=fwrite(fidx0,def_x,'double');
%     fclose(fidx0);
%     if (count~=reso^2)
%        disp('something wrong happened when saving the data!!!'); 
%     end
%     file_namey0=strcat(file_namey,'.bin',seed_s);
%     fidy0=fopen(file_namey0,'w');
%     count=fwrite(fidy0,def_y,'double');
%     fclose(fidy0);
%     if (count~=reso^2)
%        disp('something wrong happened when saving the data!!!'); 
%     end
% elseif (file_choice==1)
%     file_name1=strcat(file_name,'.dat',seed_s);
%     fid1=fopen(file_name1,'w');
%     fprintf(fid1,'%10.5e\n',def_map);
%     fclose(fid1);
%     
%     file_namex1=strcat(file_namex,'.dat',seed_s);
%     fidx1=fopen(file_namex1,'w');
%     fprintf(fidx1,'%10.5e\n',def_x);
%     fclose(fidx1);
%     
%     file_namey1=strcat(file_namey,'.dat',seed_s);
%     fidy1=fopen(file_namey1,'w');
%     fprintf(fidy1,'%10.5e\n',def_y);
%     fclose(fidy1);
% end