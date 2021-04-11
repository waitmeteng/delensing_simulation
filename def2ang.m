function [ dx_angle, dy_angle ] = def2ang(angu_size, defmap)
% produce deflection angles in x, y directions.
%--------------------------------------------------------------------------
% INPUT:
% .angu_size                     : angular size of map
% .defmap                        : deflection field
% OUTPUT: 
% .dx_angle                      : deflection angle in x direction
% .dy_angle                      : deflection angle in y direction
%--------------------------------------------------------------------------
reso=size(defmap,1);
delta_x=angu_size/reso*pi/180;
d_l=360/angu_size;
l_m=zeros(reso/2+1,reso);
lx_m=zeros(reso/2+1,reso);
ly_m=zeros(reso/2+1,reso);
for i=1:reso
    for j=1:reso/2+1
        kx=i-(reso/2+1);
        ky=j-(reso/2+1);    
        ll=round(sqrt(kx^2+ky^2)*d_l);        
        l_m(j,i)=ll;
        lx_m(j,i)=kx*d_l;
        ly_m(j,i)=ky*d_l;
    end
end
%--------------------------------------------------------------------------
def_l=fftshift(fft2(defmap))*delta_x^2; 
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
dx_angle=ifft2(ifftshift(def_lx))./(delta_x^2);
dy_angle=ifft2(ifftshift(def_ly))./(delta_x^2);

end

