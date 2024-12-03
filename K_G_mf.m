function [Kmf,Gmf] = K_G_mf(Kh,Kdry,Gdry,omegia,Viscosity,faic,alpha,Kg,Kf,Sc,x)

%====================湿骨架模量计算==============饱一种流体！！！！！

%              x=1:正常      x=2：高频极限   x=3: 低频极限

%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%  Kh:高压下软孔隙闭合而硬孔隙未闭合的骨架体积模量
%  kg:颗粒体积模量
%  wr: 软孔隙中的水占岩石孔隙中所有水的比例，称为润湿比
%  Sc :  软孔隙中的水饱和度   (最原始喷射流Sc恒为1)

Ka=1./alpha*(-3i*omegia.*Viscosity/Kf).^(1/2);
Kf_1=(1-2*besselj(1,Ka)./(Ka.*besselj(0,Ka))).*Kf;

if x==1
    Kmf=1./Kh+1./(1./(1./Kdry-1./Kh)+Sc./(1./Kf_1-1/Kg)./faic);
elseif x==2
    Kmf=1./Kh+1./(1./(1./Kdry-1./Kh)+Sc./(1./Kf-1/Kg)./faic);
else
    Kmf=Kdry;
end

Kmf=1./Kmf;

miumf=1./Gdry-4/15*(1./Kdry-1./Kmf);
Gmf=1./miumf;
end

