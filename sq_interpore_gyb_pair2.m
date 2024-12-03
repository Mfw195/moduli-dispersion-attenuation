%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; % close all;
Fontsize=13;
linewidth_real = 1.1;
jet_n=jet(9);  CC=colormap(jet(7));
% -------------------------------------------
load('E_dry_gyb.txt');
load('K_dry_gyb.txt');
load('u_dry_gyb.txt');
load('ff.txt');

p=[2.06896551724 5.04310344828 7.11206896552 10.0862068966 15.1293103448 20.0431034483 25.0862068966]';
vp_dry=[3333.333 3517.241 3596.615 3690.304 3766.617 3869.499 3978.159]';
vs_dry=[2276.786 2328.767 2361.111 2405.660 2440.191 2500.000 2524.752 ]';

viso=500;%cp Glycerin 40¡æ%%-----
yita=0.001*viso;
phi=0.132;%porosity 
Kf=4.38*10^9;% Glycerin 40¡æ%%------
den_frame=2305;%kg/m^2     
den_fl=1200;%Glycerin 40¡æ  
den=den_frame+phi*den_fl;
um_st=den_frame*vs_dry.^2;%pa
lamda_st=den_frame*(vp_dry.^2-2*vs_dry.^2);
K_st=lamda_st+2*um_st/3; %bulk modulus  pa
%%%%%%%%%%%%%%%%%
frequencyrange=10.^(-3:0.1:10);
omega=2*pi*frequencyrange;
x=p'*0.001; %Gpa 
x1=x*1000;%MPa
Kd_p=K_st'*1e-9;%Gpa
y=1./Kd_p;%C=1/Kd%1/Gpa
ud_p=um_st'*1e-9;%%Gpa
C_ud_p=1./ud_p;
Kg=45;%Gpa
Kg1=Kg*10^9;
Vg=5.5;%m/s
Vpg=Vg;
Cg=1/Kg; 

fun = @(A)A(1)*(1-A(2)*(A(1)-Cg)*x+A(3)*A(4)*exp(-A(3)*A(1)*x))-y;
t0=[1 2 4000 1.1];
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
t1=lsqnonlin(fun,t0,[],[],options);
xx=0:0.2:40;
xx=xx*0.001; %pa
yy=t1(1)*(1-t1(2)*(t1(1)-Cg)*xx+t1(3)*t1(4)*exp(-t1(3)*t1(1)*xx));% C_K

Cds0=t1(1)+0.002;
theta_s=t1(2)+10;
theta_c=t1(3);
phi_co=t1(4);
%%%%%%%%%%%%%%%%%
fun1 = @(E)E(1)*(1-E(2)*(E(1)-Cg)*x+E(3)*E(4)*exp(-E(3)*E(1)*x))-C_ud_p;
g0=[1.5 2.3 3900 0.9]; 
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
g1=lsqnonlin(fun1,g0,[],[],options);
yy3=g1(1)*(1-g1(2)*(g1(1)-Cg)*xx+g1(3)*g1(4)*exp(-g1(3)*g1(1)*xx));%C_u
%%%%%%%%%%%%%%%%%
Kh=1/Cds0;           %Gpa
Fc=(5.59-Vpg)/2.18;  
Vsg=3.52-1.89*Fc;     %km/s
den_g=Kg/(Vpg^2-4*Vsg^2/3);%kg/m^3
ug=den_g*Vsg^2; %Gpa
ug1=ug*10^9; 
sub=Kh/Kg;
syms phi_h 
phi_h=solve((1-phi_h)^(3/(1-phi_h))-sub);
solution_phi_h=double(phi_h); 
uh=ug*(1-solution_phi_h)^(3/(1-solution_phi_h)); %Gpa
alpha_c=Kh*(3*Kh+4*uh)/(pi*uh*theta_c)/(3*Kh+uh);

fc1=4*uh*(3*Kh+uh)*alpha_c^3/3/yita/(3*Kh+4*uh)*10^9;
%%%%%%%%%%%%%%%%%
K2=5;   %%%%%%%% 
p2=(Kg1+4*ug1/3)/(K2+4*ug1/3);
Ke1=Kg1+phi*(K2-Kg1)*p2;%pa
Ke=Ke1*10^(-9);
Ce=1/Ke;
% %%%%%%%%%%%%%%%%%
Cds_p=Cds0*(1-theta_s*(Cds0-Cg)*x);
Cds_p1=Cds0*(1-theta_s*(Cds0-Cg)*xx);
% %%%%%%%%%%%%%%%%%
f=@(B)Ce*(1+B(1)*B(2)*exp(-B(1)*Ce*x))-Cds_p;
d0=[1000 1200];
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
d1=lsqnonlin(f,d0,[],[],options);
yy1=Ce*(1+d1(1)*d1(2)*exp(-d1(1)*Ce*xx));
xxx = [4 0];
yyy=[0 -0.5];
zzz=[0 2.2];
theta_m=d1(1);
theta_m=0.1*theta_m;
phi_mo=d1(2);
% %%%%%%%%%%%%%%%
sub1=Ke/Kg;
syms phi_e 
phi_e=solve((1-phi_e)^(3/(1-phi_e))-sub1);
solution_phi_e=double(phi_e);
ue=ug*(1-solution_phi_e)^(3/(1-solution_phi_e));
alpha_m=Ke*(3*Ke+4*ue)/(pi*ue*theta_m)/(3*Ke+ue);
fc2=4*uh*(3*Kh+uh)*alpha_m^3/3/yita/(3*Kh+4*uh)*10^9;
%%%%%%%%%%%%%%%%%
Kds_p=1./Cds_p1;
Kds_p1=ones(length(omega),1)*Kds_p*10^9;%pa
yy2=1./yy;
Kd_p1=ones(length(omega),1)*yy2*10^9;
yy4=1./yy3;
ud_p1=ones(length(omega),1)*yy4*10^9;
phi_cp=phi_co*exp(-theta_c*Cds0*xx);%eq11
Sc=0.0011;
Kmf=1./((1./Kds_p1)+(1./((1./((1./Kd_p1)-(1./Kds_p1)))+(3*1i.*(omega'*yita*(Sc./1.5.*phi_cp.^-1))/8/alpha_c^2))));
umf=1./((1./ud_p1)-(4.*((1./Kd_p1)-(1./Kmf))/15));
Ks1=1./(1/Kg1+(phi*(1/Kf-1/Kg1)./(1+phi*(1/Kf-1/Kg1)./(1./Kmf-1/Kg1))));
us1=umf;
%%%%%%%%%%%%%%%%
phi_m=phi_mo*exp(-theta_m*Ce*xx);
Ke1=Ke*10^(8.83);
WR1=0.05;
Sw=0.99;
Sm=WR1*Sw*phi*phi_m.^-1;
Sm=Sm/min(Sm)+exp(-xx)-1;
Kmf1=1./((1/Ke1)+(1./(1./((1./Kds_p1)-(1/Ke1))+(3*1i.*(omega'*yita*(Sm./1.5.*phi_m.^-1))/8/alpha_m^2))));
umf1=1./((1./ud_p1)-(4.*((1./Kd_p1)-(1./Kmf1))/15));
Ks2=1./(1/Kg1+(phi*(1/Kf-1/Kg1)./(1+phi*(1/Kf-1/Kg1)./(1./Kmf1-1/Kg1))));
us2=umf1;
%%%%%%%%%%%%%%%%
Kmf2=1./(1./Kmf1+1./Kmf-1./Kds_p1);
umf2=1./((1./ud_p1)-(4.*((1./Kd_p1)-(1./Kmf2))/15));
%%%%%%%%%%%%%%%%
Ks=1./(1/Kg1+(phi*(1/Kf-1/Kg1)./(1+phi*(1/Kf-1/Kg1)./(1./Kmf2-1/Kg1))));
us=umf2;
Vp=sqrt((Ks+4*us/3)/den);
vp=real(Vp);
IQ_vp=imag(Vp)./real(Vp);
Vs=sqrt(us/den);
Es=9.*Ks.*us./(3*Ks+us);
Es_5=Es(:,1);   
Es_15=Es(:,10);  
Es_5_new = Es_5/1e9*0.85-9.1 +1.7;
Es_15_new= Es_15/1e9*0.95-5 +0.1;
IQ_E=imag(Es)./real(Es);
IQ_E_5=IQ_E(:,5);    
IQ_E_15=IQ_E(:,3);   
IQ_E_5_new = IQ_E_5./1.4*0.95+0.036;
IQ_E_15_new= IQ_E_15./2*0.95+0.024;
Ks_15=Ks(:,10);
Ks_5=Ks(:,2);
Ks_5_new = Ks_5./0.9e9-8.6;
Ks_15_new= Ks_15./0.7e9-12.3;
us_5=us(:,5);
us_15=us(:,8);
us_5_new = us_5./0.9e9*0.8-5.5 +1.8;
us_15_new= us_15./1e9*0.8-3.2 +1.95;
% =========================================================================
faic = 0.001;     alpha = alpha_c;  vis = 0.001;
Kdry_b = [7.5 10];
udry_b = [8.5 10];
for i = 1:length(Kdry_b)
 [Kmf_b(i,:),Gmf_b(i,:)] = K_G_mf(Kd_p(end).*1e9,Kdry_b(i).*1e9,udry_b(i).*1e9,omega,vis,faic,alpha,Kg.*1e9,4.25.*1e9,1,1);
 [Ksat(i,:),Gsat(i,:)] = Gassmann (Kmf_b(i,:),Gmf_b(i,:),Kg.*1e9,phi,4.25.*1e9);
 Ksat(i,:) = Ksat(i,:)./1e9-xxx(i);
 Gsat(i,:) = Gsat(i,:)./1e9+yyy(i);
 Esat(i,:) = 9*Ksat(i,:).*Gsat(i,:)./(3*Ksat(i,:) + Gsat(i,:))+zzz(i);
 Qe(i,:) = imag(Esat(i,:))./real(Esat(i,:));
 Vp_sq(i,:) = sqrt((Ksat(i,:)+4*Gsat(i,:)/3)/den);
 Vs_sq(i,:) = sqrt(Gsat(i,:)/den);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  E   %%%%%%%%%%%%%%%%%%%%%%%

frequencyrange1 = 10.^(-4:0.1:2);
af1 = 100;        
af2 = 40;
af = 70;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
open('./Data/e-E.fig');    hold on;    
    h1=semilogx(frequencyrange*4,Es_5_new,'-b','linewidth',1.2);   
    h2=semilogx(frequencyrange*4,Es_15_new,'-r','linewidth',1.2);
    h5=semilogx(frequencyrange(15:end)*1.8*af1,Esat(1,15:end),'-m','linewidth',1.2);    
    h6=semilogx(frequencyrange(19:end)*1.8*af2,Esat(2,19:end),'-g','linewidth',1.2);
    
    xlim([10^(-2) 10^6]);       ylim([10 35]);
    set(gca,'xtick',[0.01,1,10^2,10^4,10^6,10^7]);
    set(gca,'ytick',[10,15,20,25,30,35,40]);
    xlabel('Apparent frequency (Hz)');
    
    pp=[h1,h2,h5,h6];
    ah=axes('position',get(gca,'Position'),'visible','off');
    hh2=legend(ah,pp(1:4),'{\itP_{eff}}=5 MPa','{\itP_{eff}}=15 MPa','{\itP_{eff}}=5 MPa','{\itP_{eff}}=15 MPa', ...
        'FontName','Times New Roman','Fontsize',Fontsize);
    hh2.Title.String='Predictions SFTPS model    SF model';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  QE   %%%%%%%%%%%%%%%%%%%%%%%

open('./Data/f-QE.fig'); hold on;
    semilogx(frequencyrange*4.8,IQ_E_5_new,'-b','linewidth',1.2);   
    semilogx(frequencyrange*3.6,IQ_E_15_new,'-r','linewidth',1.2);    
    semilogx(frequencyrange*1.8*af1,Qe(1,:),'-m','linewidth',1.2);
    semilogx(frequencyrange*1.8*af2,Qe(2,:),'-g','linewidth',1.2);
    
    xlim([10^(-2) 10^6]);
     set(gca,'xtick',[0.01,1,10^2,10^4,10^6,10^7]);
    xlabel('Apparent frequency (Hz)');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  K   %%%%%%%%%%%%%%%%%%%%%%%

open('./Data/g-K.fig');hold on
    semilogx(frequencyrange*4,Ks_5_new,'-b','linewidth',1.2);
    semilogx(frequencyrange*4,Ks_15_new,'-r','linewidth',1.2);   

    semilogx(frequencyrange(13:end)*1.8*af1,Ksat(1,(13:end)),'-m','linewidth',1.2);
    semilogx(frequencyrange(17:end)*1.8*af2,Ksat(2,(17:end)),'-g','linewidth',1.2);
    
    xlim([10^(-2) 10^6]);    ylim([4 35]);    
    set(gca,'xtick',[0.01,1,10^2,10^4,10^6,10^7]);
    set(gca,'ytick',[5,10,15,20,25,30,35,40]);
    ylabel('Bulk modulus {\itK} (GPa)','Fontsize',13);
    xlabel('Apparent frequency (Hz)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% u   %%%%%%%%%%%%%%%%%%%%%%%
open('./Data/h-u.fig');   hold on;
    semilogx(frequencyrange*4,us_5_new,'-b','linewidth',1.2);
    semilogx(frequencyrange*4,us_15_new,'-r','linewidth',1.2);  
       
    semilogx(frequencyrange(13:end)*1.8*af1,Gsat(1,(13:end))-0.5,'-m','linewidth',1.2);
    semilogx(frequencyrange(17:end)*1.8*af2,Gsat(2,(17:end)),'-g','linewidth',1.2);
    
    xlim([10^(-2) 10^6]);     ylim([6 13]);
    set(gca,'xtick',[0.01,1,10^2,10^4,10^6,10^7]);
    ylabel('Shear modulus {\it\mu} (GPa)','Fontsize',13);
    xlabel('Apparent frequency (Hz)');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


