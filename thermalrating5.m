clear; close all; clc
%% dane dla interacji 
dt=1;                   % [s] krok
I(1:3600,1)=1300; % wartoœæ pr¹du
I(1:3600,2)=1200; % wartoœæ pr¹du
m=2;
NoS=length(I);          % liczba próbek 
t=(0:dt:length(I));     % wektor czasu
%% alokacja macierzy
rof=zeros(NoS,m);
kf=zeros(NoS,m);
mif=zeros(NoS,m);
R=zeros(NoS,m);
qc1=zeros(NoS,m);
qc2=zeros(NoS,m);
qcn=zeros(NoS,m);
qc=zeros(NoS,m);
qcmax=zeros(NoS,m);
qs=zeros(NoS,m);
qr=zeros(NoS,m);
Tfilm=zeros(NoS,m);
Tc=zeros(NoS,m);
dTc=zeros(NoS,m);
fi=zeros(NoS,m);
Imax=zeros(NoS,m);
Ta=zeros(NoS,m);
qc1max=zeros(NoS,m); % Alokacja macierzy
qc2max=zeros(NoS,m); % Alokacja macierzy
qcnmax=zeros(NoS,m); % Alokacja macierzy
qj=zeros(NoS,m);            % Alokacja macierzy
qrmax=zeros(NoS,m);        % Alokacja macierzy
qsmax=zeros(NoS,m);        % Alokacja macierzy
Rmax=zeros(NoS,m);    % Alokacja macier
Tfilmmax=zeros(NoS,m);% Alokacja macierzy
mifmax=zeros(NoS,m);       % Alokacja macierzy
rofmax=zeros(NoS,m); % Alokacja macierzy
kfmax=zeros(NoS,m);  % Alokacja macierzy
%% warunki pocz¹tkowe
Tc(1,1)=80;%('Podaj temperatuje pocz¹tkow¹ przewodu: ');
Tc(1,2)=Tc(1,1);
Tcmax=100;
%Ta(1:NoS,1)=40;
%Ta(1:NoS,2)=Ta(1:NoS,1);
% Vw(1:NoS,1)=0.61;
t1=1:NoS;
w=2*pi*0.04;
Ta(:,1)=40+0.5*(sin(w*t1));
fi=90;
% Vw(:,1)=abs(sin(w*t1));
Ta(1:1200,1)=23;
Ta(1201:2400,1)=30;
Ta(2401:3600,1)=42;
Ta(1:NoS,2)=Ta(1:NoS,1);
Vw(1:1200,1)=1.5;
Vw(1201:2400,1)=1.5;
Vw(2401:3600,1)=1.5;
Vw(1:NoS,2)=Vw(1:NoS,1);% prêdkoœæ wiatru
%fi(1:1200,1)=30;
%fi(1201:2400,1)=60;
%fi(2401:3600,1)=90;
%fi(1:NoS,2)=fi(1:NoS,1);%input('podaj k¹t padania wiatru wzglêdem przewodu: ');
%fi(1:NoS,2)=fi(1:NoS,1);
Kangle=1.194-cosd(fi)+0.194*cosd(2*fi)+0.368*sind(2*fi);%wind direction factor
%Qs=1023;    % promieniowanie s³oneczne
He=100; Lat=30; N=161; godzina=11; czystosc=1; alfa=0.5; 
Ctype=input('podaj typ przewodu "1" dla Drake26/7 ACSR "2" dla AFL-6 120mm2 "3" dla AFL-8 400mm2');
epsilon=0.5; delta=23.46*sind((284+N)*360/365); omega=(godzina-12)*15;
Hc=asind(cosd(Lat)*cosd(delta)+sind(Lat)*sind(delta)); As=1; Bs=1.148*10^-4; Cs=-1.108*10^-8;
if (czystosc==1)
A=-42.2391; B=63.8044; C=-1.9220; D=3.46921*10^-2; E=-3.61118*10^-4;
F=1.94318*10^-6; G=-4.07608*10^-9;
elseif (czystosc==0)
A=53.1821; B=14.2110; C=6.6138*10^-1; D=-3.1658*10^-2; E=5.4654*10^-4;
F=-4.3446*10^-6; G=1.3236*10^-8;
end
Qs=A+(B*Hc)+C*(Hc^2)+D*(Hc^3)+E*(Hc^4)+F*(Hc^5)+G*(Hc^6);
if (Ctype==1) %Drake ACSR
Srednica=28.1; mCp1=1.116*955; mCp2=0.5519*476; mCp=mCp1+mCp2;
elseif (Ctype==2) %AFL-6 120mm^2
Srednica=15.7; mCp1=0.326*955; mCp2=0.166*476; mCp=mCp1+mCp2;
elseif (Ctype==3) %AFL-8 400mm^2
Srednica=27.9; mCp1=1.084*955; mCp2=0.419*476; mCp=mCp1+mCp2; 
end
Ksolar=As+Bs*He+Cs*He^2; Qse=Ksolar*Qs;
X=sin(omega)/(sind(Lat)*cos(omega)-cos(Lat)*tan(delta));
Zc=139; szto=acosd(cosd(Hc)*cosd(Zc-90)); %78.62
Thigh=75; Tlow=25; RTlow=7.283*10^-5; RThigh=8.688*10^-5;
%% iteracja g³ówna
%Non Steady:
for k=1:m
for n=1:NoS  
Tfilm(n,k)=(Tc(n,k)+Ta(n))/2;
rof(n,k)=(1.293-1.525*10^-4*He+6.379*10^-9*He^2)/(1+0.00367*Tfilm(n,k));%1.029;
kf(n,k)=2.424*10^-2+7.477*10^-5*Tfilm(n,k)-4.407*10^-9*Tfilm(n,k)^2;  %0.0295;  % !! tu jest b³¹d 
mif(n,k)=1.458*10^-6*(Tfilm(n,k)+273)^1.5/(Tfilm(n,k)+383.4);  %2.04*10^-5;   % !! tu jest b³¹d 
% qc
qc1(n,k)=(1.01+0.0372*((Srednica*rof(n,k)*Vw(n,k))/mif(n,k))^0.52)*kf(n,k)*Kangle*(Tc(n,k)-Ta(n));%82,3
qc2(n,k)=(0.0119*((Srednica*rof(n,k)*Vw(n,k))/mif(n,k))^0.6)*kf(n,k)*Kangle*(Tc(n,k)-Ta(n));%76,9
qcn(n,k)=0.0205*rof(n,k)^0.5*Srednica^0.75*(Tc(n,k)-Ta(n))^1.25;
qc(n,k)=max([qc1(n,k),qc2(n,k),qcn(n,k)]);
% qr
qr(n,k)=0.0178*Srednica*epsilon*(((Tc(n,k)+273)/100)^4-((Ta(n)+273)/100)^4);%24,4
% qs
qs(n,k)=alfa*Qse*sind(szto)*Srednica/1000;
% R
R(n,k)=((RThigh-RTlow)/(Thigh-Tlow))*(Tc(n,k)-Tlow)+RTlow;
% dTc Tc
qj(n,k)=I(n,k)^2*R(n,k);
dTc(n,k)=((qs(n,k)+qj(n,k)-qc(n,k)-qr(n,k))/mCp)*dt;
Tc(n+1,k)=Tc(n,k)+dTc(n,k);
%Steady State 
Tfilmmax(n)=(Tcmax+Ta(n))/2;
mifmax(n)=1.458*10^-6*(Tfilmmax(n)+273)^1.5/(Tfilmmax(n)+383.4);
rofmax(n)=(1.293-1.525*10^-4*He+6.379*10^-9*He^2)/(1+0.00367*Tfilmmax(n));
kfmax(n)=2.424*10^-2+7.477*10^-5*Tfilmmax(n)-4.407*10^-9*Tfilmmax(n)^2;
qc1max(n)=(1.01+0.0372*((Srednica*rofmax(n)*Vw(n))/mifmax(n))^0.52)*kfmax(n)*Kangle*(Tcmax-Ta(n));
qc2max(n)=(0.0119*((Srednica*rofmax(n)*Vw(n))/mifmax(n))^0.6)*kfmax(n)*Kangle*(Tcmax-Ta(n));
qcnmax(n)= 0.0205*rofmax(n)^0.5*Srednica^0.75*(Tcmax-Ta(n))^1.25;
qcmax(n)=max([qc1max(n),qc2max(n),qcnmax(n)]);

qrmax(n)=0.0178*Srednica*epsilon*(((Tcmax+273)/100)^4-((Ta(n)+273)/100)^4);

qsmax(n)=alfa*Qse*sind(szto)*Srednica/1000; 

Rmax(n)=(RThigh-RTlow)/(Thigh-Tlow)*(Tcmax-Tlow)+RTlow;

Imax(n)=((qcmax(n)+qrmax(n)-qsmax(n))/Rmax(n))^0.5;
end
end
%%
plot((1:length(Tc))/60,Tc(:,1)); hold on; grid on
plot((1:length(Tc))/60,Tc(:,2)); 
title('nagrzewanie przewodu w czasie')
legend('1300 A','1200 A', 'Location', 'southeast')
xlabel('Time[s]');
ylabel('Temperature[*C]');
grid on
%%
figure
plot(dTc);
    plot(t(1:end-1),Imax(:,1)); hold on; grid on
    plot(t(1:end-1),Imax(:,2));
    title('30°')%zale¿noœæ wartoœci maksymalnej p³yn¹cego pr¹du od prêdkoœci wiatru
    xlabel('Time[s]');
    ylabel('Current[A]');
%%
% a_Drake_Imax=Imax;
% a_Drake_Tc=Tc;
% a_AFL6_Imax=Imax;
% a_AFL6_Tc=Tc;
% a_AFL8_Imax=Imax;
% a_AFL6_Tc=Tc;