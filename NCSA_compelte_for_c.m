clear all;clc
close all;
load raw_data.mat
%% 系统参数
C=1500;
fc=20000;
lambda=C./fc;
V=0.1;
%距离向参数
rB=10000;
H=30;
Rmax=300;
Rmin=7.5;

rFs=2.*rB;
rTs=1./rFs;
rs=rTs.*C./2;
Tp=10e-3;%脉冲宽度
rK=rB/Tp; 
rN_SAS=7680;
%rX2=linspace(Rmin,Rmax,rN_SAS);
rX=Rmin:rs:Rmin+rs*(rN_SAS-1);
rT=2*rX/C;
rF = ( -rN_SAS/2 : rN_SAS/2-1 )*( rFs/rN_SAS );% 距离频率向量
rc=(rX(rN_SAS/2+1)+rX(rN_SAS/2))./2;

%场景参数
M=10;                                       %子阵数量
DT=0.16;                                       %方位向天线长度
DI=0.08;                                     %子阵间距%
aX_min=-100;
aX_max=100;
r_width=10;

Lsar=lambda*rc/DT;%合成孔径
Tsar=Lsar/V;%照射时间
aK=2*V*V/lambda/rc;%教材4.38，斜视角为0
aB=aK*Tsar;
aFs=2.*aB;
aTs=1./aFs;
as=aTs.*V;

PRF=2*V/M/DI;
PRI=1/PRF;
% aN_SAS=ceil((aX_max-aX_min)/V/PRI);
% aX=linspace(aX_min,aX_max,aN_SAS);
% aT=aX/V;
% aF =( -aN_SAS/2 : aN_SAS/2-1 )*( PRF/aN_SAS );% 方位频率向量


%多子阵导致的升采样
aN_SAS=500;
PRF=PRF.*M;
PRI=1./PRF;
aN_SAS=aN_SAS*M;
aX=aX_min:as:aX_min+(aN_SAS-1)*as;
aF =( -aN_SAS/2 : aN_SAS/2-1 )*( PRF/aN_SAS );% 方位频率向量


%% 算法处理
% aimuth FFT
data_rd=fftshift(fft(raw_data,[],1),1);
% range FFT 
data_2df=fftshift(fft(data_rd,[],2),2);%行傅里叶变换

di=0;
r_ref=rc;
a = C^2-V^2;
%%  泰勒公式展开 距离方程导数
% ——————————————————————————参考距离——————————————————————————————————————
p=C^2-V^2;% 式子3-14
b1_ref=V*di+C*sqrt(r_ref.^2);% 式子3-14
b2_ref=di^2;

%式子 3-20
db1_ref=V^2;
ddb1_ref=C*V^2./sqrt(r_ref.^2);
dddb1_ref=0;
ddddb1_ref= -3*C*V^4./(r_ref.^2).^(3/2);

%式子 3-19
db2_ref = 2*V*di;
ddb2_ref = 0;
dddb2_ref = 0;
ddddb2_ref = 0;


dR0_r_ref=C./p.*(b1_ref+sqrt(b1_ref.^2+p.*b2_ref));% 式子3-16中R0,误差为0
dR1_r_ref=C./p.*(db1_ref+(2.*b1_ref.*db1_ref+p.*db2_ref)./(2.*sqrt(b1_ref.^2+p.*b2_ref)));% 式子3-21,误差小于1e-20
dR2_r_ref=C./p.*(... 
            ddb1_ref- (2.*b1_ref.*db1_ref+p.*db2_ref).^2./(4.*(b1_ref.^2+p.*b2_ref).^(3./2)) ...
            + (2.*db1_ref.^2+2.*b1_ref.*ddb1_ref+p.*ddb2_ref)./(2.*sqrt(b1_ref.^2+p.*b2_ref))...
           );% 式子3-22,误差小于1e-19
dR3_r_ref=C./p.*( dddb1_ref-...
     3.*(2.*b1_ref.*db1_ref+p.*db2_ref).*(2.*b1_ref.*ddb1_ref+2.*db1_ref.^2+p.*ddb2_ref)./(4.*(b1_ref.^2+p.*b2_ref).^(3./2))+...
     3.*(2.*b1_ref.*db1_ref+p.*db2_ref).^3./(8.*(b1_ref.^2+p.*b2_ref).^(5./2)) +...
     (6.*db1_ref.*ddb1_ref+2.*b1_ref.*dddb1_ref+p.*dddb2_ref)./(2.*(b1_ref.^2+p.*b2_ref).^(1./2))...
    );% 式子3-23,误差为0

dR4_r_ref=C./p.*(ddddb1_ref- ...
            15.*(2.*b1_ref.*db1_ref+p.*db2_ref).^4 ./(16.*(b1_ref.^2+p.*b2_ref).^(7./2))+...
            9.*(2.*b1_ref.*db1_ref+p.*db2_ref).^2.*(2.*b1_ref.*ddb1_ref+2.*db1_ref.^2+p.*ddb2_ref)./(4.*(b1_ref.^2+p.*b2_ref).^(5./2)) -...
            3.*(2.*b1_ref.*ddb1_ref+2.*db1_ref.^2+p.*ddb2_ref).^2./(4.*(b1_ref.^2+p.*b2_ref).^(3./2))-...
            (2.*b1_ref.*db1_ref+p.*db2_ref).*(6.*db1_ref.*ddb1_ref+2.*b1_ref.*dddb1_ref+p.*dddb2_ref)./((b1_ref.^2+p.*b2_ref).^(3./2))+...
            (6.*ddb1_ref.^2+8.*db1_ref.*dddb1_ref+2.*b1_ref.*ddddb1_ref+p.*ddddb2_ref)./(2.*(b1_ref.^2+p.*b2_ref).^(1./2))...
    );% 式子3-26,误差为小于1e-25

k0_r_ref=dR0_r_ref;
k1_r_ref=dR1_r_ref;
k2_r_ref=dR2_r_ref/2;
k3_r_ref=dR3_r_ref/6;
k4_r_ref=dR4_r_ref/24;

%——————————————————————————————————————距离——————————————————————————————————————
% 添加
b1=V*di+C*sqrt(rX.^2);% 式子3-14
b2=di^2;

%式子 3-20
db1=V^2;
ddb1=C*V^2./sqrt(rX.^2);
dddb1=0;
ddddb1= -3*C*V^4./(rX.^2).^(3/2);

%式子 3-19
db2 = 2*V*di;
ddb2 = 0;
dddb2 = 0;
ddddb2 = 0;

dR0=C./p.*(b1+sqrt(b1.^2+p.*b2));% 式子3-16中R0,误差小于1e-12
dR1=C./p.*(db1+(2.*b1.*db1+p.*db2)./(2.*sqrt(b1.^2+p.*b2)));% 式子3-21,误差小于1e-20
dR2=C./p.*(... 
            ddb1- (2.*b1.*db1+p.*db2).^2./(4.*(b1.^2+p.*b2).^(3./2)) ...
            + (2.*db1.^2+2.*b1.*ddb1+p.*ddb2)./(2.*sqrt(b1.^2+p.*b2))...
           );% 式子3-22,误差小于1e-18
dR3=C./p.*( dddb1-...
     3.*(2.*b1.*db1+p.*db2).*(2.*b1.*ddb1+2.*db1.^2+p.*ddb2)./(4.*(b1.^2+p.*b2).^(3./2))+...
     3.*(2.*b1.*db1+p.*db2).^3./(8.*(b1.^2+p.*b2).^(5./2)) +...
     (6.*db1.*ddb1+2.*b1.*dddb1+p.*dddb2)./(2.*(b1.^2+p.*b2).^(1./2))...
    );% 式子3-23,误差为0

dR4=C./p.*(ddddb1- ...
            15.*(2.*b1.*db1+p.*db2).^4 ./(16.*(b1.^2+p.*b2).^(7./2))+...
            9.*(2.*b1.*db1+p.*db2).^2.*(2.*b1.*ddb1+2.*db1.^2+p.*ddb2)./(4.*(b1.^2+p.*b2).^(5./2)) -...
            3.*(2.*b1.*ddb1+2.*db1.^2+p.*ddb2).^2./(4.*(b1.^2+p.*b2).^(3./2))-...
            (2.*b1.*db1+p.*db2).*(6.*db1.*ddb1+2.*b1.*dddb1+p.*dddb2)./((b1.^2+p.*b2).^(3./2))+...
            (6.*ddb1.^2+8.*db1.*dddb1+2.*b1.*ddddb1+p.*ddddb2)./(2.*(b1.^2+p.*b2).^(1./2))...
    );% 式子3-26,误差为小于1e-21

%泰勒级数展开系数
k0=dR0;
k1=dR1;
k2=dR2/2; 
k3=dR3/6;
k4=dR4/24;
%%  泰勒级数展开式 泰勒系数对距离的导数
%——————————————————————————————参考距离————————————————————————————————————————————
% 式子 3-38
 S1_ref=sqrt(C.^2.*r_ref.^2+2.*V.*di.*C.*r_ref+C.^2.*di.^2);
 S2_ref=(V.^2.*C.*r_ref+C.^2.*V.*di);
 S3_ref=(V.^4.*r_ref+di.*C.*V.^3+C.^2.*V.^2.*r_ref);
 S4_ref=(C.^2.*r_ref+V.*di.*C);
 %式子3--39
 dS1_ref=(C.^2.*r_ref+V.*di.*C)./S1_ref;
 dS2_ref=V.^2.*C;
 dS3_ref=V.^2.*C.^2+V.^4;
 
 
dk0_r_ref=C./p.*(C+C.*b1_ref./(S1_ref));%式子 3-28
dk1_r_ref=C.*V./p.*(C.*V./S1_ref-(C.*V.*r_ref+C.^2.*di).*C.*b1_ref./(S1_ref.^3));%式子 3-29



dk2_r_ref=C./p.*(-C.*V.^2./(r_ref.^2)+...
       3.*S2_ref.^2.*C.*b1_ref./(S1_ref.^5)-...
       2.*V.^2.*C.*S2_ref./(S1_ref.^3)-...
       di.*C.*V.^3./(r_ref.^2.*S1_ref)- ...
       S3_ref.*C.*b1_ref./(r_ref.*S1_ref.^3)...
    );% 式子3-30
dk2_r_ref=dk2_r_ref./2; %泰勒展开式需要处于n！

DF1_ref=-(dS2_ref.*S3_ref+S2_ref.*dS3_ref)./(r_ref.*S1_ref.^3) + S2_ref.*S3_ref.*(S1_ref.^3+3.*r_ref.*S1_ref.^2.*dS1_ref)./(r_ref.^2.*S1_ref.^6);%式子3-42
DF2_ref=3.*S2_ref.^2.*dS2_ref./S1_ref.^5 - 5.*S2_ref.^3.*dS1_ref./S1_ref.^6; %式子3-43
DF3_ref=-C.*V.^4.*(S1_ref+r_ref.*dS1_ref)./(r_ref.^2.*S1_ref.^2);%式子3-44
dk3_r_ref=3.*C./p.*(DF1_ref+DF2_ref+DF3_ref)./6;%泰勒展开式需要处于n！

DE1_ref=9.*C.*V.^4./r_ref.^4; %式子 3-48
DE2_ref=-15.*( 4.*S2_ref.^3.*dS2_ref./S1_ref.^7 - 7.*S2_ref.^4.*dS1_ref./S1_ref.^8 );%式子 3-49
DE3_ref=18.*((2.*S2_ref.*dS2_ref.*S3_ref+S2_ref.^2.*dS3_ref)./(r_ref.*S1_ref.^5) - S2_ref.^2.*S3_ref.*(S1_ref.^5+5.*r_ref.*S1_ref.^4.*dS1_ref)./(r_ref.^2.*S1_ref.^10));%式子3-50
DE4_ref=-3.*(2.*S3_ref.*dS3_ref./(r_ref.^2.*S1_ref.^3) - S3_ref.^2.*(2.*S1_ref+3.*r_ref.*dS1_ref)./(r_ref.^3.*S1_ref.^4));%式子3-51
DE5_ref=-12.*C.*V.^4.*(dS2_ref./(r_ref.*S1_ref.^3) - S2_ref.*(S1_ref.^3+3.*r_ref.*S1_ref.^2.*dS1_ref)./(r_ref.^2.*S1_ref.^6));%式子3-52
DE6_ref=-3.*C.*V.^5.*di.*(r_ref.^3.*dS1_ref+3.*r_ref.^2.*S1_ref)./(r_ref.^6.*S1_ref.^2);%式子3-53
dk4_r_ref=C./p.*(DE1_ref+DE2_ref+DE3_ref+DE4_ref+DE5_ref+DE6_ref)./(4*3*2*1);%泰勒展开式需要处于n！

%——————————————————————————————————————距离————————————————————————————————————————
% 
 S1=sqrt(C.^2.*rX.^2+2.*V.*di.*C.*rX+C.^2.*di.^2);
 S2=(V.^2.*C.*rX+C.^2.*V.*di);
 S3=(V.^4.*rX+di.*C.*V.^3+C.^2.*V.^2.*rX);
 S4=(C.^2.*rX+V.*di.*C);
 dS1=(C.^2.*rX+V.*di.*C)./S1;
 dS2=V.^2.*C;
 dS3=V.^2.*C.^2+V.^4;
 
dk0=C./p.*(C+C.*b1./(S1));%五个点误差为1e-16;其他为0

dk1=C.*V./p.*(C.*V./S1-(C.*V.*rX+C.^2.*di).*C.*b1./(S1.^3));%误差小于1E-21
dk2=C./p.*(-C.*V.^2./(rX.^2)+...
       3.*S2.^2.*C.*b1./(S1.^5)-...
       2.*V.^2.*C.*S2./(S1.^3)-...
       di.*C.*V.^3./(rX.^2.*S1)- ...
       S3.*C.*b1./(rX.*S1.^3)...
    );
dk2=dk2./2;%误差小于1E-18

DF1=-(dS2.*S3+S2.*dS3)./(rX.*S1.^3) + S2.*S3.*(S1.^3+3.*rX.*S1.^2.*dS1)./(rX.^2.*S1.^6);
DF2=3.*S2.^2.*dS2./S1.^5 - 5.*S2.^3.*dS1./S1.^6;
DF3=-C.*V.^4.*(S1+rX.*dS1)./(rX.^2.*S1.^2);
dk3=3.*C./p.*(DF1+DF2+DF3)./6;%误差小于1E-24

DE1=9.*C.*V.^4./rX.^4;
DE2=-15.*( 4.*S2.^3.*dS2./S1.^7 - 7.*S2.^4.*dS1./S1.^8 );
DE3=18.*((2.*S2.*dS2.*S3+S2.^2.*dS3)./(rX.*S1.^5) - S2.^2.*S3.*(S1.^5+5.*rX.*S1.^4.*dS1)./(rX.^2.*S1.^10));
DE4=-3.*(2.*S3.*dS3./(rX.^2.*S1.^3) - S3.^2.*(2.*S1+3.*rX.*dS1)./(rX.^3.*S1.^4));
DE5=-12.*C.*V.^4.*(dS2./(rX.*S1.^3) - S2.*(S1.^3+3.*rX.*S1.^2.*dS1)./(rX.^2.*S1.^6));
DE6=-3.*C.*V.^5.*di.*(rX.^3.*dS1+3.*rX.^2.*S1)./(rX.^6.*S1.^2);
dk4=C./p.*(DE1+DE2+DE3+DE4+DE5+DE6)./(4*3*2*1);%误差小于1E-22

%%  二维频域泰勒展开  补偿相位获取1 参考合成孔径声纳算法原理 公式3-58和3-59
A=aF'+k1.*fc./C;
B=k1./C;
%式子3-91
D1=       pi.*C./(2.*fc.*k2);
dD1=     -pi.*C./(2.*fc.^2.*k2);
ddD1=  2.*pi.*C./(2.*fc.^3.*k2);
dddD1=-6.*pi.*C./(2.*fc.^4.*k2);

D2=        pi.*k3.*C.^2./(4.*k2.^3.*fc.^2);
dD2=   -2.*pi.*k3.*C.^2./(4.*k2.^3.*fc.^3);
ddD2=   6.*pi.*k3.*C.^2./(4.*k2.^3.*fc.^4);
dddD2=-24.*pi.*k3.*C.^2./(4.*k2.^3.*fc.^5);

D3=      -pi.*C.^3.*(4.*k2.*k4-9.*k3.^2)./(32.*fc.^3.*k2.^5);
dD3=   3.*pi.*C.^3.*(4.*k2.*k4-9.*k3.^2)./(32.*fc.^4.*k2.^5);
ddD3=-12.*pi.*C.^3.*(4.*k2.*k4-9.*k3.^2)./(32.*fc.^5.*k2.^5);
dddD3=60.*pi.*C.^3.*(4.*k2.*k4-9.*k3.^2)./(32.*fc.^6.*k2.^5);

%泰勒级数展开系数，即式子15 中的
%fa_az_i,fa_rcmc_i,fa_rc_i,fa_src_i,分别对应phase_0、phase_1、phase_2、phase_3
%此处验证与推导的公式一样，因为是泰勒展开所以要处于1/(n!)
%式子3-92
phase_0=-2.*pi.*fc.*k0./C+D1.*A.^2+D2.*A.^3+D3.*A.^4; %误差小于1e-11
phase_1=-2.*pi.*k0./C + (dD1.*A.^2+dD2.*A.^3+dD3.*A.^4)+(2*D1.*B.*A + 3.*B.*A.^2.*D2 + 4.*B.*A.^3.*D3); %误差小于1e-15
phase_2=(ddD1.*A.^2+ddD2.*A.^3+ddD3.*A.^4) + 2.*(2.*dD1.*B.*A+3.*B.*A.^2.*dD2+4.*B.*A.^3.*dD3)+2.*(B.^2.*D1+3.*B.^2.*A.*D2+6.*B.^2.*A.^2.*D3);
phase_2=phase_2./2;%误差小于1e-19
phase_3=(dddD1.*A.^2+dddD2.*A.^3+dddD3.*A.^4)+3.*(2.*ddD1.*B.*A+3.*B.*A.^2.*ddD2+4.*B.*A.^3.*ddD3)+...
          6.*(B.^2.*dD1+3.*B.^2.*A.*dD2+6.*B.^2.*A.^2.*dD3)+6.*(B.^3.*D2+4.*B.^3.*A.*D3);
phase_3=phase_3./6;%误差小于1e11



%%  补偿相位获取 2 参考合成孔径声纳算法原理 公式3-58和3-59
A_ref=aF'+k1_r_ref.*fc./C;
B_ref=k1_r_ref./C;
%式 3-91
D1_ref=       pi.*C./(2.*fc.*k2_r_ref);
dD1_ref=     -pi.*C./(2.*fc.^2.*k2_r_ref);
ddD1_ref=  2.*pi.*C./(2.*fc.^3.*k2_r_ref);
dddD1_ref=-6.*pi.*C./(2.*fc.^4.*k2_r_ref);

D2_ref=        pi.*k3_r_ref.*C.^2./(4.*k2_r_ref.^3.*fc.^2);
dD2_ref=   -2.*pi.*k3_r_ref.*C.^2./(4.*k2_r_ref.^3.*fc.^3);
ddD2_ref=   6.*pi.*k3_r_ref.*C.^2./(4.*k2_r_ref.^3.*fc.^4);
dddD2_ref=-24.*pi.*k3_r_ref.*C.^2./(4.*k2_r_ref.^3.*fc.^5);

D3_ref=      -pi.*C.^3.*(4.*k2_r_ref.*k4_r_ref-9.*k3_r_ref.^2)./(32.*fc.^3.*k2_r_ref.^5);
dD3_ref=   3.*pi.*C.^3.*(4.*k2_r_ref.*k4_r_ref-9.*k3_r_ref.^2)./(32.*fc.^4.*k2_r_ref.^5);
ddD3_ref=-12.*pi.*C.^3.*(4.*k2_r_ref.*k4_r_ref-9.*k3_r_ref.^2)./(32.*fc.^5.*k2_r_ref.^5);
dddD3_ref=60.*pi.*C.^3.*(4.*k2_r_ref.*k4_r_ref-9.*k3_r_ref.^2)./(32.*fc.^6.*k2_r_ref.^5);




%式子3-92
%1阶导
phase_1_r_ref=-2.*pi.*k0_r_ref./C + (dD1_ref.*A_ref.^2+dD2_ref.*A_ref.^3+dD3_ref.*A_ref.^4)+(2*D1_ref.*B_ref.*A_ref + 3.*B_ref.*A_ref.^2.*D2_ref + 4.*B_ref.*A_ref.^3.*D3_ref);%误差小于1e-15
phase_2_r_ref=(ddD1_ref.*A_ref.^2+ddD2_ref.*A_ref.^3+ddD3_ref.*A_ref.^4) ...
               + 2.*(2.*dD1_ref.*B_ref.*A_ref+3.*B_ref.*A_ref.^2.*dD2_ref ... 
               +4.*B_ref.*A_ref.^3.*dD3_ref)+2.*(B_ref.^2.*D1_ref+3.*B_ref.^2.*A_ref.*D2_ref+6.*B_ref.^2.*A_ref.^2.*D3_ref);
phase_2_r_ref=(-2.*pi./rK+phase_2_r_ref)./2;%误差小于1e-20 二阶导
phase_3_r_ref=(dddD1_ref.*A_ref.^2+dddD2_ref.*A_ref.^3+dddD3_ref.*A_ref.^4)+3.*(2.*ddD1_ref.*B_ref.*A_ref+3.*B_ref.*A_ref.^2.*ddD2_ref+4.*B_ref.*A_ref.^3.*ddD3_ref)+...
               6.*(B_ref.^2.*dD1_ref+3.*B_ref.^2.*A_ref.*dD2_ref+6.*B_ref.^2.*A_ref.^2.*dD3_ref)+6.*(B_ref.^3.*D2_ref+4.*B_ref.^3.*A_ref.*D3_ref);
phase_3_r_ref=phase_3_r_ref./6;           
% phase_1_r_ref 对应 式16 中fa_rcmc_i
% phase_1_dr_r_ref  对应 式16 中fa_rcmc_i ,r换成r_ref
% phase_2_r_ref  对应式20分母，即fa_src_i+fa_rc
% phase_2_dr_r_ref 
%phase_3_r_ref 对应facubic,i(fa; rref )

%----------------------------------------------------------------------对距离求导--------------------------------------------------------------
%式3-110
d_A=dk1_r_ref.*fc./C;
d_A2=2.*A_ref.*d_A;
d_A3=3.*A_ref.^2.*d_A;
d_A4=4.*A_ref.^3.*d_A;
d_B=dk1_r_ref./C;
d_B2=2.*B_ref.*d_B;
%式3-111
d_B1A1=B_ref.*d_A+A_ref.*d_B;
d_B1A2=B_ref.*d_A2+A_ref.^2.*d_B;
d_B1A3=B_ref.*d_A3+A_ref.^3.*d_B;
d_B2A1=B_ref.^2.*d_A+A_ref.*d_B2;
d_B2A2=B_ref.^2.*d_A2+A_ref.^2.*d_B2;
%式3-112
d_ddD1=-pi.*(C./fc.^3).*dk2_r_ref./(k2_r_ref.^2);
d_dD1=  pi.*(C./(2.*fc.^2)).*dk2_r_ref./(k2_r_ref.^2);
d_D1=  -pi.*(C./(2.*fc)).*dk2_r_ref./(k2_r_ref.^2);
%式3-113
temp2=dk3_r_ref./(k2_r_ref.^3)-3.*k3_r_ref.*dk2_r_ref./(k2_r_ref.^4);
d_ddD2=6.*pi.*(C.^2./(4.*fc.^4)).*temp2;
d_dD2=   -pi.*(C.^2./(2.*fc.^3)).*temp2;
d_D2 =    pi.*(C.^2./(4.*fc.^2)).*temp2;
%式3-113
temp3=-(4.*dk2_r_ref.*k4_r_ref-18.*dk3_r_ref.*k3_r_ref+4.*dk4_r_ref.*k2_r_ref)./(32.*k2_r_ref.^5) ...
         -5.*(9.*k3_r_ref.^2-4.*k2_r_ref.*k4_r_ref).*dk2_r_ref./(32.*k2_r_ref.^6);
d_ddD3=12.*pi.*(C.^3./(fc.^5)).*temp3;
d_dD3= -3.*pi.*(C.^3./(fc.^4)).*temp3;
d_D3=      pi.*(C.^3./(fc.^4)).*temp3;     

%式子3-110，误差小于1E-7
phase_1_dr_r_ref=-2.*pi.*dk0_r_ref./C+ ...
                  d_dD1.*A_ref.^2+dD1_ref.*d_A2+d_dD2.*A_ref.^3+dD2_ref.*d_A3+ d_dD3.*A_ref.^4+dD3_ref.*d_A4+ ...
                  2.*D1_ref.*d_B1A1+2.*d_D1.*B_ref.*A_ref+ 3.*D2_ref.*d_B1A2+3.*d_D2.*B_ref.*A_ref.^2+ 4.*D3_ref.*d_B1A3+4.*d_D3.*B_ref.*A_ref.^3;

%式子3-110，误差小于1E-15
phase_2_dr_r_ref= d_ddD1.*A_ref.^2+ddD1_ref.*d_A2+ d_ddD2.*A_ref.^3+ ddD2_ref.*d_A3+ d_ddD3.*A_ref.^4+ ddD3_ref.*d_A4 + ....
                   2.*(2.*d_dD1.*B_ref.*A_ref+ 2.*dD1_ref.*d_B1A1+ 3.*dD2_ref.*d_B1A2+3.*d_dD2.*B_ref.*A_ref.^2+ 4.*dD3_ref.*d_B1A3+ 4.*d_dD3.*B_ref.*A_ref.^3 )+ ...
                   2*(d_B2.*D1_ref+ B_ref.^2.*d_D1+ 3.*d_B2A1.*D2_ref+ 3.*B_ref.^2.*A_ref.*d_D2+ 6.*d_B2A2.*D3_ref+ 6.*B_ref.^2.*A_ref.^2.*d_D3);
phase_2_dr_r_ref=phase_2_dr_r_ref./2;


%% 数据处理
% NCS para
Km = pi./phase_2;
Km_ref =pi./(phase_2_r_ref);% 

gama=-4*pi./C./phase_1_dr_r_ref;  % phase_1_dr_r_ref %y_i_fa,等于 2/R_rd_i
%     gama=sqrt(1-C^2*aF.^2/(4*V^2*fc^2));

gama_fa_ref=gama(aN_SAS./2);
alpha = gama_fa_ref./gama+eps; %

Ks=-pi.*phase_2_dr_r_ref./(phase_2_r_ref).^2;% KS_i_ref 式子20
Ks=Ks.*C.*gama/2;

R_rcm = phase_1./(-2*pi)*C; 
R_rcm_r_ref = phase_1_r_ref./(-2*pi)*C;%式16

q2 = Km_ref.*(alpha-1); %式子29
q3 = (1/2).*Ks.*(alpha-1); %式子30
Ym = Ks.*(2.*alpha-1)./(2.*Km_ref.^3.*(alpha-1));%式子28
Y = Ym-3./(2.*pi).*phase_3_r_ref;%式子22

% high order phase compensation
H_cubic=+2*pi/3*Y.*rF.^3;%式子21
data_filter_3_order=data_2df.*exp(+1j.*H_cubic);

% range IFFT 
data_rd_1=ifft(fftshift(data_filter_3_order,2),[],2);

% NCS 
phase_scaling=-pi.*q2.*(rT-2*r_ref/C./gama).^2-2.*pi./3.*q3.*(rT-2*r_ref/C./gama).^3; 
data_scaled=data_rd_1.*exp(1j*phase_scaling);

% range FFT 
data_2df_1=fftshift(fft(data_scaled,[],2),2);

% filtering
Href_1= -pi.*rF.^3.*Ks./(3.*alpha.*(alpha-1).*Km_ref.^3)-pi.*rF.^2./(Km_ref.*alpha)+2.*pi.*(alpha-1).*R_rcm_r_ref.*rF./(alpha.*C);
data_Href_1=data_2df_1.*exp(+1j*Href_1);

% range IFFT 
data_rd_2=ifft(fftshift(data_Href_1,2),[],2);
% azimuth processing
Href_2= pi.*Km_ref.*(alpha-1).*(2.*rX./C./gama-2.*r_ref./C./gama).^2./alpha+Ks.*(alpha-1).*(2.*rX./C./gama-2.*r_ref./C./gama).^3.*pi./(3.*alpha);
data_rd_3=data_rd_2.*exp(+1j*Href_2);
Href_3= -phase_0;
data_ac=data_rd_3.*exp(+1j*Href_3);
% azimuth IFFT 
resultData=ifft(ifftshift(data_ac,1),[],1);
M=max(max(abs(resultData)));

%toc;
% figures();
figure; imagesc(rX,aX,20*log10(abs(resultData)/M),[-40,0])
xlabel('range(m)','Fontsize',20);ylabel('azimuth(m)','Fontsize',20);
%colormap(1-gray);grid;
%% 目标单列
figure;subplot(1,3,1);
imagesc(rX,aX,20*log10(abs(resultData)/M),[-40,0])
xlim([35 65]);
xlabel('range(m)','Fontsize',20);ylabel('azimuth(m)','Fontsize',20);

subplot(1,3,2);
imagesc(rX,aX,20*log10(abs(resultData)/M),[-40,0])
xlim([135 165]);
xlabel('range(m)','Fontsize',20);ylabel('azimuth(m)','Fontsize',20);

subplot(1,3,3);
imagesc(rX,aX,20*log10(abs(resultData)/M),[-40,0])
xlim([235 265]);
xlabel('range(m)','Fontsize',20);ylabel('azimuth(m)','Fontsize',20);
%% 插值前的距离剖面
 Max_one=max(max(abs(resultData)));  %最大值
[i,j]=find(abs(resultData)==Max_one);%最大值所在行列
peak_x=rX(j);peak_y=aX(i);
fprintf('峰值点为[%d,%d]', peak_x,peak_y);
R33=abs(resultData(:,j));%抽取列,距离向
A33=abs(resultData(i,:));%抽取行，方位向
r2=20*log10(R33/max(R33));%%距离向归一化
a2=20*log10(A33/max(A33));%%方位向归一化
figure;subplot(2,2,1);plot(rX,a2);title('方位向切面');
xlabel('range(m)','Fontsize',10);ylabel('幅值-归一化','Fontsize',10);
subplot(2,2,2);plot(aX,r2);title('距离向切面');
xlabel('azimuth(m)','Fontsize',10);ylabel('幅值-归一化','Fontsize',10);
%% 计算ISLR和PSLR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%二维插值处理，D2为二维插值倍数%%%%%%%%%%%%%%%%%%
echo=resultData';%%每列有512个点，列为距离向
Max=max(max(abs(echo)));  %最大值
[X,Y]=find(abs(echo)==Max);%峰值坐标
DArea=echo(X-128:X+128,Y-128:Y+128); %%取257*257区域

D2=8;
A6= islr_pslr_interpret2d(DArea,D2);
rX_2D=rX(X-128:X+128);
aX_2D=aX(Y-128:Y+128);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%一维维插值处理，D1为插值倍数%%%%%%%%%%%%%%%%%%
 Max_one=max(max(abs(A6)));  %最大值
[i,j]=find(abs(A6)==Max_one);%最大值所在行列
A=A6(i-128:i+128,j);%抽取列,距离向
B=A6(i,j-128:j+128);%抽取行，方位向
D1=128;%插值倍数
R2=islr_pslr_interpret1d(A.',D1);
A2=islr_pslr_interpret1d(B,D1);
r2=20*log10(R2/max(R2));%%距离向归一化
a2=20*log10(A2/max(A2));%%方位向归一化



%  Max_one=max(max(abs(resultData)));  %最大值
% [i,j]=find(abs(resultData)==Max_one);%最大值所在行列
% R33=abs(resultData(:,j));%抽取列,距离向
% A33=abs(resultData(i,:));%抽取行，方位向
% r2=20*log10(R33/max(R33));%%距离向归一化
% a2=20*log10(A33/max(A33));%%方位向归一化
% figure;plot(rX,a2);title('方位向切面');
% figure;plot(aX,r2);title('距离向切面');

subplot(2,2,3);plot(a2);title('方位向切面(插值补零后)');
xlabel('range(m)','Fontsize',10);ylabel('幅值-归一化(dB)','Fontsize',10);
subplot(2,2,4);plot(r2);title('距离向切面(插值补零后)');
xlabel('azimuth(m)','Fontsize',10);ylabel('幅值-归一化(dB)','Fontsize',10);
[ISLRa,ISLRr,PSLRa,PSLRr,dIRWa,dIRWr,] = get_islr_pslr_2d(A2,R2);
dt=rT(2)-rT(1);
fs=1/dt;
IRWa=V*dIRWa/(PRF*D1*D2);
IRWr=C*dIRWr/(2*fs*D1*D2);