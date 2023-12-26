clear;
clc;
close all;

Pref=20*10^(-6);
% Pref: Reference Pressure, Air=20 /mu Pa, Water= 1 /mu Pa


% Detect Data plot(time steps,kPa)
b1=readmatrix("Wall_Pressure_SM.txt");
low = 4;                % low : Start data(row)
high1 = length(b1);       % high : End data(row)

% data variable
t1=b1(low:high1,3);
press_P1_S=b1(low:high1,5);
pressfluc_P1_S = press_P1_S-mean(press_P1_S);

fs=100000;
dt=1/fs;
N1=length(press_P1_S);
T1=dt*N1;
t_s1=0:dt:T1;
df1 = fs/N1;
ff1=[0:N1-1]'*(1/T1);

% FFT
w1=hann(N1);
press_P1_HA=w1.*pressfluc_P1_S;
press_P1_FT=fft(press_P1_HA);
abs_press_P1=2.*abs(press_P1_FT)./N1; %정규화
P_press_P1= 2.*abs_press_P1;

SPL_P1=20*log10(P_press_P1./Pref);
PSD_P1=SPL_P1-10*log10(df1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Detect Data plot(time steps,kPa)
b2=readmatrix("Wall_Pressure_WALE.txt");
low = 4;                % low : Start data(row)
high2 = length(b2);       % high : End data(row)

% data variable
t2=b2(low:high2,3);
press_P2_S=b2(low:high2,5);
pressfluc_P2_S = press_P2_S-mean(press_P2_S);

fs2=100000;
dt2=1/fs;
N2=length(press_P2_S);
T2=dt*N2;
t_s2=0:dt:T2;
df2 = fs/N2;
ff2=[0:N2-1]'*(1/T2);

% FFT
w2=hann(N2);
press_P2_HA=w2.*pressfluc_P2_S;
press_P2_FT=fft(press_P2_HA);
abs_press_P2=2.*abs(press_P2_FT)./N2; %정규화
P_press_P2= 2.*abs_press_P2;

SPL_P2=20*log10(P_press_P2./Pref);
PSD_P2=SPL_P2-10*log10(df2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detect Data plot(time steps,kPa)
b3=readmatrix("Wall_Pressure_WMLES.txt");
low = 4;                % low : Start data(row)
high3 = length(b3);       % high : End data(row)

% data variable
t3=b3(low:high3,3);
press_P3_S=b3(low:high3,5);
pressfluc_P3_S = press_P3_S-mean(press_P3_S);

fs3=100000;
dt3=1/fs;
N3=length(press_P3_S);
T3=dt*N3;
t_s3=0:dt:T3;
df3 = fs/N3;
ff3=[0:N3-1]'*(1/T3);

% FFT
w3=hann(N3);
press_P3_HA=w3.*pressfluc_P3_S;
press_P3_FT=fft(press_P3_HA);
abs_press_P3=2.*abs(press_P3_FT)./N3; %정규화
P_press_P3= 2.*abs_press_P3;

SPL_P3=20*log10(P_press_P3./Pref);
PSD_P3=SPL_P3-10*log10(df3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detect Data plot(time steps,kPa) - S_Omega
b4=readmatrix("Wall_Pressure_S_Omega.txt");
low = 4;                % low : Start data(row)
high4 = length(b4);       % high : End data(row)

% data variable
t4=b4(low:high4,3);
press_P4_S=b4(low:high4,5);
pressfluc_P4_S = press_P4_S-mean(press_P4_S);

fs4=100000;
dt4=1/fs;
N4=length(press_P4_S);
T4=dt*N4;
t_s4=0:dt:T4;
df4 = fs/N4;
ff4=[0:N4-1]'*(1/T4);

% FFT
w4=hann(N4);
press_P4_HA=w4.*pressfluc_P4_S;
press_P4_FT=fft(press_P4_HA);
abs_press_P4=2.*abs(press_P4_FT)./N4; %정규화
P_press_P4= 2.*abs_press_P4;

SPL_P4=20*log10(P_press_P4./Pref);
PSD_P4=SPL_P4-10*log10(df4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detect Data plot(time steps,kPa) - S_Omega
b5=readmatrix("Wall_Pressure_KET.txt");
low = 4;                % low : Start data(row)
high5 = length(b5);       % high : End data(row)

% data variable
t5=b5(low:high5,3);
press_P5_S=b5(low:high5,5);
pressfluc_P5_S = press_P5_S-mean(press_P5_S);

fs5=100000;
dt5=1/fs;
N5=length(press_P5_S);
T5=dt*N5;
t_s5=0:dt:T5;
df5 = fs/N5;
ff5=[0:N5-1]'*(1/T5);

% FFT
w5=hann(N5);
press_P5_HA=w5.*pressfluc_P5_S;
press_P5_FT=fft(press_P5_HA);
abs_press_P5=2.*abs(press_P5_FT)./N5; %정규화
P_press_P5= 2.*abs_press_P5;

SPL_P5=20*log10(P_press_P5./Pref);
PSD_P5=SPL_P5-10*log10(df5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detect Data plot(time steps,kPa) - S_Omega
b6=readmatrix("exp_K3.xlsx");



figure(1)   
hold on;
plot(ff1,PSD_P1,'b','linewidth',1.0);
plot(ff2,PSD_P2,'r','linewidth',1.0);
plot(ff3,PSD_P3,'g','linewidth',1.0);
plot(ff4,PSD_P4,'k','linewidth',1.0);
plot(ff5,PSD_P5,'y','linewidth',1.0);
plot(b6(:,1),b6(:,2),'linewidth',4.0);

hold off;

set(gca,'xscale','log');
xlim([50 5000])
ylim([50,130])
grid on;
title("PSD-Harforice");
xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [(dB/Hz^1^/^2], P_r_e_f=20 \muPa')
legend('K4-SM model', 'K4-WALE model','K4-WMLES model', 'K4-S-OMEGA model', 'K4-KET model','Experiment')

figure(2)

hold on;
plot(b6(:,1),b6(:,2),'r','linewidth',2.0);
plot(ff1,PSD_P1,'b','linewidth',1.0);
hold off;

%set(gca,'xscale','log');
xlim([50 5000])
ylim([50,130])
grid on;
title("PSD-Harforice-SM");
xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [(dB/Hz^1^/^2], P_r_e_f=20 \muPa')
legend('Experiment', 'K4-SM model')

figure(3)

hold on;
plot(b6(:,1),b6(:,2),'r','linewidth',2.0);
plot(ff2,PSD_P2,'b','linewidth',1.0);
hold off;

%set(gca,'xscale','log');
xlim([50 5000])
ylim([50,130])
grid on;
title("PSD-Harforice-WALE");
xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [(dB/Hz^1^/^2], P_r_e_f=20 \muPa')
legend('Experiment', 'WALE model')

figure(4)

hold on;
plot(b6(:,1),b6(:,2),'r','linewidth',2.0);
plot(ff3,PSD_P3,'b','linewidth',1.0);
hold off;

%set(gca,'xscale','log');
xlim([50 5000])
ylim([50,130])
grid on;
title("PSD-Harforice-WMLES");
xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [(dB/Hz^1^/^2], P_r_e_f=20 \muPa')
legend('Experiment', 'WMLES model')

figure(5)

hold on;
plot(b6(:,1),b6(:,2),'r','linewidth',2.0);
plot(ff4,PSD_P4,'b','linewidth',1.0);
hold off;

%set(gca,'xscale','log');
xlim([50 5000])
ylim([50,130])
grid on;
title("PSD-Harforice-S-OMEGA");
xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [(dB/Hz^1^/^2], P_r_e_f=20 \muPa')
legend('Experiment', 'S-OMEGA model')

figure(6)

hold on;
plot(b6(:,1),b6(:,2),'r','linewidth',2.0);
plot(ff5,PSD_P5,'b','linewidth',1.0);
hold off;

%set(gca,'xscale','log');
xlim([50 5000])
ylim([50,130])
grid on;
title("PSD-Harforice-KET");
xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [(dB/Hz^1^/^2], P_r_e_f=20 \muPa')
legend('Experiment', 'KET model')

