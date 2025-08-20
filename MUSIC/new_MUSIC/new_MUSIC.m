%% 
clear all;
close all;
clc;
%%
c = 1500;
M = 20;                             %阵元数
f0 = 500;
N_x = 4096;                         %信号长度
T = 1;
t = 0:1/N_x:T-1/N_x;                %信号时间
fs = N_x;                           %快拍数
l = c / f0;                         %信号波长  
d = 0.5*l;                          %阵元间距
snr = 10;                           %信噪比
Theta=[35 85];                      %两个信号的入射角度
source_number = length(Theta);      %信元数
f = fs/N_x*(0:N_x-1);               %具体频点
fl = f0; fh = f0;

%% 生成信号
rng('default');
sig1 = 2*cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(1)) / c));   % 目标信号
sig2 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(2)) / c));   % 干扰信号
noise = randn(M,length(sig1));                              % 均值为0，方差为1的高斯白噪声
SNR1 = 10^(10 / 10);   SNR2 = 10^(-30 / 10);
P1 = mean(sig1.^2);    P2 = mean(sig2.^2);                  % 计算信号功率
PN1 = P1 / SNR1;       PN2 = P2 / SNR2;                     % 计算所需的噪声功率
noise_scale1 = sqrt(PN1 / mean(noise.^2));                  % 计算噪声缩放因子
noise_scale2 = sqrt(PN2 / mean(noise.^2));
sig1 = sig1 + noise_scale1*noise;
sig2 = sig2 + noise_scale2*noise;
Sig = sig1 + sig2;

%% DOA
theta = linspace(0,180,360);
x = fft(Sig').';
I = flip(eye(M));
for i = 1 : length(theta)
    for j = find(f>=fl & f<=fh)
        R = x(:,j) * x(:,j)' / N_x;
        atheta = exp(-1j*(0:M-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
%         R = R + (I*conj(R)*I);
        [U,V] = eig(R);                
        Un = U(:,1:end-2);             % 噪声子空间
        En = U(:,end-1:end);           % 信号子空间
        Vn = Un*Un';
        Vs = En*En';
        
        Ps(i,j) = 1 / (atheta'*Vn*atheta);
        Pss(i,j) = (atheta'*Vs*atheta) / (atheta'*Vn*atheta);
    end
end
P_music = sum(Ps.' .* conj(Ps.')) / length(find(f>=fl & f<=fh));
P_music_new = sum(Pss.' .* conj(Pss.')) / length(find(f>=fl & f<=fh));


%% 画图
figure;
plot(theta,10*log10(P_music/max(P_music)));
xlabel('角度/°');ylabel('幅度/dB');
hold on;
plot(theta,10*log10(P_music_new/max(P_music_new)));
legend('MUSIC','改进MUSIC');

for i = 1:source_number
    xline(Theta(i),'--');
end

