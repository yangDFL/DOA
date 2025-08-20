%% 矩阵分解 MD算法
close all; clc; clear;
%%
c = 1500;
M = 20;                             %阵元数
m = 12;
p = M-m+1;
f0 = 500;
N_x = 4096;                         %信号长度
T = 1;
t = 0:1/N_x:T-1/N_x;                %信号时间
fs = N_x;                           %快拍数
l = c / f0;                         %信号波长  
d = 0.5*l;                          %阵元间距
Theta=[20 23 100];                  %两个信号的入射角度
source_number = length(Theta);      %信元数
f = fs/N_x*(0:N_x-1);               %具体频点
fl = f0; fh = f0;
%% 生成信号
rng('default');
sig1 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(1)) / c));   % 目标信号
sig2 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(2)) / c));   
sig3 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(3)) / c)); 
noise = randn(M,length(sig1));                              % 均值为0，方差为1的高斯白噪声
snr = 10;
SNR1 = 10^(snr / 10);  SNR2 = 10^(snr / 10); SNR3 = 10^(snr/10);
P1 = mean(sig1.^2);  P2 = mean(sig2.^2); P3 = mean(sig3.^2);                 % 计算信号功率
PN1 = P1 / SNR1;       PN2 = P2 / SNR2;  PN3 = P3 / SNR3;                  % 计算所需的噪声功率
noise_scale1 = sqrt(PN1 / mean(noise.^2));                  % 计算噪声缩放因子
noise_scale2 = sqrt(PN2 / mean(noise.^2));
noise_scale3 = sqrt(PN3 / mean(noise.^2));
sig1 = sig1 + noise_scale1*noise;
sig2 = sig2 + noise_scale2*noise;
sig3 = sig3 + noise_scale3*noise;
Sig = sig1 + sig2 + sig3;
% Sig = sig1 + sig2;
%% DOA估计
theta = linspace(0,180,180);
x = fft(Sig').';
N = 14;                         % 平滑后子阵数量
sub_M = M-N+1;                  % 子阵阵元数
for i = 1 : length(theta)
    for j = find(f>=fl & f<=fh)
        R = x(:,j) * x(:,j)' / N_x;
        atheta = exp(-1j*(0:M-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
        atheta1 = exp(-1j*(0:m-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
        atheta2 = exp(-1j*(0:N-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
        %% MMD
        l0 = M-m-1;
        J1 = fliplr(eye(M)); J2 = fliplr(eye(m));
        for ii = 0:l0
            Rm1(:,ii*M+1:(ii+1)*M) = R(ii+1:ii+m,:);
            Rm1(:,(ii+l0+1)*M+1:(ii+l0+2)*M) = J2 * conj(R(ii+1:ii+m,:)) * J1;
        end
        [U,S,V] = svd(Rm1);
        Un1 = U(:,source_number+1:end);
        %% SMD
        Rx1 = zeros(N,N);
        for ii = 1:sub_M   
            Rx1 = Rx1 + R(ii:ii+N-1, ii:ii+N-1);
        end
        R1 = Rx1 / sub_M;

        RR = J1*conj(R)*J1;
        Rx2 = zeros(N,N);
        for ii = 1:sub_M
            Rx2 = Rx2 + RR(ii:ii+N-1, ii:ii+N-1);
        end
        R2 = Rx2 / sub_M;
        Rm2 = [R1 R2];
        [U,S,V] = svd(Rm2);
        Un2 = U(:,source_number+1:end);

        p_cbf(i,j) = atheta'*x(:,j);
        p_mmd(i,j) = 1 / (atheta1'*Un1*Un1'*atheta1);
        p_smd(i,j) = 1 / (atheta2'*Un2*Un2'*atheta2);
        
    end
end

P_cbf = sum(p_cbf.' .* conj(p_cbf.')) / length(find(f>=fl & f<=fh)); 
P_mmd = sum(p_mmd.' .* conj(p_mmd.')) / length(find(f>=fl & f<=fh));
P_smd = sum(p_smd.' .* conj(p_smd.')) / length(find(f>=fl & f<=fh));

%% 画图
figure;
plot(theta,10*log10(P_cbf/max(P_cbf)));
hold on;
plot(theta,10*log10(P_mmd/max(P_mmd)));
plot(theta,10*log10(P_smd/max(P_smd)));
legend('CBF','MMD','SMD');
xlabel('角度/°');ylabel('幅度/dB');


