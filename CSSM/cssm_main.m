% CSSM算法DOA估计示例
clear; clc; close all;
% 仿真参数设置
M = 8;              % 阵元数量
d = 0.75;            % 阵元间距(m)，设为半波长
fs = 8000;          % 采样频率
T = 1;              % 信号持续时间(s)
N = fs * T;         % 采样点数
c = 1500;            % 声速(m/s)
f_center = 1000;    % 信号中心频率
B = 500;            % 信号带宽
f = f_center - B/2 : 10 : f_center + B/2;  % 频率点

% 真实DOA
doa_true = [-30, 15];  % 两个相干信号的入射角度(度)
D = length(doa_true);  % 信号源数目

% 生成相干宽带信号
t = (0:N-1)/fs;
s = zeros(D, N);
for i = 1:D
    % 生成宽带信号( chirp信号 )
    s(i,:) = chirp(t, f_center - B/2, T, f_center + B/2);
end
% 使信号相干 (第二个信号是第一个信号的延迟版本)
s(2,:) = [zeros(1, 100), s(1, 1:end-100)];

% 构造阵列流形矩阵
A = zeros(M, D);
for i = 1:D
    A(:, i) = steering_vector(M, d, c, f_center, doa_true(i));
end

% 生成接收信号 (加入噪声)
x = A * s + 2 * randn(M, N);  % M×N的接收信号矩阵

% 角度搜索范围
theta_range = -90:0.1:90;

% 计算最佳聚焦频率
% f_candidates = f_center - B/2 : 10 : f_center + B/2;  % 频率点
% [f0_opt, ~] = optimal_focus_frequency(M, d, c, f_candidates, f, D, theta_range);


% 使用CSSM算法估计DOA
[doa_est, spectrum1, spectrum2] = cssm_doa(x, f, fs, M, d, theta_range, f_center);
% [doa_est, spectrum1] = tct_doa(x, f, fs, M, d, theta_range, f_center);

% 绘制空间谱
figure;
plot(theta_range, 10*log10(spectrum1/max(spectrum1)));
hold on;
plot(theta_range, 10*log10(spectrum2/max(spectrum2)));
grid on;
xlabel('角度(度)');
ylabel('归一化功率谱(dB)');
title('CSSM算法DOA估计');

fprintf('真实DOA: %.3f %.3f\n',doa_true);
fprintf('估计DOA: %.3f %.3f\n',doa_est);



