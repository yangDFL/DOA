%%% 矢量水听器频域波束形成DOA估计算法
clear; clc; close all;

%% 参数设置
fs = 2048;          % 采样率
T = 20;              % 信号时长(s)
N = fs*T;           % 总采样点数
f0 = 200;           % 信号频率(Hz)
theta1 = -30;        % 目标1方位角(°)
theta2 = 120;       % 目标2方位角(°)
SNR = 10;           % 信噪比(dB)
win_len = 1024;     % 帧长
overlap = 0.5;      % 重叠率
nfft = 2048;        % FFT点数
theta_search = -180:0.5:180;  % 搜索角度范围

%% 生成矢量水听器信号
t = 0:1/fs:(N-1)/fs;
% 目标信号（声压通道）
s1 = cos(2*pi*200*t);
s2 = cos(2*pi*150*t + pi/4);  % 带初始相位差的第二目标
p = s1 + s2;  % 声压通道信号

% 振速通道信号（x和y方向）
vx = s1*cosd(theta1) + s2*cosd(theta2);
vy = s1*sind(theta1) + s2*sind(theta2);

% 加噪声
p = awgn(p, SNR, 'measured');
vx = awgn(vx, SNR, 'measured');
vy = awgn(vy, SNR, 'measured');

%% 分帧加窗并转换到频域
win = hanning(win_len);  % 汉宁窗
step = round(win_len*(1-overlap));
num_frames = floor((N - win_len)/step) + 1;

% 存储各通道频域信号（频率点×帧数）
P = zeros(nfft, num_frames);
Vx = zeros(nfft, num_frames);
Vy = zeros(nfft, num_frames);

for k = 1:num_frames
    idx = (k-1)*step + 1 : (k-1)*step + win_len;
    p_frame = p(idx) .* win';
    vx_frame = vx(idx) .* win';
    vy_frame = vy(idx) .* win';
    
    % FFT变换（取正频率）
    P(:,k) = fft(p_frame, nfft);
    Vx(:,k) = fft(vx_frame, nfft);
    Vy(:,k) = fft(vy_frame, nfft);
end

%% 频域波束形成
% 选择信号频率附近的频点（减少计算量）
f = (0:nfft-1)*fs/nfft;
f_idx = find(f >= 100 & f <= 250);  % 围绕目标频率的频带
num_freq = length(f_idx);

% 初始化三种算法的空间谱（角度×频率）
P_cbf = zeros(length(theta_search), num_freq);
P_mvdr = zeros(length(theta_search), num_freq);
P_music = zeros(length(theta_search), num_freq);

for fk = 1:num_freq
    f_idx_current = f_idx(fk);
    % 提取当前频率点的三通道信号（帧数×通道数）
    X = [P(f_idx_current,:).', Vx(f_idx_current,:).', Vy(f_idx_current,:).'];
    
    % 计算协方差矩阵
    R = X' * X / num_frames;  % 样本协方差矩阵
    R = (R + R')/2;  % 确保Hermitian对称
    
    % 特征值分解（用于MUSIC）
    [U, D] = eig(R);
    U_n = U(:, 1);  % 3通道，取后1个特征向量为噪声子空间
    
    % 遍历所有搜索角度
    for theta_idx = 1:length(theta_search)
        theta = theta_search(theta_idx);
        a = [1; cosd(theta); sind(theta)];  % 导向矢量
        
        % CBF
        P_cbf(theta_idx, fk) = real(a' * R * a);
        
        % MVDR（添加正则化避免矩阵奇异）
        R_inv = inv(R + 1e-6*eye(size(R)));
        P_mvdr(theta_idx, fk) = 1 / real(a' * R_inv * a);
        
        % MUSIC
        P_music(theta_idx, fk) = 1 / real(a' * U_n * U_n' * a);
    end
end

%% 结果融合（对频率加权平均）
% 以信号功率为权重（声压频谱幅度）
power_weights = abs(P(f_idx, :)) * ones(size(P,2),1) / num_frames;  % 每个频率点的功率
power_weights = power_weights / sum(power_weights);  % 归一化

% 加权融合空间谱
P_cbf_avg = P_cbf * power_weights;
P_mvdr_avg = P_mvdr * power_weights;
P_music_avg = P_music * power_weights;

%% 可视化结果
figure('Position', [100 100 1000 800]);

% 1. CBF结果
subplot(3,1,1);
plot(theta_search, P_cbf_avg, 'LineWidth', 1.5);
title('常规波束形成（CBF）空间谱', 'FontSize', 12);
xlabel('方位角(°)'); ylabel('谱值');
xlim([-180 180]); grid on;
hold on; plot([theta1 theta1], ylim, 'r--', [theta2 theta2], ylim, 'r--');

% 2. MVDR结果
subplot(3,1,2);
plot(theta_search, P_mvdr_avg, 'LineWidth', 1.5);
title('MVDR波束形成空间谱', 'FontSize', 12);
xlabel('方位角(°)'); ylabel('谱值');
xlim([-180 180]); grid on;
hold on; plot([theta1 theta1], ylim, 'r--', [theta2 theta2], ylim, 'r--');

% 3. MUSIC结果
subplot(3,1,3);
plot(theta_search, P_music_avg, 'LineWidth', 1.5);
title('MUSIC算法空间谱', 'FontSize', 12);
xlabel('方位角(°)'); ylabel('谱值');
xlim([-180 180]); grid on;
hold on; plot([theta1 theta1], ylim, 'r--', [theta2 theta2], ylim, 'r--');

sgtitle('矢量水听器频域DOA估计结果（红色虚线为真实方位）', 'FontSize', 14);
