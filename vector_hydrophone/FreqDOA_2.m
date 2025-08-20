%%% 长时间信号分帧DOA估计与方位历程图绘制
clear; clc; close all;

%% 参数设置
fs = 2048;              % 采样率
T_total = 20;           % 总时长(秒) - 长时间信号
N_total = fs * T_total; % 总采样点数
f0 = 200;               % 信号频率(Hz)
SNR = 15;               % 信噪比(dB)

% 目标运动轨迹（模拟两个目标的方位角随时间变化）
t_total = 0:1/fs:(N_total-1)/fs;
theta1 = 30 + 20*sin(2*pi*0.05*t_total);  % 目标1: 30°±20° 缓慢摆动
theta2 = 250 - 15*sin(2*pi*0.03*t_total); % 目标2: 150°±15° 缓慢摆动

% 分帧参数
frame_len = 2*fs;       % 帧长
overlap = 0.75;         % 重叠率 (75%重叠)
step = round(frame_len*(1-overlap)); % 帧移
nfft = 2048;            % FFT点数
theta_search = 0:0.5:359.5; % 搜索角度范围

%% 生成长时间矢量水听器信号
% 生成两个运动目标的信号
s1 = cos(2*pi*f0*t_total + 0.1*t_total);  % 目标1信号
s2 = cos(2*pi*f0*t_total + 0.2*t_total + pi/3);  % 目标2信号

% 声压通道信号
p = s1 + s2;

% 振速通道信号（x和y方向）
vx = zeros(size(t_total));
vy = zeros(size(t_total));
for i = 1:length(t_total)
    vx(i) = s1(i)*cosd(theta1(i)) + s2(i)*cosd(theta2(i));
    vy(i) = s1(i)*sind(theta1(i)) + s2(i)*sind(theta2(i));
end

% 添加高斯噪声
p = awgn(p, SNR, 'measured');
vx = awgn(vx, SNR, 'measured');
vy = awgn(vy, SNR, 'measured');

%% 分帧处理并进行DOA估计
num_frames = floor((N_total - frame_len)/step) + 1;  % 总帧数
t_frames = (0:num_frames-1)*step/fs;  % 每帧对应的时间

% 存储三种算法的DOA估计结果（每帧可能有两个目标）
doa_cbf = zeros(num_frames, length(theta_search));
doa_mvdr = zeros(num_frames, length(theta_search));
doa_music = zeros(num_frames, length(theta_search));

% 加窗函数
win = hanning(frame_len);

% 进度显示
progress = waitbar(0, '正在处理...');

for frame = 1:num_frames
    % 更新进度条
    waitbar(frame/num_frames, progress, sprintf('正在处理第%d/%d帧...', frame, num_frames));
    
    % 提取当前帧
    idx = (frame-1)*step + 1 : (frame-1)*step + frame_len;
    p_frame = p(idx) .* win;
    vx_frame = vx(idx) .* win;
    vy_frame = vy(idx) .* win;
    
    % FFT变换
    P = fft(p_frame, nfft);
    Vx = fft(vx_frame, nfft);
    Vy = fft(vy_frame, nfft);
    
    % 选择信号频率附近的频点
    f = (0:nfft-1)*fs/nfft;
    f_idx = find(f >= f0-10 & f <= f0+10);  % 围绕目标频率的频带
    num_freq = length(f_idx);

    % 计算三种算法的空间谱
    P_cbf = zeros(length(theta_search), num_freq);
    P_mvdr = zeros(length(theta_search), num_freq);
    P_music = zeros(length(theta_search), num_freq);

    for fk = 1:num_freq
        f_idx_current = f_idx(fk);
        % 计算当前帧的协方差矩阵
        X = [P(f_idx_current,:).', Vx(f_idx_current,:).', Vy(f_idx_current,:).'];  % 频点×通道
        R = X' * X / num_frames;  % 协方差矩阵
        R = (R + R')/2;  % 确保Hermitian对称
    
        % 特征值分解（用于MUSIC）
        [U, D] = eig(R);
        U_n = U(:, 1);  % 噪声子空间（假设2个目标）
            
        for theta_idx = 1:length(theta_search)
            theta = theta_search(theta_idx);
            a = [1; cosd(theta); sind(theta)];  % 导向矢量
            % CBF
            P_cbf(theta_idx) = real(a' * R * a);
            
            % MVDR（添加正则化）
            R_inv = inv(R + 1e-6*eye(size(R)));
            P_mvdr(theta_idx) = 1 / real(a' * R_inv * a);
            
            % MUSIC
            P_music(theta_idx) = 1 / real(a' * U_n * U_n' * a);
        end
    end
    power_weights = abs(P(f_idx, :)) * ones(size(P,2),1) / num_frames;  % 每个频率点的功率
    power_weights = power_weights / sum(power_weights);  % 归一化
    % 加权融合空间谱
    P_cbf_avg = P_cbf * power_weights;
    P_mvdr_avg = P_mvdr * power_weights;
    P_music_avg = P_music * power_weights;

    doa_cbf(frame, :) = P_cbf_avg';
    doa_mvdr(frame, :) = P_mvdr_avg';
    doa_music(frame, :) = P_music_avg';
end
close(progress);

%% 绘制方位历程图
figure;
pcolor(theta_search, ( 1:num_frames )*fs/frame_len, 10*log10(doa_cbf) );
shading interp;title('CBF');
xlabel('theta');ylabel('时间');

figure;
pcolor(theta_search, ( 1:num_frames )*fs/frame_len, 10*log10(doa_mvdr) );
shading interp;title('MVDR');
xlabel('theta');ylabel('时间');

figure;
pcolor(theta_search, ( 1:num_frames )*fs/frame_len, 10*log10(doa_music) );
shading interp;title('MUSIC');
xlabel('theta');ylabel('时间');

