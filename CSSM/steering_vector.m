function a = steering_vector(M, d, c, f, theta)
    % 构造均匀线阵的导向矢量
    % M: 阵元数, d: 阵元间距, c: 波速, f: 频率, theta: 入射角度(度)
    theta_rad = theta * pi / 180;  % 转换为弧度
    a = exp(-1i * 2 * pi * f * d * (0:M-1)' * sin(theta_rad) / c);
end