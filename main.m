% main.m
clear; clc;

% 被积函数
f = @(x) exp(x) .* sin(x);
a = -pi;
b = (7*pi)/6;

% 解析积分值
I_exact = analytical_integral(a, b);

% 高斯点数
n_values = 1:5;
% 子区间数量
m_values = 2.^(1:10);

errors = zeros(length(n_values), length(m_values));

for i = 1:length(n_values)
    n = n_values(i);
    for j = 1:length(m_values)
        m = m_values(j);
        I_num = gauss_quad(f, a, b, m, n);
        errors(i, j) = abs(I_num - I_exact);
    end
end

% 绘制误差图
plot_results(m_values, errors, n_values, '复合高斯求积法误差');

% 估计收敛阶
for i = 1:length(n_values)
    p = polyfit(log(1./m_values), log(errors(i, :)), 1);
    fprintf('n = %d 的收敛阶约为 %.2f\n', n_values(i), -p(1));
end

errors_ext = zeros(length(n_values), length(m_values));

for i = 1:length(n_values)
    n = n_values(i);
    p = 2 * n; % 高斯求积法的收敛阶为 2n
    for j = 1:length(m_values)
        m = m_values(j);
        I_ext = richardson(f, a, b, m, n, p);
        errors_ext(i, j) = abs(I_ext - I_exact);
    end
end

% 绘制误差图
plot_results(m_values, errors_ext, n_values, '理查森外推法误差');

% 估计收敛阶
for i = 1:length(n_values)
    p_fit = polyfit(log(1./m_values), log(errors_ext(i, :)), 1);
    fprintf('n = %d 的理查森外推法收敛阶约为 %.2f\n', n_values(i), -p_fit(1));
end

% 新的被积函数
f_new = @(x) exp(-x.^2) .* sin(x);
a_new = 0;
b_new = 2;

% 由于无法解析求积分，使用高精度数值积分作为参考值
I_exact_new = integral(f_new, a_new, b_new, 'AbsTol', 1e-14);

n = 5;
p = 2 * n;
errors_new = zeros(1, length(m_values));
errors_new_ext = zeros(1, length(m_values));

for j = 1:length(m_values)
    m = m_values(j);
    I_num = gauss_quad(f_new, a_new, b_new, m, n);
    I_ext = richardson(f_new, a_new, b_new, m, n, p);
    errors_new(j) = abs(I_num - I_exact_new);
    errors_new_ext(j) = abs(I_ext - I_exact_new);
end

% 绘制误差图
figure;
loglog(1./m_values, errors_new, '-o', 'DisplayName', '原始方法');
hold on;
loglog(1./m_values, errors_new_ext, '-s', 'DisplayName', '理查森外推');
xlabel('子区间长度 h');
ylabel('误差');
title('非多项式函数积分误差');
legend('Location', 'best');
grid on;

% 估计收敛阶
p_fit = polyfit(log(1./m_values), log(errors_new), 1);
p_fit_ext = polyfit(log(1./m_values), log(errors_new_ext), 1);
fprintf('原始方法收敛阶约为 %.2f\n', -p_fit(1));
fprintf('理查森外推法收敛阶约为 %.2f\n', -p_fit_ext(1));
