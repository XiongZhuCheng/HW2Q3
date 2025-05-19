function plot_results(m_values, errors, n_values, title_str)
% 绘制误差图
figure;
loglog(1./m_values, errors', '-o');
xlabel('子区间长度 h');
ylabel('误差');
legend(arrayfun(@(n) sprintf('n = %d', n), n_values, 'UniformOutput', false));
title(title_str);
grid on;
end
