function I = gauss_quad(f, a, b, m, n)
% 复合高斯求积法
% f: 被积函数句柄
% [a, b]: 积分区间
% m: 子区间数量
% n: 高斯点数（1 ≤ n ≤ 5）

% 获取标准区间 [-1,1] 上的高斯节点和权重
[xi, wi] = gauss_legendre(n);

% 子区间长度
h = (b - a) / m;
I = 0;

for k = 0:(m-1)
    xk = a + k*h;
    xk1 = xk + h;
    % 映射到当前子区间
    x = ((xk1 - xk)/2) * xi + (xk1 + xk)/2;
    I = I + ((xk1 - xk)/2) * sum(wi .* f(x));
end
end
