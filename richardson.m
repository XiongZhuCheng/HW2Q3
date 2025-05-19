function I_ext = richardson(f, a, b, m, n, p)
% 理查森外推法
% f: 被积函数句柄
% [a, b]: 积分区间
% m: 子区间数量
% n: 高斯点数
% p: 原方法的收敛阶

I_h = gauss_quad(f, a, b, m, n);
I_h2 = gauss_quad(f, a, b, 2*m, n);
I_ext = (2^p * I_h2 - I_h) / (2^p - 1);
end
