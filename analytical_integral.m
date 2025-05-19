function I = analytical_integral(a, b)
% 计算 f(x) = e^x * sin(x) 的解析积分值
F = @(x) (exp(x) .* (sin(x) - cos(x))) / 2;
I = F(b) - F(a);
end
