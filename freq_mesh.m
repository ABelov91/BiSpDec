function [h,f,arg,arg_half] = freq_mesh( N, left, right, center, width )
A = 1 - width/(right-left);

step = (right - left)/N;
arg = left:step:right;
f = zeros(1,N+1);
f_half = zeros(1,N);
arg_half = zeros(1,N);
h = zeros(1,N);
for n = 1:N+1
    f(n) = arg(n) - (A*width)*atan( (arg(n)-center)/width );
end
for n = 1:N
    temp = 0.5*(arg(n)+arg(n+1));
    arg_half(n) = temp;
    h(n) = ( 1 - A*width^2/( (temp-center)^2 + width^2 ) )*step;
    f_half(n) = ( arg_half(n) - (A*width)*atan( (arg_half(n)-center)/width ) );
end

% figure; hold on;
% plot(arg,f,'-ok')
% plot(arg_half,h,'-ok')
end