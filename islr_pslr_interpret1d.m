function output= islr_pslr_interpret1d(input,D1)
% ouput 输出数据
%input  输入数据
%D1插值倍数

B1=fftshift(fft(input));
B2=[zeros(1,floor(length(B1)*(D1-1)/2)) B1 zeros(1,floor(length(B1)*(D1-1)/2))];
output=abs(ifft(fftshift(B2)));%%方位向插值后到时域
end