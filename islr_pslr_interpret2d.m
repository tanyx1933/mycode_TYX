function output= islr_pslr_interpret2d(input,D2)
% ouput 输出数据
%input  输入数据
%D2插值倍数

 E1=fftshift(fft(input),1);%先距离向
 E2=fftshift(fft(E1.'),1);
 [row,col]=size(E2);
 A3=[zeros(row,floor(col*(D2-1)/2)) E2 zeros(row,(D2-1)*col-floor(col*(D2-1)/2))];

 A4=[zeros(D2*col,floor(row*(D2-1)/2)) A3.' zeros(D2*col,(D2-1)*row-floor(row*(D2-1)/2))];
 A4=A4.';
 A5=ifft(fftshift(A4,1));
 output=ifft(fftshift(A5.',1));%%%转回时域；
end