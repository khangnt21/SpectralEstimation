function [S_BWdB,txtcs] = ppBW_Sxx(x,Lb,phantramLb,SNRdB,kw)

%% UOC LUONG PHO VOI PHUONG PHAP BARTLETT-WELCH
%
% [S_BWdB,txt] = ppBW_Sxx(x,Lb,phantramLb,SNRdB,kw)
%
% Voi tap du lieu x cua tin hieu quan sat can uoc luong 
% pho, ta da co chieu dai cua tin hieu la Lx. Ta chia du 
% lieu thanh cac tap du lieu con, co chieu dai la Lb, 
% theo Bartlett. Neu chon Lb = Lx, tuc khong chia nho, ta 
% co phuong phap uoc luong pho Tuan hoan do. Su dung 
% thong so phantram_Lb (tu 0 den 50) cho chieu dai cua 
% phan du lieu trung lap, theo Welch. SNRdB la ty le tin 
% tren nhieu (dB). Chon cua so thong dung ap vao tin hieu 
% quan sat x bang thong so kw cho loai cua so va thong so 
% kBT cho chieu dai cua so. Cac cua so thong dung la Chu 
% nhat, Hann, Hamming, Blackmann va Bartlett, co thong so 
% kw = 1, 2, 3, 4, 5 tuong ung. Chuong trinh cho dau ra 
% la uoc luong pho S_BWdB theo phuong phap Bartlett-Welch 
% co don vi la dB va ten cua cua so ap dung.
% 
% Chuong trinh nay (ppBW_Sxx) can su dung 3 chuong trinh 
% ham sau:
%  [Wd,U,txtcs] = cuaso(Lw,kw): tinh cua so de ap vao Rx
%  [x,Sxx] = thdh(Lx,SNRdB): tao tin hieu dieu hoa x 
%  [x,Sxx] = thARMA(b,a,Lx,SNRdB): tao tin hieu ARMA x 

% Viet cho giao trinh: 
% Xu ly tin hieu ngau nhien, Dai hoc Quoc gia Ha Noi, 2024
% Tac gia: Nguyen Linh Trung, Huynh Huu Tue
% ========================================================

% Tinh du lieu cua so Wd va thong so dieu chinh U.
% Wd co chieu dai Lb bang chieu dai tin hieu con. 
[Wd,U,txtcs] = cuaso(Lb,kw);   

% Ap dung Welch: Cho phep cac tap du lieu con trung lap. 
% Tinh so mau du lieu trung lap dd. 
dd = phantramLb*Lb/100; 
D = floor(Lb-dd);

% Ap dung Bartlett: Chia tap du lieu thanh k tap cong va
% tinh tuan hoan do cua chung
k = 0;
S_hat = zeros(1,Lb); 
while  k*D <= length(x)-Lb
    xcon = x(k*D+1:k*D+Lb).*Wd; % tin hieu con thu k+1
    S_hat = S_hat+(abs(fft(xcon)).^2)/Lb;
    k = k+1;
end

% Ap dung Bartlett: Uoc luong pho cua tap du lieu x
S_hat = S_hat/(U*k); 
S_BWdB = 10*log10(S_hat);  % theo thang dB