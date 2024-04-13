function [S_BT,S_BTdB,Rx] = ppBT_Sxx(x,Lb,phantramLb,SNRdB,kw,kBT)

%% UOC LUONG PHO VOI PHUONG PHAP BLACKMAN-TUKEY
%
% [S_BT,S_BTdB,Rx] = ppBT_Sxx(x,Lb,phantramLb,SNRdB,kw,kBT)
%
% Voi tap du lieu x cua tin hieu quan sat can uoc luong
% pho, ta da co chieu dai cua tin hieu la Lx. Neu ket hop 
% voi phuong phap Bartlett, ta chia du lieu thanh cac tap 
% du lieu con, co chieu dai la Lb. Neu chon Lb = Lx, tuc 
% khong chia nho, ta co phuong phap uoc luong pho Tuan 
% hoan do. Neu ket hop them Welch, ta su dung thong so 
% phantram_Lb (tu 0 den 50) cho chieu dai cua phan du 
% lieu trung lap. Neu da chon Lb = Lx, thi phantramLb 
% duoc gan bang 0. SNRdB la ty le tin tren nhieu (dB). 
% Chon cua so thong dung ap vao ham tu tuong quan Rx bang 
% thong so kw cho loai cua so va thong so kBT cho chieu 
% dai cua so. Cac cua so thong dung la Chu nhat, Hann, 
% Hamming, Blackmann va Bartlett, co thong so kw = 1, 2, 
% 3, 4, 5 tuong ung. Co 2 cach chon chieu dai kBT cua so 
% kw ap vao Rx trong phuong phap Blackman-Tukey; kBT = 1 
% cho truong hop chieu dai cua so bang 1/5 lan chieu dai 
% cua Rx (theo goi y cua: SM Kay, Modern spectrum 
% estimation: Theory and application, Prentice-Hall, 
% 1998; trang 81), kBT = 2 cho chieu dai bang chieu dai 
% cua Rx. Chuong trinh cho dau ra la uoc luong pho S_BTdB 
% theo phuong phap Blackman-Tukey co don vi la dB va uoc 
% luong ham tu tuong quan Rx.
% 
% Chuong trinh nay (ppBT_Sxx) can su dung 3 chuong trinh 
% ham sau:
%  [Wd,U,txtcs] = cuaso(Lw,kw): tinh cua so de ap vao Rx
%  [x,Sxx] = thdh(Lx,SNRdB): tao tin hieu dieu hoa x
%  [x,Sxx] = thARMA(b,a,Lx,SNRdB): tao tin hieu ARMA x

% Viet cho giao trinh: 
% Xu ly tin hieu ngau nhien, Dai hoc Quoc gia Ha Noi, 2024
% Tac gia: Nguyen Linh Trung, Huynh Huu Tue
% ========================================================
%%
% Tinh cua so Wd de ap vao ham tu tuong quan Rx
Lr = 2*Lb-1; % chieu dai (so le) cua ham tu tuong quan Rx
if kBT == 1  % tao cua so co chieu dai bang 1/5 cua Rx
    Lw = floor(Lr/5);
    Lw = Lw + (floor(Lw/2) == Lw/2); % Lw la so le
    [Wd,U] = cuaso(Lw,kw); % tao cua so loai kw
    D0 = .5*(Lr-1) + 1;     % vi tri dinh cua Rx
    d0 = .5*(Lw-1) + 1;    % vi tri dinh cua cua so Wd
    D0d0 = D0 - d0;         % khoang cach giua 2 dinh 
    % Them chuoi 0 vao Wd sao cho Wd co cung chieu dai 
    % voi Rx va co dinh trung voi dinh cua Rx
    Wd = [zeros(1,D0d0),Wd,zeros(1,D0d0)];   
elseif kBT == 2  % tao cua so co cung chieu dai voi Rx
    [Wd,U,txtcs] = cuaso(Lr,kw);   
end

% Ap dung them Welch: Cho phep cac tap du lieu con trung 
% lap. Xac dinh so mau du lieu trung lap dd giua 2 tap 
% con lien ke. 
if Lb == size(x)        % Lb = Lx, khong ap dung Bartlett
   phantramLb = 0;
end
dd = phantramLb*Lb/100; % dau vao phantramLb = 0-50
D = floor(Lb-dd);       % du lieu khong trung lap

% Ap dung Blackman-Tukey: Tinh uoc luong pho cua moi tap 
% du lieu con bang cach lay bien doi Fourier (FFT) cua 
% uoc luong cua ham tu tuong quan Rx cua tap du lieu con.
% Trong do, ta co ap dung them Bartlett, chia tap du lieu 
% x thanh k tap con.
S_BT = zeros(1,Lr);
k = 0;
while  k*D <= length(x)-Lb
    % Tao tap du lieu con thu k, k*D+1 la diem bat dau
    X = x(k*D+1:k*D+Lb);
    % Uoc luong ham tu tuong quan cua X
    Rx = conv(X,X(end:-1:1))/length(X); 
    % Uoc luong pho bang FFT cua Rx, co ket hop cua so
    R_hat = Rx.*Wd;
    S_BT = S_BT+abs(fft(R_hat));
    k = k+1;
end

% Ap dung them Bartlett: Uoc luong pho cua tap du lieu x
S_BT = S_BT/k; 
S_BTdB = 10*log10(S_BT);  % theo thang dB