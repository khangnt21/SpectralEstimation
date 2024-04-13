%% SO SANH UOC LUONG PHO BARTLETT-WELCH VA BLACKMAN-TUKEY
%
% Chuong trinh co 2 goi y (prompt) yeu cau ta nhap loai 
% tin hieu quan sat x (thong so th) va cach tinh chieu 
% dai cua so dung de dieu chinh uoc luong tu tuong quan 
% Rx trong phuong phap Blackman-Tukey (thong so kBT). Ta 
% chon th = 1 de thi nghiem voi tin hieu dieu hoa, hoac 
% th = 2 voi tin hieu ARMA. Chon kBT = 1 neu muon cua so 
% dai bang 1/5 Rx (theo goi y cua Stephen M. Kay, Modern 
% spectral estimation: Theory and application, Prenctice-
% Hall, 1998, trang 81), hoac kBT = 2 neu muon cua so dai 
% bang Rx. 
% Cac thong so con lai khac da duoc chon san trong chuong 
% trinh. Tin hieu quan sat x co chieu dai Lx, co sai so 
% quan sat do nhieu cong voi ty le tin tren nhieu theo 
% don vi dB la SNRdB. Neu ta muon tin hieu sach (khong co 
% nhieu quan sat), chon gia tri du lon cua SNRdB (chang 
% han SNRdB = 1000), tuong duong voi SNR = vo han. Chia 
% nho x thanh cac tin hieu con co chieu dai Lb (phuong 
% phap Bartlett-Welch). Thong so phantram_Lb (co gia tri 
% tu 0 den 50) xac dinh phan du lieu trung lap giua hai 
% tin hieu con lien ke. Cua so dieu chinh tin hieu con la 
% loai cua so kw, co chieu dai Lb; kw = 1 (Chu nhat), 2 
% (Hann), 3 (Hamming), 4 (Blackmann) va 5 (Bartlett). 
% Ta cung dung thong so kw cho loai cua so dieu chinh uoc 
% luong tu tuong quan Rx trong phuong phap Blackman-
% Tukey; luc do chieu dai Rx la 2*Lb-1. 
%
% Chuong trinh lenh nay (sosanh_BW_BT) su dung them 5 
% chuong trinh ham sau:
%   [Wd,U,txtcs] = cuaso(Lw,kw)
%   [x,Sxx] = thdh(A,fnu,Lx,SNRdB)         
%   [x,Sxx] = thARMA(a,b,Lx,SNRdB)
%   [S_BWdB,SxxdB] = ppBW_Sxx(x,Lb,phantramLb,SNRdB,th,kw)
%   [S_BTdB,SxxdB] = ppBT_Sxx(x,Lb,phantramLb,SNRdB,th,kw,kBT)

% Viet cho giao trinh: 
% Xu ly tin hieu ngau nhien, Dai hoc Quoc gia Ha Noi, 2024
% Tac gia: Nguyen Linh Trung, Huynh Huu Tue
% ========================================================

%% Phan I: Nhap du lieu va tao tin hieu {x(n)}
% Tin hieu {x(n)} co chieu dai Lx voi SNRdB cho truoc. 
% Chu y su khac biet giua thang truc hoanh, cho tan so 
% so nu cua Sxx, cua tin hieu dieu hoa va tin hieu ARMA.

% Nhap thong tin dau vao cho th va kBT
th = input('Chon tin hieu ta muon phan tich (th = 1 la dieu hoa, th = 2 la ARMA): th = ')
kBT = input('Chon chieu dai cua so dieu chinh uoc luong tu tuong quan Rx (kBT = 1 la dai bang 1/5 Rx, = 2 la dai bang Rx): kBT = ')

% Chon chieu dai tin hieu
Lx = 200000; % 200, 500, 1000, 10000, 20000, ...

% Chon chieu dai tin hieu con
Lb = 200; % 200, 256, 500, 1000, 5000   
phantramLb = 50;    % 0-50%: trung lap giua 2 tin hieu con
% Nhieu quan sat
%SNRdB = 1000;       % khong co nhieu
%SNRdB = 10;        % co it nhieu
SNRdB = 0;         % co nhieu tuong duong tin hieu
%SNRdB = -10;       % co nhieu lon hon tin hieu

% Chon loai cua so
%kw = 1; cuaso = ['Chu nhat'];
kw = 2; cuaso = ['Hann'];
%kw = 3; cuaso = ['Hamming'];
%kw = 4; cuaso = ['Blackman'];
%kw = 5; cuaso = ['Bartlett'];

% Chon loai va tao tin hieu dung de thi nghiem. Chuong
% trinh co 2 loai tin hieu: dieu hoa va ARMA. Voi loai 
% dieu hoa, ta xet 2 truong hop: co 1 tan so (hinh sin) 
% va co 3 tan so. Voi loai ARMA, xet 2 truong hop: co pho 
% chua 3 dinh (de so sanh voi truong hop tin hieu dieu 
% hoa co 3 tan so), co pho thong thap. Tin hieu nao khong 
% dung thi khoa cac lenh lien quan.
if th == 1      % tin hieu dieu hoa
    % Tin hieu hinh sin
%    A = 1;
%    fnu = 0.123;
    % Tin hieu dieu hoa chua 3 tan so
    A = [1 1 1]';
    fnu = [0.1 0.12345 0.2]';
    [x,Sxx] = thdh(A,fnu,Lx,SNRdB); % tao tin hieu
    txt = ['dieu hoa'];
elseif th == 2  % tin hieu ARMA
    % Tin hieu ARMA co pho chua 3 dinh
%    b = 1;
%    a = [1.0000 -2.1248 2.2574 -1.7883 2.2125 -2.0411 0.9415];
    % Tin hieu ARMA co pho thong thap
    b = [0.0464 0.1829 0.2572 0.1549];
    a = [1 -0.8664 0.6630 -0.1514];
    txt = ['ARMA'];
    [x,Sxx] = thARMA(b,a,Lx,SNRdB); % tao tin hieu
end
SxxdB = 10*log10(Sxx);

%% Phan II: Uoc luong pho cong suat bang 2 phuong phap: 
%% Bartlett-Welch va Blackman-Tukey

%  Phuong phap Bartlett-Welch
S_BWdB = ppBW_Sxx(x,Lb,phantramLb,SNRdB,kw);
nBW = 1:length(S_BWdB)/2;
tBW = (nBW-1)/length(S_BWdB);

% Phuong phap Blackman-Tukey (co su dung them Welch)
S_BTdB = ppBT_Sxx(x,Lb,phantramLb,SNRdB,kw,kBT);
nBT = 1:length(S_BTdB)/2;
tBT = (nBT-1)/length(S_BTdB);
 
if th == 1      % tin hieu dieu hoa
    nSx = 1:length(SxxdB)/2;
    tSx = (nSx-1)/length(SxxdB);
elseif th==2    % tin hieu ARMA
    nSx = 1:length(SxxdB);
    tSx = .5*(nSx-1)/length(SxxdB);
end

%% Phan III: Hien thi ket qua, qua 2 do thi. Do thi 1 
%% bieu dien rieng le 3 pho mat do cong suat: chinh xac 
%% (khong co nhieu, lay truc tiep tu chuong trinh tao tin 
%% hieu), uoc luong theo Bartlett-Welch, uoc luong theo 
%% Blackman-Tukey. Do thi 2 gop ca ba pho trong do thi 1 
%% de de so sanh.

% Do thi 1
figure
subplot(311); 
plot(tSx,SxxdB(nSx),'k','linewidth',1.2);
legend('Chinh xac')
xlabel('Tan so chuan hoa \nu') 
ylabel('Pho cong suat (dB)')
if th == 1 
    title({['Uoc luong pho cong suat bang Bartlett-Welch va Blackman-Tukey'] 
        ['tin hieu = ', txt, '; Lx = ', num2str(Lx), '; Lb = ', num2str(Lb), '; phantramLb = ', num2str(phantramLb), '; SNR = ', num2str(SNRdB), ' dB']
        ['cua so ', cuaso, '; \nu = [',num2str(fnu'),']']})
else
    title({['Uoc luong pho cong suat bang 2 phuong phap Bartlett-Welch va Blackman-Tukey'] 
        ['tin hieu = ', txt, '; Lx = ', num2str(Lx), '; Lb = ', num2str(Lb), '; phantramLb = ', num2str(phantramLb), '; SNR = ', num2str(SNRdB), ' dB']
                ['cua so ', cuaso]})
end
subplot(312);
plot(tBW,S_BWdB(nBW),'b','linewidth',1.5);
%xlabel('Tan so chuan hoa \nu') 
%ylabel('Pho cong suat (dB)')
legend('Bartlett-Welch')
subplot(313); 
plot(tBT,S_BTdB(nBT),'r','linewidth',1.5);
if kBT == 1
    legend('Blackman-Tukey; Lw = Lr/5')
elseif kBT == 2
    legend('Blackman-Tukey; Lw = Lr')
end
xlabel('Tan so chuan hoa \nu')
ylabel('Pho cong suat (dB)')

% Do thi 2
figure
plot(tSx,SxxdB(nSx),'k:','linewidth',1.5); 
hold
plot(tBW,S_BWdB(nBW),'b-.',tBT,S_BTdB(nBT),'r-','linewidth',1.5); 
hold
if kBT == 1
    legend(['Chinh xac'], ['Bartlett-Welch'], ['Blackman-Tukey; Lw = Lr/5'])
else
    legend(['Chinh xac'], ['Bartlett-Welch'], ['Blackman-Tukey; Lw = Lr'])
end
if th == 1 
    title({['Uoc luong pho cong suat tin hieu ', txt, ' bang Bartlett-Welch va Blackman-Tukey']
        ['Lx = ', num2str(Lx), '; Lb = ', num2str(Lb), '; phantramLb = ', num2str(phantramLb), '; SNR = ', num2str(SNRdB), ' dB', '; cua so ', cuaso]
        ['tan so: \nu = [',num2str(fnu'),']']})
else
    title({['Uoc luong pho cong suat tin hieu ', txt, ' bang Bartlett-Welch va Blackman-Tukey']
        ['Lx = ', num2str(Lx), '; Lb = ', num2str(Lb), '; phantramLb = ', num2str(phantramLb), '; SNR = ', num2str(SNRdB), ' dB','; cua so ', cuaso]})
end
xlabel('Tan so chuan hoa \nu')
ylabel('Pho cong suat (dB)')