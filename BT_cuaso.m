%% ANH HUONG CUA LOAI CUA SO LEN PHUONG PHAP BLACKMAN-TUCKEY
%
% Tinh va ve tuan hoan do cua mot tin hieu dieu hoa. Ta 
% co dinh tan so chuan hoa fnu va ty le tin tren nhieu 
% SNRdB (theo thang dB), thay doi chieu dai cua so Lb,
% thay doi cac loai cua so, cach chon do dai cua so Rx.

% Su dung trong bao cao cua Nguyen The Khang
% ========================================================
%%
%% Phan I: Nhap du lieu va tao tin hieu {x(n)}
% Tin hieu {x(n)} co chieu dai Lx voi SNRdB cho truoc. 
% Chu y su khac biet giua thang truc hoanh, cho tan so 
% so nu cua Sxx, cua tin hieu dieu hoa va tin hieu ARMA.

% Nhap thong tin dau vao cho th va kBT
th = input('Chon tin hieu ta muon phan tich (th = 1 la dieu hoa loai 1, th = 2 la dieu hoa loai 2, th = 3 la ARMA): th = ');
kBT = input('Chon chieu dai cua so dieu chinh uoc luong tu tuong quan Rx (kBT = 1 la dai bang 1/5 Rx, = 2 la dai bang Rx): kBT = ');
snr = input('Chon muc do nhieu quan sat (snr = 1 la SNRdB = 10dB, snr = 2 la SNRdB = -10dB): snr = ');
% Chon chieu dai tin hieu
Lx = 20000; % 200, 500, 1000, 10000, 20000, ...

% Chon chieu dai tin hieu con
Lb = 200; % 200, 256, 500, 1000, 5000   
phantramLb = 50;    % 0-50%: trung lap giua 2 tin hieu con

% Cai dat tham so nhieu quan sat
if snr == 1
    SNRdB = 10;
elseif snr == 2 
    SNRdB = -10;
end
% SNRdB = 0;
% Chon loai va tao tin hieu dung de thi nghiem. 
% Chuong trinh co 2 loai tin hieu: dieu hoa va ARMA. 
% Voi loai dieu hoa, ta xet 2 truong hop: 
%   co 1 tan so (hinh sin) 
%   co 3 tan so. 
% Voi loai ARMA, xet 2 truong hop: 
%   co pho chua 3 dinh (de so sanh voi truong hop tin hieu dieu 
% hoa co 3 tan so), 
%   co pho thong thap. 
% Tin hieu nao khong dung thi khoa cac lenh lien quan.

if th == 1      % tin hieu dieu hoa loai 1
    % Tin hieu dieu hoa chua 3 tan so
    A = [1 1 1]';
    fnu = [0.1 0.22345 0.3]';
    [x,Sxx] = thdh(A,fnu,Lx,SNRdB); % tao tin hieu
    txt = ['dieu hoa 1'];
elseif th == 2      % tin hieu dieu hoa loai 2
    % Tin hieu dieu hoa chua 3 tan so
    A = [1 1 1]';
    fnu = [0.1 0.12345 0.2]';
    [x,Sxx] = thdh(A,fnu,Lx,SNRdB); % tao tin hieu
    txt = ['dieu hoa 2'];
elseif th == 3  % tin hieu ARMA
    % Tin hieu ARMA co pho thong thap
    b = [0.0464 0.1829 0.2572 0.1549];
    a = [1 -0.8664 0.6630 -0.1514];
    txt = ['ARMA'];
    [x,Sxx] = thARMA(b,a,Lx,SNRdB); % tao tin hieu
end
SxxdB = 10*log10(Sxx);

if th == 1 || th == 2  % tin hieu dieu hoa
    nSx = 1:length(SxxdB)/2;
    tSx = (nSx-1)/length(SxxdB);
elseif th==3    % tin hieu ARMA
    nSx = 1:length(SxxdB);
    tSx = .5*(nSx-1)/length(SxxdB);
end
if kBT == 1
    txtlen = ['Rx/5'];
elseif kBT == 2 
    txtlen = ['Rx'];
end

%% Phan II: Quan sat pho cong suat uoc luong va danh gia
figure;
subplot(2,3,1)
plot(tSx,SxxdB(nSx),'k','linewidth',1.5); grid;
legend('Chinh xac')
xlabel('Tan so chuan hoa \nu') 
ylabel('Pho cong suat (dB)')
title({['Uoc luong pho cong suat tin hieu ', txt, ' bang Blackman-Tukey;']
        ['Lx = ', num2str(Lx), '; Lb = ', num2str(Lb), '; phantramLb = ', num2str(phantramLb), '; SNR = ', num2str(SNRdB), ' dB;']
        ['tan so: \nu = [',num2str(fnu'),']', ' ;do dai cua so = ', txtlen]})
% Chon cac loai cua so
% kw = 1; % cua so Chu nhat
% kw = 2; % cua so Hann
% kw = 3; % cua so Hamming
% kw = 4; % cua so Blackman
% kw = 5; % cua so Bartlett
kww = [1 2 3 4 5];
mat = zeros(1,399);
for kw = kww
% Phuong phap Blackman-Tukey (co su dung them Welch)
    if kw == 1 
            txtcs = ['Chu nhat'];
        elseif kw==2 
            txtcs = ['Hann'];
        elseif kw==3 
            txtcs = ['Hamming']; 
        elseif kw==4
            txtcs = ['Blackman'];
        elseif kw==5
            txtcs = ['Bartlett']; 
    end
    subplot(2,3,kw+1)
    [S_BT,S_BTdB,Rx] = ppBT_Sxx(x,Lb,phantramLb,SNRdB,kw,kBT);
    nBT = 1:length(S_BT)/2;
    tBT = (nBT-1)/length(S_BT);
%     plot(tSx,SxxdB(nSx),'k:','linewidth',1.5); 
    plot(tBT,S_BTdB(nBT),'r-','linewidth',1.5); grid;
    if kBT == 1
        legend(['', txtcs])
    else
        legend(['', txtcs])
    end
    xlabel('Tan so chuan hoa \nu')
    ylabel('Pho cong suat (dB)')
    mat = [mat; S_BTdB];
end
figure;
% plot(tSx,Sxx(nSx),'k','linewidth',1.5); hold on;
mat(1,:) = [];

plot(tBT,mat(1,nBT),'Color', [0 0.4470 0.7410],'linewidth',1.5); hold on;
plot(tBT,mat(2,nBT),'Color', [0.8500 0.3250 0.0980],'linewidth',1.5); hold on;
plot(tBT,mat(3,nBT),'Color', [0.9290 0.6940 0.1250],'linewidth',1.5); hold on;
plot(tBT,mat(4,nBT),'Color', [0.4940 0.1840 0.5560],'linewidth',1.5); hold on;
plot(tBT,mat(5,nBT),'Color', [0.4660 0.6740 0.1880],'linewidth',1.5); hold off;
legend(['Chu nhat'],['Han'],['Hamming'],['Blackman'],['Barlett'])
grid;
xlabel('Tan so chuan hoa \nu')
ylabel('Pho cong suat (dB)')
title({['Uoc luong pho cong suat tin hieu ', txt, ' bang Blackman-Tukey;']
        ['Lx = ', num2str(Lx), '; Lb = ', num2str(Lb), '; phantramLb = ', num2str(phantramLb), '; SNR = ', num2str(SNRdB), ' dB;']
        ['tan so: \nu = [',num2str(fnu'),']', '; do dai cua so = ', txtlen]})
