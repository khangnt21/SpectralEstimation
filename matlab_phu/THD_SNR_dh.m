%% ANH HUONG CUA NHIEU QUAN SAT LEN TUAN HOAN DO CUA 
% TIN HIEU DIEU HOA
%
% Tinh va bieu dien Tuan hoan do cua mot tin hieu dieu 
% hoa. Ta co dinh chieu dai Lb, tan so chuan hoa fnu 
% cua tin hieu, loai cua so, thay doi ty le tin tren nhieu 
% SNRdB (theo thang dB). Trong truong hop khong co nhieu 
% quan sat (SNRdB = vo cuc), chon gia tri du lon (nhu 
% SNRdB = 1000).

% Viet cho giao trinh: 
% Xu ly tin hieu ngau nhien, Dai hoc Quoc gia Ha Noi, 2024
% Tac gia: Nguyen Linh Trung, Huynh Huu Tue
% ========================================================
%%
% Chon mot trong hai tin hieu dieu hoa: co 1 tan so (hinh 
% sin) va co 3 tan so. Khoa cac dong lenh cua tin hieu 
% khong chon.
A = 1;          % Tin hieu dieu hoa 1 tan so (hinh sin)
fnu = .121;     
%A = [1 1 1]';  % Tin hieu dieu hoa 3 tan so
%fnu = [.1 .121 .2]';  
Lb = 500;       % chieu dai tin hieu (= chieu dai cua so)

% Chon cac muc SNRdB de danh gia anh huong cua nhieu quan 
% sat: nhieu lon, nhieu tuong duong tin hieu, nhieu it, 
% khong co nhieu.
SNRdBall = [-10 10 0 1000];

% Chon loai cua so
kw = 1; % cua so Chu nhat
%kw = 2; % cua so Hann
%kw = 3; % cua so Hamming
%kw = 4; % cua so Blackman
%kw = 5; % cua so Bartlett

%% Uoc luong pho bang Tuan hoan do va hien thi ket qua
figure
nn = 0:Lb/2;
nt = nn/Lb;
for k = 1:length(SNRdBall)
    % Tinh thong so nhieu
    SNRdB = SNRdBall(k);
    S = sum(A.^2)/2;
    N0 = S*10^(-SNRdB/10);
    sigma = sqrt(N0);

    % Tao tin hieu quan sat co nhieu cong
    xsach = sum([diag(A)*cos(2*pi*fnu*(0:Lb-1) + 2*pi*rand(length(A),1)); zeros(1,Lb)]);
    x = xsach + sigma*randn(1,Lb);

    % Tinh du lieu cua so
    [Wd,U,txtcs] = cuaso(Lb,kw);
    
    % Tinh tuan hoan do
    Pxx = (abs(fft(x.*Wd).^2))/(U*Lb);
    PxxdB = 10*log10(Pxx);
    
    % Hien thi ket qua
    subplot(2,2,k);
    plot(nt,PxxdB(1+nn));
    axis([0 .5 min(PxxdB)-10 max(PxxdB)+10])
    xlabel('Tan so chuan hoa, \nu')
    ylabel('Tuan hoan do, S_{xx} (dB)')
    SNRstr = [num2str(SNRdBall(k)) ' dB'];
    if SNRdB >= 1000
        SNRstr = '\infty';
    end
    legend(['SNR = ', SNRstr])
    if length(A) == 1
        title({['Tuan hoan do Pxx (dB) cua tin hieu hinh sin co tan so']
            ['\nu = ', num2str(fnu(1)), '; Lb = ', num2str(Lb), '; cua so ', txtcs]})
    else
        title({['Tuan hoan do Pxx (dB) cua tin hieu dieu hoa co 3 tan so']
            ['\nu = ', num2str(fnu(1)), ', ', num2str(fnu(2)), ', ', num2str(fnu(3)), '; Lb = ', num2str(Lb), '; cua so ', txtcs]})
    end
end