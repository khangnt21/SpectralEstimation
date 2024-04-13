%% ANH HUONG CUA LOAI CUA SO LEN TUAN HOAN DO
%
% Tinh va ve tuan hoan do cua mot tin hieu dieu hoa. Ta 
% co dinh tan so chuan hoa fnu va ty le tin tren nhieu 
% SNRdB (theo thang dB), thay doi chieu dai cua so Lb.

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

% Tinh thong so nhieu
SNRdB = 10;
S = sum(A.^2)/2;
N0 = S*10^(-SNRdB/10);
sigma = sqrt(N0);
SNRstr = [num2str(SNRdB) ' dB'];
if SNRdB >= 1000
    SNRstr = '\infty';
end

% Chon cac loai cua so
% kw = 1; % cua so Chu nhat
% kw = 2; % cua so Hann
% kw = 3; % cua so Hamming
% kw = 4; % cua so Blackman
% kw = 5; % cua so Bartlett
% kww = [1 2]; 
kww = [1 5]; 
%% Tinh va hien thi tuan hoan do theo cac loai cua so
figure
for k = 1:length(kww) % cho cac loai cua so quan tam
    % Tao tin hieu quan sat co nhieu cong
    xsach = sum([diag(A)*cos(2*pi*fnu*(0:Lb-1) + 2*pi*rand(length(A),1)); zeros(1,Lb)]);
    x = xsach + sigma*randn(1,Lb);

    % Tinh du lieu cua so
    [Wd,U,txtcs] = cuaso(Lb,kww(k)); 

    % Tinh tuan hoan do
    Pxx = (abs(fft(x.*Wd).^2))/(U*Lb);
    PxxdB = 10*log10(Pxx);

    % Hien thi tuan hoan do theo thang tuyen tinh
    nn = 0:Lb/2;
    nt = nn/Lb;
    subplot(length(kww),2,2*k-1);
    plot(nt,Pxx(1+nn));
    axis([0 .5 min(Pxx)-5 max(Pxx)+5])
    xlabel('Tan so chuan hoa, \nu')
    ylabel('Tuan hoan do, Pxx')
    legend(['Cua so = ', txtcs])
    if length(A) == 1
        title({['Tuan hoan do Pxx cua tin hieu hinh sin co tan so']
            ['\nu = ', num2str(fnu(1)), '; Lb = ', num2str(Lb), '; SNR = ', SNRstr]})
    else
        title({['Tuan hoan do Pxx cua tin hieu dieu hoa co 3 tan so']
            ['\nu = ', num2str(fnu(1)), ', ', num2str(fnu(2)), ', ', num2str(fnu(3)), '; Lb = ', num2str(Lb), '; SNR = ', SNRstr]})
    end
    
    % Hien thi tuan hoan do theo thang dB
    subplot(length(kww),2,2*k);
    plot(nt,PxxdB(1+nn));
    axis([0 .5 min(PxxdB)-10 max(PxxdB)+10])
    xlabel('Tan so chuan hoa, \nu')
    ylabel('Tuan hoan do, P_{xx} (dB)')
    legend(['Cua so = ', txtcs])
    if length(A) == 1
        title({['Tuan hoan do Pxx (dB) cua tin hieu hinh sin co tan so']
            ['\nu = ', num2str(fnu(1)), '; Lb = ', num2str(Lb), '; SNR = ', SNRstr]})
    else
        title({['Tuan hoan do Pxx (dB) cua tin hieu dieu hoa co 3 tan so']
            ['\nu = ', num2str(fnu(1)), ', ', num2str(fnu(2)), ', ', num2str(fnu(3)), '; Lb = ', num2str(Lb), '; SNR = ', SNRstr]})
    end
end