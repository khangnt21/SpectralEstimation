%% ANH HUONG CUA CHIEU DAI TIN HIEU LEN TUAN HOAN DO
%
% Tinh va bieu dien Tuan hoan do cua mot tin hieu dieu 
% hoa. Ta co dinh tan so fnu cua tin hieu va ty le tin
% tren nhieu SNRdB, thay doi chieu dai Lb cua tin hieu. 
% Trong truong hop khong co nhieu quan sat (SNRdB = vo 
% cuc), chon gia tri du lon (nhu SNRdB = 1000).

% Viet cho giao trinh: 
% Xu ly tin hieu ngau nhien, Dai hoc Quoc gia Ha Noi, 2024
% Tac gia: Nguyen Linh Trung, Huynh Huu Tue
% ========================================================
%%
% Chon mot trong hai tin hieu dieu hoa: co 1 tan so (hinh 
% sin) va co 3 tan so. Khoa cac dong lenh cua tin hieu 
% khong chon.
% A = 1;          % Tin hieu dieu hoa 1 tan so (hinh sin)
% fnu = .121;     
A = [1 1 1]';  % Tin hieu dieu hoa 3 tan so
fnu = [.1 .121 .2]';  
% Chon cac chieu dai tin hieu can danh gia
Lb_all = [100 500 1000 1500]; 

% Chon loai cua so
%kw = 1; % cua so Chu nhat
kw = 2; % cua so Hann
%kw = 3; % cua so Hamming
%kw = 4; % cua so Blackman
%kw = 5; % cua so Bartlett

% Tinh thong so nhieu
SNRdB = -10;
S = sum(A.^2)/2;
N0 = S*10^(-SNRdB/10);
sigma = sqrt(N0);
SNRstr = [num2str(SNRdB) ' dB'];
if SNRdB >= 1000
    SNRstr = '\infty';
end

%% Uoc luong pho bang Tuan hoan do va hien thi ket qua
figure
%Lbb = Lb_all(1:2); % chon 2 chieu dai trong Lb_all
Lbb = Lb_all(3:4); % chon 2 chieu dai trong Lb_all

for k = 1:length(Lbb) % cho 2 chieu dai khac nhau
    % Tao tin hieu quan sat co nhieu cong
    Lb = Lbb(k);   
    xsach = sum([diag(A)*cos(2*pi*fnu*(0:Lb-1) + 2*pi*rand(length(A),1)); zeros(1,Lb)]);
    x = xsach + sigma*randn(1,Lb);

    % Tinh du lieu cua so
    [Wd,U,txtcs] = cuaso(Lb,kw);
    
    % Tinh tuan hoan do
    nn = 0:Lb/2;
    nt = nn/Lb;
    Pxx = (abs(fft(x.*Wd).^2))/(U*Lb);
    PxxdB = 10*log10(Pxx);
    
    % Hien thi ket qua theo thang tuyen tinh
    subplot(length(Lbb),2,2*k-1);
    plot(nt,Pxx(1+nn));
    axis([0 .5 min(Pxx)-5 max(Pxx)+5])
    xlabel('Tan so chuan hoa, \nu')
    ylabel('Tuan hoan do, P_{xx}')
    legend(['Lb = ',num2str(Lb)])
    if length(A) == 1
        title({['Tuan hoan do Pxx cua tin hieu hinh sin co tan so']
            ['\nu = ', num2str(fnu(1)), '; SNR = ', SNRstr, '; cua so ', txtcs]})
    else
        title({['Tuan hoan do Pxx cua tin hieu dieu hoa co 3 tan so']
            ['\nu = ', num2str(fnu(1)), ', ', num2str(fnu(2)), ', ', num2str(fnu(3)), '; SNR = ', SNRstr, '; cua so ', txtcs]})
    end
    
    % Hien thi ket qua theo thang dB
    subplot(length(Lbb),2,2*k);
    plot(nt,PxxdB(1+nn));
    axis([0 .5 min(PxxdB)-10 max(PxxdB)+10])
    xlabel('Tan so chuan hoa, \nu')
    ylabel('Tuan hoan do, P_{xx} (dB)')
    legend(['Lb = ', num2str(Lb)])
    if length(A) == 1
        title({['Tuan hoan do Pxx (dB) cua tin hieu hinh sin co tan so']
            ['\nu = ', num2str(fnu(1)), '; SNR = ', SNRstr, '; cua so ', txtcs]})
    else
        title({['Tuan hoan do Pxx (dB) cua tin hieu dieu hoa co 3 tan so']
            ['\nu = ', num2str(fnu(1)), ', ', num2str(fnu(2)), ', ', num2str(fnu(3)), '; SNR = ', SNRstr, '; cua so ', txtcs]})
    end
end