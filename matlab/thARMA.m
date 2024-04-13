function [x,Sxx] = thARMA(b,a,Lx,SNRdB)

%% TAO TIN HIEU ARMA
%
% [x,Sxx] = thARMA(b,a,Lx,SNRdB)
%
% Chon cac thong so de tao tin hieu ARMA chua Q tri khong
% va P tri cuc. Ta co the thay doi 2 vec-to hang 
% b = [b0 b1 ... bQ] va a = [1  a1 ... aP] 
% cua ham chuyen
%          b0 + b1*z^(-1) + ... + bQ*z^(-Q)
%   H(z) = --------------------------------
%          1  + a1*z^(-1) + ... + aP*z^(-P)
% ma ta dung de tao nhung tap du lieu khac nhau, bang 
% cach kich thich H(z) voi mot nhieu trang o dau vao. Dau
% ra cua H chinh la tin hieu sach ARMA ta muon tao ra. Lx
% la chieu dai cua tin hieu va SNRdB la ty le tin tren 
% nhieu theo don vi (dB).

% Viet cho giao trinh: 
% Xu ly tin hieu ngau nhien, Dai hoc Quoc gia Ha Noi, 2024
% Tac gia: Nguyen Linh Trung, Huynh Huu Tue
% ========================================================

% Tao tin hieu sach (quan sat chinh xac)
xsach = filter(b,a,randn(1,Lx));   

% Tao tin hieu quan sat, co sai so la nhieu cong
S = mean(xsach.^2);             % Cong suat tin hieu
N0 = S*(10^(-.1*SNRdB));        % Cong suat nhieu
sigma = sqrt(N0);               
x = xsach + sigma*randn(1,Lx); 

% Tinh Sxx ly tuong tai Lf tan so chuan hoa (nu) tu 0 
% den 0.5
Lf = 250;   % so diem tan so so goc Wk = k*pi/Lf, 
            % k = 0,...,Lf-1
H = freqz(b,a,Lf);  % Dap ung tan so H(W) tai cac diem 
                    % nu = Wk/(2*pi) 
Sxx = abs(H).^2;