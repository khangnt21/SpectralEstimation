function [x,Sxx] = thdh(A,fnu,Lx,SNRdB)

%% TAO TIN HIEU DIEU HOA
%
% [x,Sxx] = thdh(A,fnu,Lx,SNRdB)
%
% Ham nay chon cac thong so de tao tin hieu dieu hoa chua 
% K tan so. A la vec-to cot chua K bien do cua K tan so 
% va fnu la vec-to cot chua K tan so cua tin hieu dieu 
% hoa. Ta co the thay doi A va fnu de co nhung tap du 
% lieu khac nhau. Dung quen la A va fnu phai co cung so 
% chieu K. Cac thanh phan cua tin hieu dieu hoa deu co 
% pha goc bat ky. Lx la chieu dai cua tin hieu. SNRdB la 
% ty le tin tren nhieu, theo don vi dB. Neu muon khong co
% tin hieu quan sat, chon SNRdB du lon (SNRdB >= 1000). 
% Dau ra la tin hieu dieu hoa co chua nhieu (x) va pho 
% cong suat ly tuong (Sxx) cua tin hieu dieu hoa khong
% co nhieu quan sat.

% Viet cho giao trinh: 
% Xu ly tin hieu ngau nhien, Dai hoc Quoc gia Ha Noi, 2024
% Tac gia: Nguyen Linh Trung, Huynh Huu Tue
% ========================================================

% Tao tin hieu sach (quan sat chinh xac)
xsach = sum(diag(A)*cos(2*pi*fnu*(0:Lx-1) + 2*pi*rand(length(A),1)));

% Tao tin hieu quan sat x, co sai so la nhieu cong
S = .5*sum(A.^2);           % Cong suat tin hieu
N0 = S*(10^(-.1*SNRdB));    % Cong suat nhieu
sigma = sqrt(N0);   
x = xsach + sigma*randn(1,Lx);

% Tinh Sxx cua {x(n)} su dung tin hieu dai Ns ma ta muon
Ns = Lx;
%Ns = min(Lx,2000); % mo lenh nay neu muon han che so 
                    % luong mau pho
Sxx = (abs(fft(x(1:Ns))).^2)/Ns;    % tinh tuan hoan do