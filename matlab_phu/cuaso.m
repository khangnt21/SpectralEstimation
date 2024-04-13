function [Wd,U,txtcs] = cuaso(Lw,kw)

%% TAO CUA SO THONG DUNG
%
% [Wd,U,txtcs] = cuaso(Lw,kw)
%
% Thong so kw dung de xac dinh loai cua so nao duoc tao 
% ra va Lw la chieu dai tuong ung cua cua so do. Co 5 cua 
% so thong dung va kw tuong ung la Chu nhat (kw = 1), 
% Hann (kw = 2), Hamming (kw = 3), Blackmann (kw = 4) va 
% Bartlett (kw = 5). Dau ra la ten o dang text cua cua so 
% duoc chon txtcs, chuoi gia tri Wd cua cua so duoc chon, 
% he so dieu chinh U. 

% Viet cho giao trinh: 
% Xu ly tin hieu ngau nhien, Dai hoc Quoc gia Ha Noi, 2024
% Tac gia: Nguyen Linh Trung, Huynh Huu Tue
% ========================================================
%%
N  = Lw-1;
nm = 0:N;
M  = N-1;

% Tao ma tran Mw chua gia tri 5 cua so thong dung
Mw = [ones(1,N+1);  % Chu nhat
    0.5*(1-cos(2*pi*nm/N)); % Hann
    0.54 - 0.46*cos(2*pi*nm/N);   % Hamming    
    0.42 - 0.5*cos(2*pi*nm/N) + 0.08*cos(4*pi*nm/N); % Blackman
    (2/N)*[0:N/2,(floor(N/2)-(0.5*N <=floor(0.5*N))):-1:0]]; % Bartlett
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

Wd = Mw(kw,:);      % chuoi gia tri cua so thu kw
U = mean(Wd.^2);    % he so dieu chinh cua cua so