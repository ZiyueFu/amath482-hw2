% Amath482 HW2
% Fiona Fu
clc; clear; close all
% figure(1)
% [y, Fs] = audioread('GNR.m4a');
% tr_gnr = length(y)/Fs; % record time in seconds
% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Sweet Child O'' Mine');
% % p8 = audioplayer(y,Fs); playblocking(p8);
% 
% % GNR Guitar
% L = 14;  n = length(y);
% t2 = linspace(0,L,n+1); t = t2(1:n);
% k = (1/L)*[0:n/2-1 -n/2:-1];
% ks = fftshift(k);
% 
% % Spectrogram
% a = 20;
% tau = 0:0.1:14;
% 
% for j = 1:length(tau)
%     g = exp(-a*(t - tau(j)).^2); % Window function
%     Yg = g.*transpose(y); %Signal filter
%     Ygt = fft(Yg);
%     Ygt_spec(:,j) = fftshift(abs(Ygt));
% end
% figure(2)
% pcolor(tau,ks,log(abs(Ygt_spec)+ 1))
% shading interp
% set(gca,'ylim',[200 500],'Fontsize',12)
% colormap(hot)
% colorbar
% xlabel('Time (t)'), ylabel('Frequency (Hz)')
% title('Spectrogram of GNR - Guitar')



% Floyd Bass

figure(1)
[y1, Fs] = audioread('Floyd.m4a');
% y1 = y1((1:2635920),1);
y1 = y1((1:878640),1);
% y1 = y1((1757281:2635920),1);
tr_gnr = length(y1)/Fs; % record time in seconds
plot((1:length(y1))/Fs,y1);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Comfortably Numb');
% p8 = audioplayer(y1,Fs); playblocking(p8);

% FT of Floyd (divided into 3 clips, 20s each)
L = 20;  n = length(y1);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (1/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Spectrogram
a = 20;
tau = 0:0.1:20;

for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2); % Window function
    Yf = g.*transpose(y1); %Signal filter
    Yft = fft(Yf);
    Yft_spec(:,j) = fftshift(abs(Yft));
end

figure(2)
pcolor(tau,ks,log(abs(Yft_spec)+ 1))
shading interp
set(gca,'ylim',[150 1000],'Fontsize',12)
colormap(hot)
colorbar
xlabel('Time (t)'), ylabel('Frequency (Hz)')
title('Spectrogram of Floyd - Guitar(0-20s)')

% Floyd Bass

figure(1)
[y1, Fs] = audioread('Floyd.m4a');
% y1 = y1((1:2635920),1);
y1 = y1((1:878640),:);
% y1 = y1((1757281:2635920),1);
tr_gnr = length(y1)/Fs; % record time in seconds
plot((1:length(y1))/Fs,y1);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Comfortably Numb');
% p8 = audioplayer(y1,Fs); playblocking(p8);

% FT of Floyd (divided into 3 clips, 20s each)
L = 20;  n = length(y1);
t2 = linspace(0,L,n+1); t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

% Spectrogram
a = 5;
tau = 0:0.1:20;
freq_max = [];
for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2); % Window function
    Yf = g.*transpose(y1); %Signal filter
    Yft = fft(Yf);
    [f_max, ind] = max(Yft);
    freq_max = [freq_max;k(ind)];
    Yft_spec(:,j) = abs(fftshift(Yft));
end
figure(2)
plot(tau, abs(freq_max./(2*pi)))
xlabel('Time(s)'); ylabel('Frequency (Hz)');
title('Floyd-Bass isolated');

figure (3)
pcolor(tau, ks, Yft_spec)
shading interp
set(gca,'ylim',[150 1000],'Fontsize',12)
colormap(hot)
colorbar
xlabel('Time (t)'), ylabel('Frequency (Hz)')
title('Spectrogram of Floyd - Guitar(0-20s)')
