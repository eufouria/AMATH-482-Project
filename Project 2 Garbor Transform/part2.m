clear all; close all; clc

piano-----------------------------------------------------------
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs;  % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)');
p8 = audioplayer(y,Fs); playblocking(p8);
n = length(y);
t2=linspace(0,tr_piano,n+1); t=t2(2:n+1); 
k=(2*pi/tr_piano)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);
v1 = y';
tslide = 0:0.25:tr_piano;

Sgt_spec = zeros(length(tslide),n);
music = [];
for j=1:length(tslide)
   g=exp(-40*(t-tslide(j)).^2); 
   v1filter = g.*v1;
   v1filtert= fft(v1filter); 
   Sgt_spec(j,:) = abs(fftshift(v1filtert)); 
       
   index = find(v1filtert == max(v1filtert(:)));
   music = [music, abs(k(index))/(2*pi)];
end
pcolor(tslide,ks/(2*pi),abs(Sgt_spec).'),
shading interp
colormap(hot)
ylim([0 500])

title(['Spectrogram of Piano sound'])
plot(tslide, music, 'Linewidth', 2), hold on

plot((1:length(y))/Fs,y*1000)
ylim([0 400])
xlabel('Time[sec]')
ylabel('Frequency')
title('Frequency of the Piano')
321 288 256 288 319 320 320 284 284 284 320 320 320 319 288 258 288
321 320 323 319 287 287 320 287 254
E4 D4 C4 D4 E4 E4 E4 D4 D4 D4 E4 E4 E4 E4 D4 C4 D4 E4 E4 E4 E4 D4 D4
E4 D4 C4


%recording----------------------------------------------------------------------------


[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs;  % record time in seconds
% plot((1:length(y))/Fs,y);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Mary had a little lamb (recorder)');
%p8 = audioplayer(y,Fs); playblocking(p8);
n = length(y);
t2=linspace(0,tr_rec,n+1); t=t2(2:n+1); 
k=(2*pi/tr_rec)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);
v2 = y';
tslide=0:0.25:tr_rec;
Sgt_spec = zeros(length(tslide),n);
music = [];
for j=1:length(tslide)
   g=exp(-40.*(t-tslide(j)).^2); 
   v2filter = v2.*g;
   v2filtert= fft(v2filter); 
   Sgt_spec(j,:) = fftshift(abs(v2filtert)); 
       
   index = find(v2filtert == max(v2filtert(:)));
   music = [music, abs(k(index))/(2*pi)];
end
% pcolor(tslide,ks/(2*pi),abs(Sgt_spec).'), 
% shading interp 
% colormap(hot)
% ylim([0 2000])
% title(['Record Spectrogram '])
plot(tslide, music, 'Linewidth', 2), hold on
plot((1:length(y))/Fs,y*4000)
xlabel('Time')
ylabel('Frequency')
title('Frequency of the Recorder')
% % 1050 901 804 908 1005 1006 1007 900 900 900 1052 1052 1052 1052 900 806 894
% % 1036 1036 1036 1036 849 911 994 845 794
% % B A G A B B B A A A B B  B B A G A B B B B  A A B A G 


