clear all; close all; clc
load handel
v = y';
v1 = v;
v = v(1:length(v) - 1);
%p8 = audioplayer(v,Fs);playblocking(p8);
L=9; n = length(v);
t2=linspace(0,L,n+1); t=t2(1:n); 
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);
tslide=0:0.1:9;
%% comparing 3 windows
j = 45
tslide = 0:0.1:9;
gaussian
g=exp(-10*(t-tslide(j)).^2);
mexican hat
mxh = (1 - (t-tslide(j)).^2).*exp(-(t-tslide(j)).^2/2); 
shannon filter
s = zeros(1, length(t));
s = zeros(1, length(t));
index_of_tslidej = round(tslide(45)/9*length(t)); % index of tslide
index_width = round(0.5/9*length(t));
if index_of_tslidej-index_width <= 0 
    s(1:index_of_tslidej+index_width) = ones();
elseif index_of_tslidej+index_width > length(t);
    s(index_of_tslidej-index_width:length(t)) = ones();
else
    s(index_of_tslidej-index_width:index_of_tslidej+index_width) = ones();
end;
plot(t, v); hold on 
a1 = plot(t, mxh, 'linewidth',[2]); L1 = 'Mexican Hat';
a2 = plot(t, g,  'linewidth',[2]); L2 = 'Gaussian';
a3 = plot(t, s,  'linewidth',[2]); L3 = 'Shannon';
ylim([-1.5 1.5])
ylabel('Amplitude')
xlabel('Time [sec]')
legend([a1;a2;a3], L1, L2, L3)

%% exploring changing width
a_vector = [0.1 10 100 10000]; %different width
tslide=0:0.1:9;
Sgt_spec = zeros(length(tslide),n);
for i=1:length(a_vector)
   a = a_vector(i);
   for j=1:length(tslide)
       g=exp(-a.*(t-tslide(j)).^2); 
       vfilter = g.*v; 
       vfiltert= fft(vfilter); 
       Sgt_spec(j,:) = fftshift(abs(vfiltert)); % We don't want to scale it
   end
   
   subplot(2, 2, i)
   pcolor(tslide,ks/(2*pi),Sgt_spec.'), 
   xlabel('Time [sec]')
   ylabel('Frequency')
   shading interp 
   colormap(hot)
   title(['Gabor Transform translation = 0.1, width = ', num2str(a)])
   saveas(gca, ['gabor_trans.png'])
end

%% exploring change translation 
t_step = [0.005]; % large step mean undersampling, small step mean oversampling
for i=1:length(t_step)
a = 100;
step = t_step(i);
tslide=0:step:9;
Sgt_spec = zeros(length(tslide),n);
for j=1:length(tslide)
g=exp(-a.*(t-tslide(j)).^2); 
vfilter = g.*v; 
vfiltert= fft(vfilter); 
Sgt_spec(j,:) = fftshift(abs(vfiltert)); % We don't want to scale it
end
pcolor(tslide,ks/(2*pi),Sgt_spec.')
xlabel('Time [sec]')
ylabel('Frequency')
shading interp 
colormap(hot)
title(['Gabor Transform width = 100, translation =',num2str(step)])
saveas(gcf, ['translation=',num2str(step),'.png'])
end



%% different window: mexican hat && shanon
mexican hat


a = 100;
tslide=0:0.1:9;
for j=1:length(tslide)
    filter = (1 - (t-tslide(j)).^2).*exp(-a*(t-tslide(j)).^2/2);
    vfilter = v.*filter;
    vfiltert= fft(vfilter); 
    Sgt_spec(j,:) = fftshift(abs(vfiltert)); 
end
pcolor(tslide,ks/(2*pi),Sgt_spec.'), 
shading interp 
colormap(hot)
title(['Gabor Transform with Mexican Hat Wavelet a = 100 and translation = 0.1'])
saveas(gcf, ['mexican_hat.png'])

%% shannon filter
Sgt_spec = zeros(length(tslide),n);
for j=1:length(tslide)
    s = zeros(1, length(t));
    index_of_tslidej = round(tslide(j)/9*length(t)); % index of tslide
    index_width = round(0.5/9*length(t));
    if index_of_tslidej-index_width <= 0 && index_of_tslidej+index_width > length(t)
        s(1:length(t)) = ones();
    elseif index_of_tslidej-index_width <= 0 
        s(1:index_of_tslidej+index_width) = ones();
    elseif index_of_tslidej+index_width > length(t)
        s(index_of_tslidej-index_width:length(t)) = ones();
    else
        s(index_of_tslidej-index_width:index_of_tslidej+index_width) = ones();
    end
    filter = abs(t-tslide(j)) < 1/2;
    vfilter = v.*s;
    vfiltert= fft(vfilter); 
    Sgt_spec(j,:) = fftshift(abs(vfiltert)); 
end
pcolor(tslide,ks/(2*pi),Sgt_spec.'), 
shading interp 
colormap(hot)
title(['Gabor Transform with Shannon'])
saveas(gca, 'shannon.png')

    
   



