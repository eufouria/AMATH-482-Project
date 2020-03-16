clear; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

utn = zeros(n, n, n);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    utn = utn + fftn(Un);
end

%averaging
utn_ave = abs(fftshift(utn./20));% averaging 
index = find(utn_ave == max(utn_ave(:)));%get indices
x_ss = Kx(index);
y_ss = Ky(index);
z_ss = Kz(index); 
    
% filter
tau = 0.5;
coord = zeros(20, 3);
filter = exp(-1*tau*((Kx - x_ss).^2 + (Ky - y_ss).^2 +(Kz - z_ss).^2));
filter = fftshift(filter);
for i=1:20
    Un(:,:,:)=reshape(Undata(i,:),n,n,n);
    utn = fftn(Un);   
    unft=filter.*utn;    
    unf = ifftn(unft);
    index = find(unf == max(unf(:)));
    xc = X(index);
    yc = Y(index);
    zc = Z(index);
    coord(i,:) = [xc yc zc];
end
plot3(coord(:,1), coord(:,2), coord(:,3), 'Linewidth', [2]), grid on; hold on;
plot3(coord(:,1), coord(:,2), coord(:,3), 'ro'); hold on;
plot3(coord(20,1), coord(20,2), coord(20,3), 'g*', 'Linewidth', [2])
xlim([-12, 12]);
ylim([-12, 12]);
zlim([-12, 12]);
xlabel('x');
ylabel('y');
zlabel('z');













