clear all; close all; clc
load('cam1_2.mat')
load('cam2_2.mat')
load('cam3_2.mat')

for i  = 1:size(vidFrames1_2, 4)
    X = double(rgb2gray(vidFrames1_2(:,:,:,i)));
    X_filter = X(200:480,300:400);
    flash = max(X_filter(:));
    [x,y] = find(X_filter >= flash*9/10);
    x1(i)=mean(y);
    y1(i)=mean(x);
end

for i  = 1:size(vidFrames2_2, 4)
    X = double(rgb2gray(vidFrames2_2(:,:,:,i))); 
    X_filter = X(50:480,200:400);
    flash = max(X_filter(:));
    [x,y] = find(X_filter >= flash*9/10);
    x2(i)=mean(y);
    y2(i)=mean(x);
end  

for i = 1:size(vidFrames3_2, 4)
    X = double(rgb2gray(vidFrames3_2(:,:,:,i)));
    X_filter = X(200:350, 250:500);
    flash = max(X_filter(:));
    [x,y] = find(X_filter >= flash*9/10);
    x3(i)=mean(y);
    y3(i)=mean(x);
end

index1 = find(y1(1:30) == min(y1(1:30)));
index2 = find(y2(1:30) == min(y2(1:30)));
index3 = find(y3(1:30) == min(y3(1:30)));
x1_shift = x1(index1:end);
y1_shift = y1(index1:end);
x2_shift = x2(index2:end);
y2_shift = y2(index2:end);
x3_shift = x3(index3:end);
y3_shift = y3(index3:end);

min_frame = min([max(size(x1_shift)) max(size(x2_shift)) max(size(x3_shift))]);
t=1:min_frame;
subplot(3,1,1)
plot(t,x1_shift(1:min_frame),t,y1_shift(1:min_frame))
xlabel('Frame'); ylabel('Location')
legend('Horizontal','Vertical','Location','NorthEast')
title('camera 1')

subplot(3,1,2)
plot(t,x2_shift(1:min_frame),t,y2_shift(1:min_frame))
xlabel('Frame'); ylabel('Location')
legend('Horizontal','Vertical','Location','NorthEast')
title('camera 2')

subplot(3,1,3)
plot(t,x3_shift(1:min_frame),t,y3_shift(1:min_frame))
xlabel('Frame'); ylabel('Location')
legend('Horizontal','Vertical','Location','NorthEast')
title('camera 3')


data = [x1_shift(1:min_frame); y1_shift(1:min_frame); x2_shift(1:min_frame); y2_shift(1:min_frame); x3_shift(1:min_frame); y3_shift(1:min_frame)];
[m, n] = size(data);
mn = mean(data, 2);
data = data-repmat(mn,1,n);
[u, s, v] = svd(data/sqrt(n-1), 'econ');
figure(1)
plot(diag(s).^2/sum(diag(s.^2)), 'o', 'Linewidth', 2)
xlabel('Principal component'); ylabel('Percentage of Energy');
ylim([0 1])
Y=u'*data; % produce the principal components projection
figure(2)
plot(Y(1,:)); hold on
plot(Y(2,:))
plot(Y(3,:))
plot(Y(4,:))
plot(Y(5,:))
plot(Y(6,:))
legend('Component 1','Component 2','Component 3','Component 4','Component 5','Component 6','Location','NorthEast')
xlabel('Frame'); 



