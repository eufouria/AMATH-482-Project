clear all; close all; clc


load('test_data.mat');
load('train_data.mat');


%% prepare testset

feature = 42;
[U,S,V,w,reddy,bichphuong,bray] = dc_trainer(data1,data2,data3, feature);

mean(reddy);
mean(bichphuong);
mean(bray);


TestMat = U'*TestSet;  % PCA projection
pval = w'*TestMat; % LDA projection

figure(2)
for k=1:3
    subplot(3,3,3*k-2)
    plot(1:24,V(1:24,k),'ko-'), title('Reddy')
    subplot(3,3,3*k-1)
    plot(25:48,V(25:48,k),'ko-'), title('BichPhuong')
    subplot(3,3,3*k)
    plot(49:72,V(49:72,k),'ko-'), title('bray') 
end

t1 = length(bray);
t2 = 1;
while bray(t1)>bichphuong(t2)
    t1 = t1-1;
    t2 = t2+1;
end
threshold1 = (bray(t1)+bichphuong(t2))/2
t1 = length(bichphuong);
t2 = 1;
while bichphuong(t1)>reddy(t2)
    t1 = t1-1;
    t2 = t2+1;
end
threshold2 = (bichphuong(t1)+reddy(t2))/2
count = 0;
for i = 1:size(pval, 2)
    if pval(i) < threshold1
        if label(i) == 2
            count = count + 1;
        end
    elseif pval(i) > threshold2
        if label(i) == 0
            count = count + 1;
        end
    else
        if label(i) == 1
            count = count + 1;
        end
    end
end
accuracy_rate = count/size(label,2)

figure(1)
twos = zeros(1, size(bray, 2)) + 2;
%h = plot(reddy, zeros(size(reddy, 2)),'dr',noo, ones(size(noo, 2)),'ob',onlyc, twos, '*g','Linewidth',2)

scatter(reddy, zeros(1, size(reddy, 2)),'d', 'Linewidth', 2)
hold on
scatter(bichphuong, ones(1, size(bichphuong,2)),'o', 'Linewidth', 2)
scatter(bray, twos, '*', 'Linewidth', 2)
xline(threshold1,'--r');
xline(threshold2,'--k');
legend('Reddy', 'Bich Phuong', 'Bray','Threshold 1', 'Threshold 2', 'Location','northwest')





%%MUSIC_TRAINER
function [U,S,V,w,sortelec,sortinst, sorthhop] = dc_trainer(elec,inst,hhop,feature)
    [m1, n1] = size(elec);
    [m3, n2] = size(inst);
    [m3, n3] = size(hhop);
    [U,S,V] = svd([elec inst hhop],'econ');
    
    music = S*V'; % projection onto principal components
    U = U(:,1:feature);
    elecs = music(1:feature,1:n1);
    insts = music(1:feature,n1+1:n1+n2);
    hhops = music(1:feature,n1+n2+1:n1+n2+n3);
    
    m1 = mean(elecs,2);
    m2 = mean(insts,2);
    m3 = mean(hhops,2);
    m = (m1 + m2 + m3)/3;
    
    Sw = 0; % within class variances
    for k=1:n1
        Sw = Sw + (elecs(:,k)-m1)*(elecs(:,k)-m1)';
    end
    for k=1:n2
        Sw = Sw + (insts(:,k)-m2)*(insts(:,k)-m2)';
    end
    for k=1:n3
        Sw = Sw + (hhops(:,k)-m3)*(hhops(:,k)-m3)';
    end
    
    Sb = n1*(m1-m)*(m1-m)' + n2*(m2-m)*(m2-m)' + n3*(m3-m)*(m3-m)'; % between class 
    
    [V2,D] = eig(Sb,Sw); % linear discriminant analysis
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);
    
    velec = w'*elecs;
    vinst = w'*insts; 
    vhhop = w'*hhops; 
    
    sortelec = sort(velec);
    sortinst = sort(vinst);
    sorthhop = sort(vhhop);
   
end

