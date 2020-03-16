clear all; close all; clc
%% data processing
%chad crouch
inst1 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Lobo_Loco/Not_my_Brain/Lobo_Loco_-_01_-_Brain_ID_1270.mp3";
inst2 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/My_Yearnings/Sometimes_it_Rains_ID_1206.mp3";
inst3 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_I_Oaks_Bottom_Instrumental/Chad_Crouch_-_Drumming_In_The_Rain_Instrumental.mp3";
inst4 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_I_Oaks_Bottom_Instrumental/Chad_Crouch_-_Ballad_Of_The_Blackbirds_Instrumental.mp3";
inst5 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_I_Oaks_Bottom_Instrumental/Chad_Crouch_-_The_Chorus_Ceases_Instrumental.mp3";
inst6 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_II_Reed_Canyon_Instrumental/Chad_Crouch_-_01_-_Headwaters_Instrumental.mp3";
inst7 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_II_Reed_Canyon_Instrumental/Chad_Crouch_-_Arrival_Of_The_Geese_Instrumental.mp3";
inst8 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_II_Reed_Canyon_Instrumental/Chad_Crouch_-_Children_By_The_Creek.mp3";
inst9 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_II_Reed_Canyon_Instrumental/Chad_Crouch_-_The_Family_Instrumental.mp3";
inst10 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_II_Reed_Canyon_Instrumental/Chad_Crouch_-_05_-_Song_Sparrow_Serenade_Instrumental.mp3";
inst11 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/ccCommunity/Chad_Crouch/Field_Report_Vol_I_Oaks_Bottom_Instrumental/Chad_Crouch_-_04_-_The_Spring_Instrumental.mp3";

%ketsa
elec1 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_12_-_Green_Man.mp3";
elec2 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_11_-_Slow_Vibing.mp3";
elec3 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_07_-_The_Stork.mp3";
elec4 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequency/Ketsa_-_08_-_Multiverse.mp3";
elec5 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_13_-_Mission_Ready.mp3";
elec6 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_10_-_Memories_Renewed.mp3";
elec7 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_09_-_Life_Illusion.mp3";
elec8 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_04_-_Enclosed.mp3";
elec9 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_03_-_Dusty_Hills.mp3";
elec10 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_01_-_Dreaming_Days.mp3";
elec11 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_06_-_Day_Trips.mp3";
elec12 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/Creative_Commons/Ketsa/Raising_Frequecy/Ketsa_-_05_-_Crescents.mp3";

%Yung Katz
hhop1 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/August_2019/Yung_Kartz_-_04_-_One_Way.mp3";
hhop2 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/August_2019/Yung_Kartz_-_03_-_Intranet.mp3";
hhop3 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/August_2019/Yung_Kartz_-_02_-_Stranger.mp3";
hhop4 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/August_2019/Yung_Kartz_-_01_-_Jeopardy.mp3";
hhop5 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/June_2019/Yung_Kartz_-_08_-_Too_Much_Ice.mp3";
hhop6 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/June_2019/Yung_Kartz_-_06_-_Too_Grimy.mp3";
hhop7 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/June_2019/Yung_Kartz_-_05_-_Starz.mp3";
hhop8 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/June_2019/Yung_Kartz_-_07_-_Straight_Shot.mp3";
hhop9 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/June_2019/Yung_Kartz_-_02_-_Young_and_Reckless.mp3";
hhop10 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/July_2019/Yung_Kartz_-_07_-_Frontline.mp3";
hhop11 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/July_2019/Yung_Kartz_-_08_-_Landline.mp3";
hhop12 = "https://files.freemusicarchive.org/storage-freemusicarchive-org/music/no_curator/Yung_Kartz/July_2019/Yung_Kartz_-_05_-_Magic.mp3";

instrument = [inst1 inst2 inst3 inst4 inst5 inst6 inst7 inst8 inst9 inst10 inst11];
electric = [elec1 elec2 elec3 elec4 elec5 elec6 elec7 elec8 elec9 elec10 elec11];
hiphop = [hhop1 hhop2 hhop3 hhop4 hhop5 hhop6 hhop7 hhop8 hhop9 hhop10 hhop11];
data1 = [];
data2 = [];
data3 = [];
for i = 1:11
    [y, Fs] = webread(electric(i));
    [m,n]=size(y);
    dt=1/Fs;
    t=dt*(0:m-1);
    idx1 = (t>10) & (t<15);
    y1 = y(idx1);
    data1 = [data1; y1];
end
for i = 1:11
    [y, Fs] = webread(instrument(i));
    [m,n]=size(y);
    dt=1/Fs;
    t=dt*(0:m-1);
    idx1 = (t>10) & (t<15);
    y1 = y(idx1);
    data2 = [data2; y1];
end
for i = 1:11
    [y, Fs] = webread(hiphop(i));
    [m,n]=size(y);
    dt=1/Fs;
    t=dt*(0:m-1);
    idx1 = (t>10) & (t<15);
    y1 = y(idx1);
    data3 = [data3; y1];
end
%% classification

feature = 33;
[U,S,V,w,sortelec,sortinst, sorthhop] = dc_trainer(data1',data2',data3', feature);

plot(sortelec,zeros(11),'dr','Linewidth',2)
hold on
plot(sortinst,ones(11),'ob','Linewidth',2)
hold on
twos = zeros(11) + 2;
plot(sorthhop,twos,'*g','Linewidth',2)
ylim([0 2])

threshold1 = (max(sortelec) + max(sortinst))/2;

t1 = length(sortinst);
t2 = 1;
while sortinst(t1)>sorthhop(t2)
    t1 = t1-1;
    t2 = t2+1;
end
threshold2 = (sortinst(t1)+sorthhop(t2))/2;

for i = 1:11
    [y, Fs] = webread(electric(i));
    [m,n]=size(y);
    dt=1/Fs;
    t=dt*(0:m-1);
    idx = (t>10) & (t<15);
    y = y(idx);
    TestMat = U'*y';  % PCA projection
    pval = w'*TestMat; % LDA projection

    if pval < threshold1
        disp('elec')
    elseif pval > threshold2
        disp('hhop')
    else
        disp('inst')
    end
end



%% dc_trainer


function dcData = dc_wavelet(dcfile)
    [m,n] = size(dcfile); 
    pxl = sqrt(m);
    nw = m/4; % wavelet resolution
    dcData = zeros(nw,n);
    
    for k = 1:n
       X = im2double(reshape(dcfile(:,k),pxl,pxl)); 
       [~,cH,cV,~]=dwt2(X,'haar');
       cod_cH1 = rescale(abs(cH));
       cod_cV1 = rescale(abs(cV));
       cod_edge = cod_cH1+cod_cV1;
       dcData(:,k) = reshape(cod_edge,nw,1);
    end
end

function [U,S,V,w, sortelec,sortinst, sorthhop] = dc_trainer(elec,inst,hhop,feature)
    
    [U,S,V] = svd([elec inst hhop],'econ');
    
    music = S*V'; % projection onto principal components
    U = U(:,1:feature);
    elecs = music(1:feature,1:11);
    insts = music(1:feature,12:22);
    hhops = music(1:feature,23:33);
    
    m1 = mean(elecs,2);
    m2 = mean(insts,2);
    m3 = mean(hhops,2);
    m = (m1 + m2 + m3)/3;
    
    Sw = 0; % within class variances
    for k=1:11
        Sw = Sw + (elecs(:,k)-m1)*(elecs(:,k)-m1)';
    end
    for k=1:11
        Sw = Sw + (insts(:,k)-m2)*(insts(:,k)-m2)';
    end
    for k=1:11
        Sw = Sw + (hhops(:,k)-m3)*(hhops(:,k)-m3)';
    end
    
    Sb = (m1-m)*(m1-m)' + (m2-m)*(m2-m)' + (m2-m)*(m1-m)'; % between class 
    
    [V2,D] = eig(Sb,Sw); % linear discriminant analysis
    [~,ind] = max(abs(diag(D)));
    w = V2(:,ind); w = w/norm(w,2);
    
    
    velec = w'*elecs;
    vinst = w'*insts; 
    vhhop = w'*hhops; 
        
    sortinst = sort(vinst);
    sortelec = sort(velec);
    sorthhop = sort(vhhop);
   
end

