%% load data
clc
clear all
close all
T = load('spikes_gratings/data_monkey3_gratings');
TALL = load('spikes_gratings/S_monkey3');
TSOME = load('spikes_gratings/Snew_monkey3');
Unum = [106 88 112];
%% set parameters
monkeyNum = 3;
newNum = size(TSOME.S(1).mean_FRs,1);
stimNum = 12;
unitNum = Unum(monkeyNum);
%%
Counts = zeros(stimNum, newNum, 200);
for k = 1 : 12
    for j = 1 : 200
        temp = TSOME.S(k).trial(j).counts;
        for i = 1 : newNum
            Counts(k, i, j) = sum(temp(i, :));
        end
    end
end
%%
Vectors = zeros(newNum, stimNum*200);
for i = 1 : newNum
    for j = 1 : stimNum
        Vectors(i, (j-1)*200+1 : j*200) = (Counts(j, i, :)-mean(Counts(j, i, :)))/sqrt(var(Counts(j, i, :)));
    end
end
%% find index in all neurons
idxInAll = zeros(newNum, 1);
t2 = TALL.S(1).mean_FRs;
t1 = TSOME.S(1).mean_FRs;
for j =  1 : newNum
    for i = 1 : unitNum
        if(isequal(t1(j, 1:50), t2(i, 1:50)))
            idxInAll(j, 1) = i;
        end
    end
end
%%
mainMAP = T.data.MAP;
mainChannels = T.data.CHANNELS;
data = zeros(newNum*(newNum-1)/2, 3);
cnt = 0;
for i = 1 : newNum
    for j = i+1 : newNum
        cnt = cnt + 1;
        %rsc
        rscTemp = corrcoef(Vectors(i, :), Vectors(j, :));
        rsc = rscTemp(2, 1);
        %rsig
        TPlotMat1 = zeros(1, stimNum);
        for k = 1: stimNum
            temp = mean(TALL.S(k).mean_FRs,2);
            TPlotMat1(k) = temp(idxInAll(i)); 
        end
        TPlotMat2 = zeros(1, stimNum);
        for k = 1: stimNum
            temp = mean(TALL.S(k).mean_FRs,2);
            TPlotMat2(k) = temp(idxInAll(j)); 
        end
        rsigTemp = corrcoef(TPlotMat1, TPlotMat2);
        rsig = rsigTemp(2, 1);
        %distance
        [d1, d2] = find(mainMAP == mainChannels(idxInAll(j), 1));
        [d3, d4] = find(mainMAP == mainChannels(idxInAll(i), 1));
        distance = 0.4*sqrt((d1-d3)^2+(d2-d4)^2);
        data(cnt, :) = [rsc rsig distance];
    end
end
%% first plot
close all
L = find(data(:, 2)>=0.5);
tempData = data(L, :);
numbers = zeros(1,8);
stds = zeros(1,8);
for i = 1 : 8
    Distance1 = find(tempData(:, 3) > 0.5*(i-1) + 0.25 & tempData(:, 3) <= 0.5*i + 0.25);
    numbers(i) = mean(tempData(Distance1, 1));
    stds(i) = sqrt(var(tempData(Distance1, 1)))/sqrt(length(Distance1));
end
% plot(0.5:0.5:4, numbers)
errorbar(0.5:0.5:4, numbers, stds)
L = find(data(:, 2)<0.5 & data(:, 2)>=0);
tempData = data(L, :);
for i = 1 : 8
    Distance1 = find(tempData(:, 3) > 0.5*(i-1) + 0.25 & tempData(:, 3) <= 0.5*i + 0.25);
    numbers(i) = mean(tempData(Distance1, 1));
    stds(i) = sqrt(var(tempData(Distance1, 1)))/sqrt(length(Distance1));
end
hold on
% plot(0.5:0.5:4, numbers)
errorbar(0.5:0.5:4, numbers, stds)

L = find(data(:, 2)<0 & data(:, 2)>=-0.5);
tempData = data(L, :);
for i = 1 : 8
    Distance1 = find(tempData(:, 3) > 0.5*(i-1) + 0.25 & tempData(:, 3) <= 0.5*i + 0.25);
    numbers(i) = mean(tempData(Distance1, 1));
    stds(i) = sqrt(var(tempData(Distance1, 1)))/sqrt(length(Distance1));
end
hold on
% plot(0.5:0.5:4, numbers)
errorbar(0.5:0.5:4, numbers, stds)
L = find(data(:, 2)<-0.5);
tempData = data(L, :);

for i = 1 : 8
    Distance1 = find(tempData(:, 3) > 0.5*(i-1) + 0.25 & tempData(:, 3) <= 0.5*i + 0.25);
    numbers(i) = mean(tempData(Distance1, 1));
    stds(i) = sqrt(var(tempData(Distance1, 1)))/sqrt(length(Distance1));
end
hold on
% plot(0.5:0.5:4, numbers)
errorbar(0.5:0.5:4, numbers, stds)

legend("$r_{sig}\ge0.5$","$0.5>r_{sig}\ge0$", "$0>r_{sig}\ge-0.5$","$-0.5>r_{sig}$", 'interpreter', 'latex')
xlabel("Distance between electrodes (mm)", 'interpreter', 'latex')
ylabel("Spike count correlation ($r_{sc}$)",'interpreter', 'latex')
title(" Dependence of $r_{sc}$ on distance for monkey 3",'interpreter', 'latex')


%% second plot
figure
L = find(data(:, 3)>0 & data(:, 3)<=1);
tempData = data(L, :);
numbers = zeros(1,10);
stds = zeros(1,10);

for i = 1 : 10
    rsig1 = find(tempData(:, 2) >= 0.2*(i-1)-1 & tempData(:, 2) < 0.2*(i)-1);
    numbers(i) = mean(tempData(rsig1, 1));
    stds(i) = sqrt(var(tempData(rsig1, 1)))/sqrt(length(rsig1));
end
errorbar(-0.9:0.2:0.9, numbers, stds)
hold on
L = find(data(:, 3)>1 & data(:, 3)<=2);
tempData = data(L, :);
for i = 1 : 10
    rsig1 = find(tempData(:, 2) >= 0.2*(i-1)-1 & tempData(:, 2) < 0.2*(i)-1);
    numbers(i) = mean(tempData(rsig1, 1));
    stds(i) = sqrt(var(tempData(rsig1, 1)))/sqrt(length(rsig1));
end
errorbar(-0.9:0.2:0.9, numbers, stds)
hold on

L = find(data(:, 3)>2 & data(:, 3)<=3);
tempData = data(L, :);
for i = 1 : 10
    rsig1 = find(tempData(:, 2) >= 0.2*(i-1)-1 & tempData(:, 2) < 0.2*(i)-1);
    numbers(i) = mean(tempData(rsig1, 1));
    stds(i) = sqrt(var(tempData(rsig1, 1)))/sqrt(length(rsig1));
end
errorbar(-0.9:0.2:0.9, numbers, stds)
hold on
L = find(data(:, 3)>3 & data(:, 3)<=4);
tempData = data(L, :);
for i = 1 : 10
    rsig1 = find(tempData(:, 2) >= 0.2*(i-1)-1 & tempData(:, 2) < 0.2*(i)-1);
    numbers(i) = mean(tempData(rsig1, 1));
    stds(i) = sqrt(var(tempData(rsig1, 1)))/sqrt(length(rsig1));
end
errorbar(-0.9:0.2:0.9, numbers, stds)
hold on
L = find(data(:, 3)>4);
tempData = data(L, :);
for i = 1 : 10
    rsig1 = find(tempData(:, 2) >= 0.2*(i-1)-1 & tempData(:, 2) < 0.2*(i)-1);
    numbers(i) = mean(tempData(rsig1, 1));
    stds(i) = sqrt(var(tempData(rsig1, 1)))/sqrt(length(rsig1));
end
errorbar(-0.9:0.2:0.9, numbers, stds)

legend("$1>distance\ge0$","$2>distance\ge1$","$3>distance\ge2$", "$4>distance\ge3$", "$distance\ge4$", 'Location','north','interpreter', 'latex')
xlabel("Orientation tuning similarity ($r_{sig}$)",'interpreter', 'latex')
ylabel("Spike count correlation ($r_{sc}$)",'interpreter', 'latex')
title(" Dependence of $r_{sc}$ on $r_{signal}$ for monkey 3",'interpreter', 'latex')



%%
close all
figure
NewMAP = zeros(7, 10);
for i = 1 : 7
    for j = 1 : 10
        rr = find(data(:, 3) > 0.5*(i-1) + 0.25 & data(:, 3) <= 0.5*i + 0.25 & data(:, 2) >= 0.2*(j-1)-1 & data(:, 2) < 0.2*(j)-1);
        NewMAP(i, j) = mean(data(rr, 1));
    end
end

im = image([0.5 3.5],[-0.9 0.9],NewMAP(:, :)');
im.CDataMapping = 'scaled';
c = colorbar;
c.Label.String = sprintf('Spike count correlation');
set(gca,'YDir','normal')
xlabel("Distance between electrodes (mm)", 'interpreter', 'latex')
ylabel("Orientation tuning similarity ($r_{sig}$)",'interpreter', 'latex')
title(" Dependence of $r_{sc}$ on both tuning similarity ($r_{signal}$) and distance for monkey 3",'interpreter', 'latex')






