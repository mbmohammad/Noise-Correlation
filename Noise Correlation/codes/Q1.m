%% load data
clc
clear all
close all
TALL = load('spikes_gratings/S_monkey1');
TSOME = load('spikes_gratings/Snew_monkey1');
Unum = [106 88 112];
%% set parameters
monkeyNum = 1;
newNum = size(TSOME.S(1).mean_FRs,1);
stimNum = 12;
unitNum = Unum(monkeyNum);
%% most active neuron for each stimulus
MFRs = zeros(newNum, stimNum);
for i = 1 : stimNum
    temp = TSOME.S(i).mean_FRs;
    MFRs(:, i) = mean(temp,2);
end
%%
MFRs1 = zeros(unitNum, stimNum);
for i = 1 : stimNum
    temp = TALL.S(i).mean_FRs;
    MFRs1(:, i) = mean(temp,2);
end
%% find index in all neurons
idxInAll(12, 1) = 0;
for j =  1 : stimNum
    [~, idxMFR] = max(MFRs(:,j));
    t1 = TSOME.S(1).mean_FRs;
    for i = 1 : unitNum
        t2 = TALL.S(1).mean_FRs;
        if(isequal(t1(idxMFR, 1:50), t2(i, 1:50)))
            idxInAll(j, 1) = i;
        end
    end
end
%% plot most active neurons tuning curve
for i = 1 : stimNum
finalPlotTuningCurve(idxInAll(i), stimNum, TALL, monkeyNum, (i-1)*30)
end

%% functions
function plotTuningCurve(A, a, b, bestToWhich)

x = 0 : 30 : 330;
figure = plot(x, A);
xlabel("degrees of drifting gratings")
ylabel("meanResponse")
title([sprintf("Tuning curve for unit %d of monkey %d",a , b)
    sprintf("Most active neuron in response to %d", bestToWhich)])
saveas(figure,sprintf("%d%d",b, bestToWhich/30+1),'png')
end

function finalPlotTuningCurve(idxInAll, stimNum, TALL, monkeyNum, bestToWhich)
TPlotMat = zeros(1, stimNum);
for i = 1: stimNum
    temp = mean(TALL.S(i).mean_FRs,2);
    TPlotMat(i) = temp(idxInAll); 
end
plotTuningCurve(TPlotMat, idxInAll, monkeyNum, bestToWhich);
end






