%% load data
clc
clear all
close all
T = load('spikes_gratings/data_monkey1_gratings');
TALL = load('spikes_gratings/S_monkey1');
TSOME = load('spikes_gratings/Snew_monkey1');
Unum = [106 88 112];
%% set parameters
monkeyNum = 1;
newNum = size(TSOME.S(1).mean_FRs,1);
stimNum = 12;
unitNum = Unum(monkeyNum);
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
%%
[MFRsMax, MFRsMaxIndex] = max(MFRs,[],2);
%% MAP
mainMAP = T.data.MAP;
mainChannels = T.data.CHANNELS;
MAP = zeros(10);
MAP1 = zeros(10);
for i =  1 : 10
    for j = 1 : 10
        for k = 1 : unitNum
            if(mainMAP(i, j) == mainChannels(k, 1))
                for l = 1 : newNum
                    if(idxInAll(l) == k)
                        if(MAP(i, j) == 0)
                            MAP(i, j) = MFRsMaxIndex(l);
                            MAP1(i, j) = k;
                            fprintf('i = %d, j = %d, l = %d, MFRsMaxIndex(l) = %d\n',i, j, l, MFRsMaxIndex(l));
                        elseif(MAP(i, j) ~= 0)
                            for ii = 1 : newNum
                                if(idxInAll(ii) == MAP1(i, j))
                            if(MFRs1(k, MFRsMaxIndex(l)) > MFRs1(MAP1(i, j), MFRsMaxIndex(ii)))
                                MAP(i, j) = MFRsMaxIndex(l);
                                MAP1(i, j) = k;                            
                                fprintf('i = %d, j = %d, l = %d, MFRsMaxIndex(l) = %d\n',i, j, l, MFRsMaxIndex(l));
                            end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
%%
figure
im = image(30*MAP(:, :));
im.CDataMapping = 'scaled';
colorbar
title("Prefered orientations for Monkey 1", 'interpreter', 'latex')
xlabel("X", 'interpreter', 'latex')
ylabel("Y", 'interpreter', 'latex')



%% functions



