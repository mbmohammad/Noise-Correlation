%% load data
clc
clear all
close all
T = load('spikes_gratings/data_monkey1_gratings');
TALL = load('spikes_gratings/S_monkey1');
Unum = [106 88 112];
%% set parameters
monkeyNum = 1;
stimNum = 12;
unitNum = Unum(monkeyNum);
%%
MFRs1 = zeros(unitNum, stimNum);
for i = 1 : stimNum
    temp = TALL.S(i).mean_FRs;
    MFRs1(:, i) = mean(temp(:, 1:35),2);
end
%%
[MFRsMax, MFRsMaxIndex] = max(MFRs1,[],2);
%% MAP
mainMAP = T.data.MAP;
mainChannels = T.data.CHANNELS;
MAP = zeros(10);
MAP1 = zeros(10);
%%
for i =  1 : 10
    for j = 1 : 10
        for k = 1 : unitNum
            if(mainMAP(i, j) == mainChannels(k, 1))
                if(MAP(i,j) == 0)
                   MAP(i,j) =  MFRsMaxIndex(k);
                   MAP1(i,j) = k;
                elseif(MAP(i, j)~= 0)
                    temp = mean(MFRs1, 2);
                    if(temp(k) > temp(MAP1(i,j)))
                        MAP(i,j) =  MFRsMaxIndex(k);
                        MAP1(i,j) = k;
                    end
                end
            end
        end
    end 
end
%%
im = image(MAP(:, :));
im.CDataMapping = 'scaled';
colorbar
%%
for i =  1: 12
t = mean(TALL.S(i).mean_FRs);
figure
plot(t)
end
















