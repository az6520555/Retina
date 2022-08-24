%% if the data is sorted
filename='0319_HMM_UD_G20_7min_Br50_Q100_8_sort.mat';
path_sort=['E:\Chou\ssss\'];
% ================== new sorting sequence ========================
    ss=[29,30,28,27,22,21,14,20,...
        13,6,12,5,19,11,4,3,10,...
        18,2,9,1,8,17,7,16,15,26,...
        25,23,24,32,31,33,34,39,...
        40,47,41,48,55,49,56,42,...
        50,57,58,51,43,59,52,60,...
        53,44,54,45,46,35,36,38,37];

% ================= for sorted data in .mat files ==========================
    Spikes = cell(1,60);
    load([path_sort,filename]);
    temp_spikes={};
    for h=1:60
        if h<11
            temp_spikes{ss(h)} = eval(['adc00',int2str(h-1)]);
        else
            temp_spikes{ss(h)} = eval(['adc0',int2str(h-1)]);
        end
    end
    for i=1:60
        if isempty(temp_spikes{i})==1
            continue
        end
        for j=1:length(temp_spikes{i}(:,1))
            if temp_spikes{i}(j,3)==0 % this determine which unit we choose
                Spikes{i}=[Spikes{i} temp_spikes{i}(j,1)];
            end
        end
    end


%% OnOff response raster plot
% clear
cd('E:\Chou\ssss')
filename='0319_HMM_UD_G20_7min_Br50_Q100_8.mat'
load(filename)
for i=1:60
    if isempty(Spikes{i})==1
        Spikes{i}=0;
    end
end
figure(87);
LineFormat.Color = [1 0.3 0.3];
plotSpikeRaster(Spikes,'PlotType','vertline','LineFormat',LineFormat)
hold on