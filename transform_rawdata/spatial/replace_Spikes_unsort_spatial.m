% replace Spikes of unsort checkerboard
clear 
path_merge='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210331\merge\'; % open original data
path_spike=['F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Spatial stimuli\20210331\data\']; % spikes data
cd(path_merge)
mkdir AddSpikes
savepath=[path_merge,'AddSpikes\']; % new data
cd(savepath)

filenames={'0224_Checkerboard_30Hz_27_30min_Br50_Q100','0224_Checkerboard_30Hz_27_15min_Br50_Q100'};

for i=1:length(filenames)
    load([path_spike,filenames{i}])
    clearvars a_data Infos TimeStamps
    load([path_merge,'merge_',filenames{i}])
    
    if length(TimeStamps)==2
        for ch=1:60
            Spikes{ch}=Spikes{ch}(Spikes{ch}>TimeStamps(1) & Spikes{ch}<TimeStamps(2)) - TimeStamps(1);
        end
    end
    reconstruct_spikes=Spikes;
    save([filenames{i},'_AddSpikes.mat'],'bin_pos','BinningInterval','data_name','diode_BT','pass','pwd','TimeStamps','reconstruct_spikes')
end