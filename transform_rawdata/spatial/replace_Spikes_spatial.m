% replace Spikes by sorted data
clear 
cd('E:\Chou\chou20210323\merge\'); % open original data
path_sort=['E:\Chou\chou20210323\sort\']; % sorted spikes data
savepath=['E:\Chou\chou20210323\sort_complete_data\']; % new data
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
unit_number=3;  

for z = 1:n_file
    file = all_file(z).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    for u=1:unit_number
        clearvars Spikes
        try
            load([path_sort,filename(7:end-4),'_unit',num2str(u),ext]) % neglect "merge_"(6) and ".mat"(4)
        catch
            continue
        end
        % cut pre-stimulus duration usin TimeStamps
        if length(TimeStamps)==2
            for ch=1:60
                Spikes{ch}=Spikes{ch}(Spikes{ch}>TimeStamps(1) & Spikes{ch}<TimeStamps(2)) - TimeStamps(1);
            end
        end
        
        sorted_spikes=Spikes;
        if exist('reconstruct_spikes','var') & exist('type','var')
            save([savepath,filename(1:end-4),'_sort_unit',num2str(u),'.mat'],'sorted_spikes','bin_pos','TimeStamps','reconstruct_spikes','diode_BT','BinningInterval','type')
        else
            save([savepath,filename(1:end-4),'_sort_unit',num2str(u),'.mat'],'sorted_spikes','bin_pos','TimeStamps','diode_BT','BinningInterval')
        end
    end
end
    
    