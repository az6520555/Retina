% replace Spikes by sorted data
clear 
cd('\\192.168.0.102\Public\Retina\Chou\Exp\20220901\SplitData\'); % open original data
path_sort=['\\192.168.0.102\Public\Retina\Chou\Exp\20220901\sort\']; % sorted spikes data
savepath=['\\192.168.0.102\Public\Retina\Chou\Exp\20220901\sorted_data_final\']; % new data
mkdir(savepath)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
n_file = length(all_file); 
unit_number=4;  

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
            load([path_sort,filename(1:end-4),'_unit',num2str(u),ext])
        catch
            continue
        end
        save([savepath,filename(1:end-4),'_sort_unit',num2str(u),'.mat'],'a_data','Spikes','TimeStamps')
    end
end
    
    