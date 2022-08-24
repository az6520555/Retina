path=['F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20190131'];
cd(path);
all_file = dir('*.mat');
for i=7
    load(all_file(i).name)
    for j=1:60
        if length(Spikes{j})<100
            Spikes{j}=[];
        end
    end
end