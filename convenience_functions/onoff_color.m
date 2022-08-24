function onoff_color(file)
    %% change axes color to classify on, on-off or off cells; use for multiplot
            rr =[9,17,25,33,41,49,...
          2,10,18,26,34,42,50,58,...
          3,11,19,27,35,43,51,59,...
          4,12,20,28,36,44,52,60,...
          5,13,21,29,37,45,53,61,...
          6,14,22,30,38,46,54,62,...
          7,15,23,31,39,47,55,63,...
            16,24,32,40,48,56];

    load(file)
    for jj=1:length(On_channel)
        subplot(8,8,rr(On_channel(jj)));hold on
        box on
        ax=gca;
        ax.XColor=[0.6 0 0]; % red
        ax.YColor=[0.6 0 0]; % red
        ax.LineWidth = 1;
    end
    for jj=1:length(Off_channel)
        subplot(8,8,rr(Off_channel(jj)));hold on
        box on
        ax=gca;
        ax.XColor=[0 0.6 0]; % green
        ax.YColor=[0 0.6 0]; % green
        ax.LineWidth = 1;
    end
    for jj=1:length(OnOff_channel)
        subplot(8,8,rr(OnOff_channel(jj)));hold on
        box on
        ax=gca;
        ax.XColor=[0.6 0 0.6]; % magenta
        ax.YColor=[0.6 0 0.6]; % magenta
        ax.LineWidth = 1;
    end
end