function[]=mcd_reader(varargin)
if nargin > 0 
	eval(varargin{1}); %execute the argument 
else
    begin;
end;
end

function begin()
%     cd D:\Michael
    close all
    global h
    h = figure(1);
    % set(h,'units','normalize','position',[0 0 1 1]);

    % uicontrol('Style','text', 'Units','normalized', ...
    %           'Position',[0.02  0.18 0.07 0.05],...
    %           'String', 'time');
    % uicontrol('Style','edit', 'Units','normalized', ...
    %           'Position',[0.02  0.18 0.07 0.05],...
    %           'String', 'time');

    uicontrol('Style','slider','units','normalize','position',[0.1,0.05,0.8,0.05],'tag','slide',...
        'Callback','mcd_reader(''move_slide'')');
    
    uicontrol('Style','text','units','normalize','string','0','position',[0.078,0.01,0.06,0.04],'tag','starttime','fontunits','normalized')
    parentcolor = get(get(findobj('tag','starttime'), 'parent'),'color');
    set(findobj('tag','starttime'),'BackgroundColor',parentcolor)
    uicontrol('Style','text','units','normalize','string','1','position',[0.88,0.01,0.06,0.04],'tag','endtime','BackgroundColor',parentcolor,'fontunits','normalized')
%     set(findobj('tag','endtime'),'BackgroundColor',get(get(findobj('tag','starttime'), 'parent'), 'color'))
    uicontrol('Style','text','units','normalize','string','','position',[0.135,0.01,0.06,0.04],'tag','nowtime','BackgroundColor',parentcolor,'fontunits','normalized')
%     set(findobj('tag','nowtime'),'BackgroundColor',get(get(findobj('tag','starttime'), 'parent'), 'color'))
    
    uicontrol('Style','text','units','normalize','string','trial','position',[0.0,0.94,0.08,0.04],'BackgroundColor',parentcolor,'fontunits','normalized')
    uicontrol('Style','edit','units','normalize','string','1',    'position',[0.09 ,0.94,0.07,0.04],'tag','trial','fontunits','normalized',...
              'Callback','mcd_reader(''next_trial(''''change'''')''  )')
    uicontrol('Style','text','units','normalize','string','plot_time','position',[0.0,0.90,0.08,0.04],'BackgroundColor',parentcolor,'fontunits','normalized')
    uicontrol('Style','edit','units','normalize','string','1000'     ,'position',[0.09 ,0.90,0.07,0.04],'tag','show_time','fontunits','normalized',...
              'Callback','mcd_reader(''move_slide'')')
    uicontrol('Style','text','units','normalize','string','ms'       ,'position',[0.165,0.90,0.025,0.04],'BackgroundColor',parentcolor,'fontunits','normalized')
    uicontrol('Style','text','units','normalize','string','yrange','position',[0.0,0.86,0.08,0.04],'BackgroundColor',parentcolor,'fontunits','normalized')
    uicontrol('Style','popupmenu','units','normalized','tag','yrange','position',[0.09,0.86,0.1,0.04],'fontunits','normalized','String',{'500','1000','1500','2000','3000'},'fontunits','normalized',...
              'Callback','mcd_reader(''move_slide'')','value',2)
    
    
    % hListener = addlistener(hSlider,'Value','PostSet',@(x,y) f('123'));
    uicontrol('Style','pushbutton', 'Units','normalized', 'position',[0.03,0.12,0.1,0.07],'string','select file','Callback','mcd_reader(''select_file'')','fontunits','normalized');
    
    uicontrol('Style','pushbutton', 'Units','normalized', 'position',[0.03,0.05,0.05,0.05],'string','<<',...
              'Callback', 'mcd_reader(''next_trial(''''last'''')''  )');
    uicontrol('Style','pushbutton', 'Units','normalized', 'position',[0.92,0.05,0.05,0.05],'string','>>',...
              'Callback', 'mcd_reader(''next_trial(''''next'''')'')');
    
    %xtick label
    uicontrol('Style','text','units','normalized','position',[0.175, 0.94,0.05,0.02],'string','0','tag','xtick1','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.225, 0.94,0.05,0.02],'string','500','tag','xtick3','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.275, 0.94,0.05,0.02],'string','0','tag','xtick1','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.325, 0.94,0.05,0.02],'string','500','tag','xtick3','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.375, 0.94,0.05,0.02],'string','0','tag','xtick1','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.425, 0.94,0.05,0.02],'string','500','tag','xtick3','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.475, 0.94,0.05,0.02],'string','0','tag','xtick1','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.525, 0.94,0.05,0.02],'string','500','tag','xtick3','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.575, 0.94,0.05,0.02],'string','0','tag','xtick1','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.625, 0.94,0.05,0.02],'string','500','tag','xtick3','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.675, 0.94,0.05,0.02],'string','0','tag','xtick1','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
    uicontrol('Style','text','units','normalized','position',[0.725, 0.94,0.05,0.02],'string','500','tag','xtick3','BackgroundColor',parentcolor,'fontunits','normalized','fontsize',1)
          
    uicontrol('style','checkbox','units','normalized','position',[0.81,0.9,0.04,0.05],'tag','filter','BackgroundColor',parentcolor,'fontunits','normalized','Callback','mcd_reader(''move_slide'')')
    uicontrol('style','text',    'units','normalized','position',[0.84,0.9,0.04,0.04],'string','filter','BackgroundColor',parentcolor,'fontunits','normalized')
    uicontrol('style','checkbox','units','normalized','position',[0.81,0.85,0.04,0.05],'tag','analog','BackgroundColor',parentcolor,'fontunits','normalized','Callback','mcd_reader(''move_slide'')')
    uicontrol('style','text',    'units','normalized','position',[0.84,0.85,0.06,0.04],'string','analog','BackgroundColor',parentcolor,'fontunits','normalized')
    uicontrol('style','checkbox','units','normalized','position',[0.81,0.95,0.04,0.05],'tag','period_on','BackgroundColor',parentcolor,'fontunits','normalized','Callback','mcd_reader(''next_trial(''''change'''')''  )')
    uicontrol('style','text',    'units','normalized','position',[0.84,0.95,0.06,0.04],'string','period','BackgroundColor',parentcolor,'fontunits','normalized')
    uicontrol('Style','edit',    'units','normalize', 'string','1',    'position',[0.89,0.95,0.06,0.04],'tag','period','fontunits','normalized',...
              'Callback','mcd_reader(''next_trial(''''change'''')''  )')
    for i = 1:64
        if any(i==[1,8,57])
            continue
        end
        j = (mod(i-1,8))*8+floor((i-1)/8)+1;
        h(i) = axes('units','normalized','position',[0.1+1/10*floor((i-1)/8),0.95-1/10*(mod(i-1,8)+1+0.1),1/10.3,1/10.3]);
    %     title(num2str(i))
        set(gca,'xticklabel',[],'yticklabel',[])
        set(gca,'xtick',[0,1,2,3],'xlim',[0 4])
    end
    h(1:6) = h(2:7); h(7:54) = h(9:56); h(55:61) = h(58:64);h(62:64) = [];
end

function select_file()
    clear global starts a_data
    global file SamplingRate a_data starts h t1 t2
    clear starts
%     file = ('D:\Michael\c2d_rest0=5_rest1=4_dur=0.5_bar.mcd');
    [filename,pathname] = uigetfile('*.mcd');
    cd(pathname)
    if ~ischar(filename)
        return
    end
    file = strcat(pathname,filename);
    AllDataInfo =datastrm(file);
        SamplingRate = getfield(AllDataInfo,'MillisamplesPerSecond2');
        SamplingRate = SamplingRate(1)/1000;
        t1=getfield(AllDataInfo,'sweepStartTime');
        t2=getfield(AllDataInfo,'sweepStopTime') ;
    StartStopVector = [t1 , t2];
    a_data = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Analog Raw Data');
    a_data = a_data.data;
    if size(a_data,1) >1
        a_data = a_data(2,:);
    end
    if isempty(a_data)
        a_data = zeros(1,t2/1000*SamplingRate);
    end
    [pks,TimeStamps]=findpeaks(diff(a_data),'MINPEAKHEIGHT',5*std(diff(a_data))); 
    if isempty(pks); TimeStamps=[]; 
    elseif mean(pks<50); TimeStamps=[]; end;
    if ~isempty(TimeStamps)
        if std(diff(TimeStamps)) <30000     % too uniform, may be only trigger on starts
            starts = TimeStamps;
        elseif ~isempty(diff(TimeStamps))
            starts = find_starts(TimeStamps);
        end
    end
    starts = floor(starts/SamplingRate*1000);
    set(findobj('tag','slide'),'value',0);
    move_slide;
%     StartStopVector = [starts(1) , starts(1)+SamplingRate]
%     AllData = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Electrode Raw Data');
%     datas = AllData.data;
%     for i = 1:60
%         plot(h(i),[0:1/SamplingRate:1], datas(i,starts(1):starts(1)+SamplingRate));
%         set(h(i),'xtick',[],'ytick',[])
%     end
%     
%     
    set(findobj('tag','endtime'),'string',(starts(2)-starts(1))/1000);
end

function move_slide()
    global file SamplingRate starts t1 t2 h a_data 
    
    t_start = get(findobj('tag','slide'),'value');
    t_show = str2num(get(findobj('tag','show_time'),'string'));  % in ms   
%     t_show = t_show ;    % in points

    AllDataInfo =datastrm(file);
    % StartStopVector = [get(findobj('tag','t_start'),'value'),get(findobj('tag','t_end'),'value')];
    trial = str2num(get(findobj('tag','trial'),'string'));
    if isempty(starts)
        starts=[t1,t2];
        if trial ~=1
            trial=1;
            disp('it seems there is no analog data, so trial is set to 1');
        end
    end
    if ~get(findobj('tag','period_on'),'value')
        period = starts(2) - starts(1);
    else
        period = str2num(get(findobj('tag','period'),'string'));
        period = period*1000;
        starts = 0:period:t2;
    end
    
    % if get(findobj('tag','breaklimit'),'value') == 0;
    %     StartStopVector = [start
    if isempty(trial) || ~exist('trial','var')
        trial = 2;
    end
    if trial>length(starts)
        fprintf('only %d trials\n',length(starts));
        return
    end
    if (1-t_start)*period <t_show/1000
        t_start = 1 - t_show/period/1000;
    end

    if isempty(t_show) || ~exist('t_show','var')
        t_show = 1000;    % all units in ms for startend variable
    end
    StartStopVector = [starts(trial)+period*t_start , starts(trial)+period*t_start+t_show]
    AllData = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Electrode Raw Data');
    datas = AllData.data;
    
    t_monitor = [period*t_start, period*t_start+t_show]/1000*SamplingRate;
    yrange = get(findobj('tag','yrange'),'string');
    yrange = str2num(yrange{get(findobj('tag','yrange'),'value')});
    
    
    a_data_plot = nextdata(AllDataInfo,'startend',StartStopVector,'originorder','on' ,'streamname','Analog Raw Data');
    a_data_plot = a_data_plot.data;
    if ~isempty(a_data_plot)
        plot(h(61),[t_monitor(1)/SamplingRate:1/SamplingRate:(t_monitor(2)-1)/SamplingRate], a_data_plot(1:t_show/1000*SamplingRate));
        set(h(61),'xticklabel',[],'yticklabel',[],'xlim',[t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate],'ylim',[min(a_data_plot)-1000,max(a_data_plot)+1000])
    end
    
    
    if get(findobj('tag','analog'),'value')==0
        if get(findobj('tag','filter'),'value')==0
            for i = 1:60
                plot(h(i),[t_monitor(1)/SamplingRate:1/SamplingRate:(t_monitor(2)-1)/SamplingRate], datas(i,1:t_show/1000*SamplingRate));
                v_zero = mean(datas(i,:));
                set(h(i),'xtick',linspace(t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate,5))
                set(h(i),'xticklabel',[],'yticklabel',[],'xlim',[t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate],'ylim',[v_zero-yrange/2,v_zero+yrange/2])
            end
        else
            [b,a] = butter(2,200/10000,'high'); % set butter filter
            for i = 1:60
                FilterData = filter(b,a,datas(i,:));
                plot(h(i),[t_monitor(1)/SamplingRate:1/SamplingRate:(t_monitor(2)-1)/SamplingRate], FilterData(1:t_show/1000*SamplingRate));
                v_zero = mean(FilterData(:));
                set(h(i),'xtick',linspace(t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate,5))
                set(h(i),'xticklabel',[],'yticklabel',[],'xlim',[t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate],'ylim',[v_zero-yrange/2,v_zero+yrange/2])
            end
        end
    else
        if get(findobj('tag','filter'),'value')==0
            for i = 1:60
                [h0,~,h2]=plotyy(h(i),[t_monitor(1)/SamplingRate:1/SamplingRate:(t_monitor(2)-1)/SamplingRate], datas(i,1:t_show/1000*SamplingRate),[t_monitor(1)/SamplingRate:1/SamplingRate:(t_monitor(2)-1)/SamplingRate], a_data_plot(1:t_show/1000*SamplingRate));
                v_zero = mean(datas(i,:));
                set(h2,'color',[0.3 0.3 0.3])
                set(h0,'xtick',linspace(t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate,5))
                set(h0,'xticklabel',[],'yticklabel',[])
                set(h0(1),'xlim',[t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate],'ylim',[v_zero-yrange/2,v_zero+yrange/2])
                set(h0(2),'xticklabel',[],'yticklabel',[],'xlim',[t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate],'ylim',[min(a_data_plot)-3000,max(a_data_plot)+3000])
            end
        else
            [b,a] = butter(2,200/10000,'high'); % set butter filter
            for i = 1:60
                FilterData = filter(b,a,datas(i,:));
                [h0,~,h2]=plotyy(h(i),[t_monitor(1)/SamplingRate:1/SamplingRate:(t_monitor(2)-1)/SamplingRate], FilterData(1:t_show/1000*SamplingRate),[t_monitor(1)/SamplingRate:1/SamplingRate:(t_monitor(2)-1)/SamplingRate], a_data_plot(1:t_show/1000*SamplingRate));
                v_zero = mean(FilterData(:));
                set(h2,'color',[0.3 0.3 0.3])
                set(h0,'xtick',linspace(t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate,5))
                set(h0,'xticklabel',[],'yticklabel',[])
                set(h0(1),'xlim',[t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate],'ylim',[v_zero-yrange/2,v_zero+yrange/2])
                set(h0(2),'xticklabel',[],'yticklabel',[],'xlim',[t_monitor(1)/SamplingRate,(t_monitor(2)-1)/SamplingRate],'ylim',[min(a_data_plot)-3000,max(a_data_plot)+3000])
            end
        end
    end



    
    
    %     t_start
    set(findobj('tag','nowtime'),'string',sprintf('%.2f',((period*t_start)/1000)),'position',[0.135+0.68*t_start,0.01,0.04,0.04])
    objs = findobj('tag','xtick1');
    for i = 1:length(objs)
        set(objs(i),'string',sprintf('%.2f',t_monitor(1)/SamplingRate))
    end
    objs = findobj('tag','xtick3');
    for i = 1:length(objs)
        set(objs(i),'string',sprintf('%.2f',(t_monitor(2)+t_monitor(1))/2/SamplingRate))
    end
end

function next_trial(s)
global starts t2  a_data SamplingRate
    trial = floor(str2num(get(findobj('tag','trial'),'string')));
    if strcmp(s,'next') & trial<length(starts)
        trial = trial +1;
    elseif strcmp(s,'last') & trial > 1
        trial = trial-1;
    elseif strcmp(s,'change')
        if trial < 1
            trial = 1;
        elseif trial > length(starts)
            trial = length(starts);
        end
    end
    if trial<0 |trial > length(starts)
        error('something wrong in next_trial(s)')
    end
    set(findobj('tag','trial'),'string',num2str(trial));
%     manual_period = get(findobj('tag','period_on'),'value');
    if get(findobj('tag','period_on'),'value')
        period = str2num(get(findobj('tag','period'),'string'));
        period = period*1000;
        if trial>t2/period
            trial = floor(t2/period);
            set(findobj('tag','trial'),'string',num2str(trial));
            set(findobj('tag','endtime'),'string',period/1000);
        end
        set(findobj('tag','endtime'),'string',num2str(period/1000));
    else
        [pks,TimeStamps]=findpeaks(diff(a_data),'MINPEAKHEIGHT',5*std(diff(a_data))); 
        if isempty(pks); TimeStamps=[]; 
        elseif mean(pks<50); TimeStamps=[]; end;
        if ~isempty(TimeStamps)
            if std(diff(TimeStamps)) <1000     % too uniform, may be only trigger on starts
                starts = TimeStamps;
            elseif ~isempty(diff(TimeStamps))
                starts = find_starts(TimeStamps);
            end
        end
        starts = floor(starts/SamplingRate*1000);
        set(findobj('tag','endtime'),'string',num2str((starts(2)-starts(1))/1000));
        
    end
    move_slide()
end
% 
% uicontrol('Style','pushbutton', 'Units','normalized', ...
%           'Position',[0.05  0.8 0.25 0.15],...
%           'String', 'set init by MAP','Callback','gui_bar(''get_init_map'');');
