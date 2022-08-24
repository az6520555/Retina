clear
close all
path='F:\§Úªº¶³ºÝµwºÐ\Retina exp\exp data\Sorted_final_data\20200419'
rate=20000;
cd(path)
all_file = dir('*.mat') ; % change the type of the files which you want to select, subdir or dir.
files=[85 88 97 94 91]-30;


for z = 1:length(files)
    file = all_file(files(z)).name ;
    [pathstr, name, ext] = fileparts(file);
    directory = [pathstr,'\'];
    filename = [name,ext];
    load([filename]);
    name(name=='_')='-';
    
    a=[];
    a=a_data(1,:);
    a=smooth(a,20);
    
    TS=TimeStamps;
    t=1/rate:1/rate:length(a)/rate;
    datatemp=[];
    datatemp=a(t>TS(1) & t<TS(2));
    datashape=[];
    datashape=datatemp-mean(datatemp);

    aaa=[];
    lags=[];
    [aaa,lags]=autocorr(datashape,'NumLags',100000);
    [M ind]=min(abs(aaa-0.5));
    corrtime(z)=ind/rate;
end