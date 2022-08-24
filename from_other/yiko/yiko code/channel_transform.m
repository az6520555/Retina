clear all
%%offline number to Tina number(mea number)
rNumber=[59,42,46,7,10,52,...
    55,38,21,25,28,31,14,56,...
    13,34,17,3,49,35,18,39,...
    9,48,30,51,60,22,4,43,...
    6,27,45,24,15,53,11,32,...
    2,41,58,12,26,40,57,36,...
    20,37,54,50,47,44,1,19,...
        16,33,29,8,5,23];

offline_ch=[6 ];  %offline channel number
series=[];
series(1,1:length(offline_ch))=offline_ch;
count=1;
for ch=offline_ch
series(2,count)=find(rNumber==ch); %Tina's number
count=count+1;
end
series(3,:)=sort(series(2,:));
series(3,:)


%Tina number to offline number
Tina_ch=[ 31 9 56 51 49 42 28 ];
count=1;
this=[];
for ch=Tina_ch
this(1,count)=rNumber(ch); %Tina's number
count=count+1;
end
this

sort(this)



