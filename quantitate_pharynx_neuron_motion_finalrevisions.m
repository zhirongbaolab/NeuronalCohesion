%{
%m1=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_09082019s1e1workspace.mat');
m1=load('D:\kris_neuronscreen\figures_3_31_2020\9082019s1e1_reloaded.mat');
m3=load('D:\kris_neuronscreen\figures_3_31_2020\reloaded_mutant_10232019s2e1.mat');

%m3=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_10232019s2e1workspace.mat');


w1=load('D:\kris_neuronscreen\figrues_1_31_2020\wt_e1workspace.mat');
w2=load('D:\kris_neuronscreen\figrues_1_31_2020\wt_e2workspace.mat');
w3=load('D:\kris_neuronscreen\figrues_1_31_2020\wt_e3workspace.mat');

%m2=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_10232019s1e1workspace.mat');
m2=load('D:\kris_neuronscreen\figures_3_31_2020\1023s1e1withadditionalediting.mat');
%}

%wt list for all embryos
editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABprpappapa'
'ABplpappaaa'
%'ABprpapppap'
     
'ABarapapppa'
'ABarapapppp'    
'ABarapapapp'
'ABalpappppa'    
'ABalpappppp' 
};

anisotropy=1;
%pull z values of 4 first cells on list

dvmutn=[]
dvmutp=[]
dvwtn=[]
dvwtp=[]


dvmutnc=[]
dvmutpc=[]
dvwtnc=[]
dvwtpc=[]

[percellmotion,cells]=getpercellmotionforembryo(w1,editedcellssuls);

%wt
for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        if ~isempty(target)
        if target<=4
            dvwtn=[dvwtn;distance(motion1(1,:)',motion1(end,:)')];
            dvwtnc=[dvwtnc;cumulativeMotion(motion1)];
        else
            dvwtp=[dvwtp;distance(motion1(1,:)',motion1(end,:)')];
            dvwtpc=[dvwtpc;cumulativeMotion(motion1)];
        end
        end
end
end
[percellmotion,cells]=getpercellmotionforembryo(w2,editedcellssuls);

%wt
for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        if ~isempty(target)
        if target<=4
            dvwtn=[dvwtn;distance(motion1(1,:)',motion1(end,:)')];
               dvwtnc=[dvwtnc;cumulativeMotion(motion1)];
        else
            dvwtp=[dvwtp;distance(motion1(1,:)',motion1(end,:)')];
               dvwtpc=[dvwtpc;cumulativeMotion(motion1)];
        end
        end
end
end
[percellmotion,cells]=getpercellmotionforembryo(w3,editedcellssuls);
%wt
for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        if ~isempty(target)
        if target<=4
            dvwtn=[dvwtn;distance(motion1(1,:)',motion1(end,:)')];
               dvwtnc=[dvwtnc;cumulativeMotion(motion1)];
        else
            dvwtp=[dvwtp;distance(motion1(1,:)',motion1(end,:)')];
             dvwtpc=[dvwtpc;cumulativeMotion(motion1)];
        end
        end
end
end






%mutant
%m1=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_09082019s1e1workspace.mat');
%m2=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_10232019s1e1workspace.mat');
%m3=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_10232019s2e1workspace.mat');

anisotropy=1/.16;

%0908s1e1


editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABprpappapa'
'ABplpappaaa'
     
'MSpapapap'
'ABarpaapppa'
'ABalapappa'
'MSaapapp'
'MSaapaapp'
};
[percellmotion,cells]=getpercellmotionforembryo(m1,editedcellssuls);
dvmutn1=[];
dvmutp1=[];
dvmutn1c=[];
dvmutp1c=[];
for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        if target<=4
            dvmutn1=[dvmutn1;distance(motion1(1,:)',motion1(end,:)')];
            dvmutn1c=[dvmutn1c;cumulativeMotion(motion1)];
        else
            dvmutp1=[dvmutp1;distance(motion1(1,:)',motion1(end,:)')];
            dvmutp1c=[dvmutp1c;cumulativeMotion(motion1)];
        end
end
end
%figure; hist(dvmutn1);
%figure; hist(dvmutp1);

%cells for 1023 s1 e1 
%this as to be replaced****
editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABprpappapa'
'ABplpappaaa'
     
'ABalpappppa'
'MSaaapapa'
'MSaaapapp'
'ABalpappppp'
'MSaaaaapp'
};
[percellmotion,cells]=getpercellmotionforembryo(m2,editedcellssuls);
dvmutn2=[];
dvmutp2=[];

dvmutn2c=[];
dvmutp2c=[];

for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        if target<=4
            dvmutn2=[dvmutn2;distance(motion1(1,:)',motion1(end,:)')];
             dvmutn2c=[dvmutn2c;cumulativeMotion(motion1)];
        else
            dvmutp2=[dvmutp2;distance(motion1(1,:)',motion1(end,:)')];
             dvmutp2c=[dvmutp2c;cumulativeMotion(motion1)];
        end
end
end

%figure; hist(dvmutn2);
%figure; hist(dvmutp2);

%cells for 1023 s2 e1
editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABprpappapa'
'ABplpappaaa'
     
'ABaraaapapa'
'MSaapappa'
'MSaapaapp'
'ABaraapaapp'
'ABaraappapp'
};
[percellmotion,cells]=getpercellmotionforembryo(m3,editedcellssuls);
dvmutn3=[];
dvmutp3=[];
dvmutn3c=[];
dvmutp3c=[];

for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        if target<=4
            dvmutn3=[dvmutn3;distance(motion1(1,:)',motion1(end,:)')];
             dvmutn3c=[dvmutn3c;cumulativeMotion(motion1)];
        else
            dvmutp3=[dvmutp3;distance(motion1(1,:)',motion1(end,:)')];
             dvmutp3c=[dvmutp3c;cumulativeMotion(motion1)];
        end
end
end
%figure; hist(dvmutn3);
%figure; hist(dvmutp3);

dvmutn=[dvmutn1;dvmutn2;dvmutn3];
dvmutp=[dvmutp1;dvmutp2;dvmutp3];


dvmutnc=[dvmutn1c;dvmutn2c;dvmutn3c];
dvmutpc=[dvmutp1c;dvmutp2c;dvmutn3c];

%dvmutn=[dvmutn2;dvmutn3];
%dvmutp=[dvmutp2;dvmutn3];


%figure; hist(dvwtn);
%figure; hist(dvwtp);
%figure; hist(dvmutn);
%figure; hist(dvmutp);

figure; errorbar([mean(dvwtn),mean(dvwtp),mean(dvmutn),mean(dvmutp)],[std(dvwtn),std(dvwtp),std(dvmutn),std(dvmutp)]);
title 'displacement'

%figure; errorbar([mean(dvwtnc),mean(dvwtpc),mean(dvmutnc),mean(dvmutpc)],[std(dvwtnc),std(dvwtpc),std(dvmutnc),std(dvmutpc)]);
%title 'cumulative distance'

data=[mean(dvwtp),mean(dvmutp),mean(dvwtn),mean(dvmutn)];
stdev=[std(dvwtp),std(dvmutp),std(dvwtn),std(dvmutn)];
x=[1,2,3,4]
figure;
bar(x,data);
hold on
er=errorbar(x,data,stdev);
er.LineStyle='none';
title 'displacement'

%[h,p,c]=ttest(dvwtp,dvmutp)
%unpaired
[h,p,c]=ttest2(dvwtp(1:end-1),dvmutp)
[h,p,c]=ttest2(dvwtn,dvmutn)

%paired
[h,p,c]=ttest(dvwtp(1:end-1),dvmutp)
[h,p,c]=ttest(dvwtn,dvmutn)


%[h,p,c]=vartest2([dvwtp(1:end-1)],dvmutp)
%[h,p,c]=vartest2(dvwtn,dvmutn)



%try just smdds
newneuronwt=dvwtn([1,2,5,6,9,10]);
newneruonmt=dvmutn([1,2,5,6,9,10]);
[h,p,c]=ttest2(newneuronwt,newneruonmt)

%now try pharynx subset too

data=[mean(dvwtp),mean(dvmutp),mean(newneuronwt),mean(newneruonmt)];
stdev=[std(dvwtp),std(dvmutp),std(newneuronwt),std(newneruonmt)];
x=[1,2,3,4]
figure;
bar(x,data);
hold on
er=errorbar(x,data,stdev);
er.LineStyle='none';
title 'displacement'

return
