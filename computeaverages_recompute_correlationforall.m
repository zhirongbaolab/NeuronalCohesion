%coallate and average correlation maps

m1=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_09082019s1e1workspace.mat');
m2=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_10232019s1e1workspace.mat');
m3=load('D:\kris_neuronscreen\figrues_1_31_2020\mutant_10232019s2e1workspace.mat');
w1=load('D:\kris_neuronscreen\figrues_1_31_2020\wt_e1workspace.mat');
w2=load('D:\kris_neuronscreen\figrues_1_31_2020\wt_e2workspace.mat');
w3=load('D:\kris_neuronscreen\figrues_1_31_2020\wt_e3workspace.mat');

%test recomputation/anisotropy effect
%reorder
%recolor

editedcells={
'SMDDR'
'SIBVR'
'SIADR'
'AIYR'
'ABprpaapppp'
'SMDDL'
'SIBVL'
'SIADL'
'AIYL'
'ABplpaapppp'
};


m1=recomputeCorrelation(m1);
m2=recomputeCorrelation(m2);
m3=recomputeCorrelation(m3);
w1=recomputeCorrelation(w1);
w2=recomputeCorrelation(w2);
w3=recomputeCorrelation(w3);
%}

figure;
subplot(1,3,1);
imagesc(w1.corrv);
subplot(1,3,2);
imagesc(w2.corrv);
subplot(1,3,3);
imagesc(w3.corrv);

figure;
subplot(1,3,1);
imagesc(m1.corrv);
subplot(1,3,2);
imagesc(m2.corrv);
subplot(1,3,3);
imagesc(m3.corrv);

wildaverage=(w1.corrv+w2.corrv+w3.corrv)./3;

mutantaverage=(m1.corrv+m2.corrv+m3.corrv)./3;

%blotting above diagonal
%{
for i=1:length(wildaverage)
    for j=1:length(wildaverage)
        if i<=j
            wildaverage(i,j)=1;
            mutantaverage(i,j)=1;
        end
    end
end
%}

figure;imagesc(wildaverage);
yticks(linspace(1,length(editedcells),length(editedcells)));

yticklabels(w1.editedcells)

figure;imagesc(mutantaverage);
yticks(linspace(1,length(editedcells),length(editedcells)));

yticklabels(editedcells)
%statstest
allwt=[];
allmutant=[];

for i=1:length(wildaverage)
    for j=1:length(wildaverage)
        if i<j
          allwt=[allwt;w1.corrv(i,j);w2.corrv(i,j);w3.corrv(i,j)];
          allmutant=[allmutant;m1.corrv(i,j);m2.corrv(i,j);m3.corrv(i,j)];
        end
    end
end

%figure; hist (allwt)
%figure;hist(allmutant);

[h,p,c]=ttest(abs(allmutant),abs(allwt))
%unpaired is worse
[h,p,c]=ttest2(allmutant,allwt)
%[h,p,c]=vartest2(allmutant,allwt)
%p=anova1([allmutant',allwt'])


