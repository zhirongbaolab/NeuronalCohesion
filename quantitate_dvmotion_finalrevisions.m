
%pull z values of 4 first cells on list

dvmut=[]
dvmutalt=[]
for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        if target<=4
            dvmut(target)=motion1(1,3)-motion1(end,3);
            dvmutalt(target)=sum(abs(motion1(1:end-1,3)-motion1(2:end,3)));
        end
end
end


dvwt=[]
dvwtalt=[]
for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        if target<=4
            dvwt(target)=motion1(1,3)-motion1(end,3);
            [motion1(1,:),motion1(end,:)]
             dvwtalt(target)=sum(abs(motion1(1:end-1,3)-motion1(2:end,3)));

        end
end
end

%these two embs are opposite so invert one
dvwttest=dvwt*-1;
mean(dvwttest)
mean(dvmut)

[h,p,c]=ttest(allmutant,allwt)
[h,p,c]=vartest2(allmutant,allwt)
p=anova1([allmutant',allwt'])