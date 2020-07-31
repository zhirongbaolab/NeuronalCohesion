function [percellmotion,cells]=getpercellmotionforembryo(loaded_dataset,editedcellssuls)

cells=loaded_dataset.cells;
%want to parse these not ones in original file
%editedcellssuls=loaded_dataset.editedcellssuls;
tstart=loaded_dataset.tstart;
anisotropy=loaded_dataset.anisotropy;


percellmotion={};
%for all cells

%parse data into target cells 
%missing=1321
for i=1:length (cells)
%if in time range
istarget=~isempty(find(strcmp(editedcellssuls,cells{i}.name), 1));
isaliveduringtime=cells{i}.endtime>tstart;
if i==1321
'test'
end
    if(istarget& isaliveduringtime)
        if i==1321
'test'
end
        %note cells can end early by tracking error/death or
        %can start late because of nonterminal period
        %technically could be nonterminal (some are) and then divide 
        %hence be mising some at end but dont expect this 
        if (cells{i}.endtime-length(cells{i}.exp)>tstart)
            %cell is born during time window copy relevant portion of
            %parent and concatenate
            parentbit=[];
            for j=1:length(cells)
                %is parent
                if(strcmp(cells{j}.name,cells{i}.name(1:end-1)))
                    %compute start point in parent
                    startindex=cells{j}.endtime-length(cells{j}.exp);%start time of j
                    startindex=tstart-startindex+1;
                    parentbit=cells{j}.pos(startindex:end,:);
                end
            end
            ['found ',cells{i}.name]
            percellmotion{i}=[parentbit;cells{i}.pos];
      
        else
        ['found ',cells{i}.name]
        %if born at beginning of period there is some to trim off
        startindex=cells{i}.endtime-length(cells{i}.exp);%start time of j
        startindex=tstart-startindex+1;
        
        percellmotion{i}=cells{i}.pos(startindex:end,:);
        end
    else
        percellmotion{i}=[];
    end
end
cellcolors={[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],'g','g','g','g','g','g'};



figure
hold on

for i=1:400
    scatter3(cells{1,i}.pos(1,1),cells{1,i}.pos(1,2),cells{1,i}.pos(1,3)*anisotropy,'.k');
end
for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
        motion1=smoothdata(motion1);
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        
 %       plot3( motion1(:,1), motion1(:,2), motion1(:,3),'Color',cellcolors{target},'LineWidth',4);
     maxwidth=16; %version for fig 2
     maxwidth=14; %version for fig 3
 for j=1:length(motion1)-1
  
     widthc=maxwidth*j/(length(motion1)-1);
     if j>1
         plot3( motion1(j-1:j+1,1), motion1(j-1:j+1,2), motion1(j-1:j+1,3),'Color',cellcolors{target},'LineWidth',widthc)
    
     else
    plot3( motion1(j:j+1,1), motion1(j:j+1,2), motion1(j:j+1,3),'Color',cellcolors{target},'LineWidth',widthc)
     end
 end
     
               scatter3(motion1(end,1),motion1(end,2),motion1(end,3),3000,cellcolors{target},'.');
    %text(motion1(end,1),motion1(end,2),motion1(end,3),editedcells{target},'FontSize',17);

end
end
axis equal
axis vis3d
title(loaded_dataset.emb);




end

