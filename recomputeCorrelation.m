function datarev=recomputeCorrelation(datarev)
%recompute correlation
cells=datarev.cells;
%want to parse these not ones in original file
%editedcellssuls=loaded_dataset.editedcellssuls;
tstart=datarev.tstart;
anisotropy=datarev.anisotropy;


%chunk list of cells

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


%chunk list of cells
editedcellssuls={
'ABprpapaaaa'
'ABprpapaapp'
'ABprpapaapa'
'ABprpapaaap'
'ABprpaapppp'
'ABplpapaaaa'
'ABplpapaapp'
'ABplpapaapa'
'ABplpapaaap'
'ABplpaapppp'

};

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

corrv=[];
for i=1:length(percellmotion)
    if ~isempty(percellmotion{1,i})
        for j=1:length(percellmotion)
            if ~isempty(percellmotion{1,j})&&j~=i
                motion1=percellmotion{1,i}; 
                motion2=percellmotion{1,j};
                motion1(:,3)=motion1(:,3).*anisotropy;
                motion2(:,3)=motion2(:,3).*anisotropy;
                motion1=smoothdata(motion1);
                motion2=smoothdata(motion2); 
                v1=motion1(1:end-1,:)-motion1(2:end,:);
                v2=motion2(1:end-1,:)-motion2(2:end,:);
                 ind1=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        ind2=(find(strcmp(editedcellssuls,cells{j}.name), 1));
        minlength=min(length(v1),length(v2));
        %if one is shorter bc of tracking hiccup chop off end
        coormat=corrcoef(v1(1:minlength,:),v2(1:minlength,:));
                corrv(ind1,ind2)=coormat(1,2);
            else
            
            end
        end
        
    end
end
for i=1:length(corrv)
    corrv(i,i)=1;
end

%store only thing we need from here
datarev.corrv=corrv;

end

