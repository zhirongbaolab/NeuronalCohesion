%smddl/r rmev
%wt pairing
pairings={'ABplpapaaaa' 'ABalpappppa'    'ABalpappppp'    'ABarapapppa'
  'ABprpapaaaa' 'ABarapapppp'    'ABarapapppa'    'ABarapapapp'
   'ABplpappaaa' 'ABalpappppa'    'ABalpappppp'    'ABarapapapp'
   }

 pairings={'ABplpapaaaa' 'ABalpappppa'
        'ABprpapaaaa' 'ABarapapppa'}
    
%manual pairing for mutant 1
%{
pairings={'ABplpapaaaa' 'ABaraappapp'    'ABalpappapa'    'ABaraapaaap'
        'ABprpapaaaa' 'MSpapaaap'    'ABarapapppp'   'MSpapaaaa'}
    
    %mt single
 pairings={'ABplpapaaaa' ''
        'ABprpapaaaa' ''}
%}
    

%{
pairings={'ABprpapppap'   'MSaappaa'       'MSpapappp'      'MSpapaaap'      'MSaappap'   
    'ABplpapaaaa' 'ABalpappppa'    'ABalpappppp'    'ABarapapppa'    'MSaapaaaa'  
    'ABprpapaaaa' 'ABarapapppp'    'ABarapapppa'    'ABarapapapp'    'ABarapaaaap'
    'ABprpapppaa' 'ABarapapppp'    'ABarapapapp'    'ABalpappapp'    'MSaappaa'   
    'ABplpappaaa' 'ABalpappppa'    'ABalpappppp'    'ABalpappapa'    'ABarapapapp'
    }
%}

%filter (as should have done on correlation (and now is)) to get equivalent time period
%of pairs
allpositions={};
for i=1:size(pairings,1)
    for k=1:size(pairings,2)
        for j=1:length(cells)
            istarget=~isempty(find(strcmp(pairings{i,k},cells{j}.name), 1));
            isaliveduringtime=cells{j}.endtime>tstart;
            
            if(istarget& isaliveduringtime)
              
                %note cells can end early by tracking error/death or
                %can start late because of nonterminal period
                %technically could be nonterminal (some are) and then divide
                %hence be mising some at end but dont expect this
                if (cells{j}.endtime-length(cells{j}.exp)>tstart)
                    
                    %cell is born during time window copy relevant portion of
                    %parent and concatenate
                    parentbit=[];
                    for psearch=1:length(cells)
                        %is parent
                        if(strcmp(cells{psearch}.name,cells{j}.name(1:end-1)))
                            %compute start point in parent
                            startindex=cells{psearch}.endtime-length(cells{psearch}.exp);%start time of j
                            startindex=tstart-startindex+1;
                            parentbit=cells{psearch}.pos(startindex:end,:);
                        end
                    end
                    
                    allpositions{i,k}=[parentbit;cells{j}.pos];
                    
                else
                    'here 3'
                    %if born at beginning of period there is some to trim off
                    startindex=cells{j}.endtime-length(cells{j}.exp);%start time of j
                    startindex=tstart-startindex+1;        
                    allpositions{i,k}=cells{j}.pos(startindex:end,:);
                end
                
            end
        end
    end
end

minlengths=[];
for i=1:size(pairings,1)
    minlengths(i)=length(allpositions{i,1});
    if (size(pairings,1)>1)
    for j=2:size(pairings,2)
        minlengths(i)=min(minlengths(i),length(allpositions{i,j}));
    end
    end
end

distances={}
for i=1:size(pairings,1)
    temp=allpositions{i,2}(1:minlengths(i),:);
    for k=3:size(pairings,2)
        temp=temp+allpositions{i,k}(1:minlengths(i),:);
    end
    temp=temp./(size(pairings,2)-1); %average position
    disvals=[];
    for j=1:minlengths(i)
        disvals(j)=distance_anisotropic(allpositions{i,1}(j,:)',temp(j,:)',[1,1,anisotropy]);
    end
    distances{i}=disvals;
end
%}
%{
%average of distance to 3 instead of distance to average of 3
distances={}
for i=1:size(pairings,1)
    temp=0;
    for k=2:size(pairings,2)
       % temp=temp+allpositions{i,k}(1:minlengths(i),:);
         temp=temp+distance_anisotropic(allpositions{i,1}(j,:)',allpositions{i,k}(1:minlengths(i),:)',[1,1,anisotropy]);
    end
    temp=temp./(size(pairings,2)-1); %average position
      distances{i}=temp;
end
%}
figure
for i=1:size(pairings,1)
    subplot(size(pairings,1),1,i);
   
    plot(distances{i});
    axis([0,70,0,50])
    hold on
    %plot(distanceswt{i},'g');
end

