
%not this script (no just below but belwo that is whats used to compute
%supplemental fig neuron/pharynx adjacencies over time

%assumes main testdriver_handnamed_mutant_wt.m has been run
%to compute voranois on top of this coallates neuron-pharyx adjacency
load('allneurons.mat');
load('allpharynx.mat');
perneuroncontact=zeros(length(allneurons),tend-tstart+1);
failurecasetallyn=zeros(length(allneurons),tend-tstart+1);

failurecasetallyp=zeros(length(allneurons),tend-tstart+1);
%for every neuron for every tp extract contact graph
for i=1:length(allneurons)
    for t=tstart:tend
        tnames=embdat(t).names;
        indneu=find(strcmp(tnames,allneurons{i}));
        %if not found try for its parent which should be treated as
        %equivalent
        if isempty(indneu)
            indneu=find(strcmp(tnames,allneurons{i}(1:end-1)));
        end
        
        
          if (isempty(indneu))
                failurecasetallyn(i,t-tstart+1)=failurecasetallyn(i,t-tstart+1)+1;
                end
        
        contactmade=0;
        for k=1:length (allpharynx)
            indpha=find(strcmp(tnames,allpharynx{k}));
            %if not found try for its parent which should be treated as
            %equivalent
            if isempty(indpha)
                indpha=find(strcmp(tnames,allpharynx{k}(1:end-1)));
            end
            if ~isempty(indneu)&&~isempty(indpha)&&~isempty(find(graphs{t}{indpha}==indneu,1))
                contactmade=contactmade+1;
            else
                if (isempty(indpha))
                failurecasetallyp(i,t-tstart+1)=failurecasetallyp(i,t-tstart+1)+1;
                end
            end
            
        end
        
        perneuroncontact(i,t)=contactmade;
        
%        failurecasetally(i,t)=failurecase;
    end
end

realcontact=perneuroncontact(:,tstart:tend);


return


%below is code to compute supplemental figure of contacts over time

perneuroncontact=zeros(length(allneurons),tend-tstart+1);
failurecasetallyn=zeros(length(allneurons),tend-tstart+1);
failurecasetallyp=zeros(length(allneurons),tend-tstart+1);
%for every neuron for every tp extract contact graph
for i=1:length(allneurons)
    for t=tstart:tend
        tnames=embdat(t).names;
        indneu=find(strcmp(tnames,allneurons{i}));
        %if not found try for its parent which should be treated as
        %equivalent
        if isempty(indneu)
            indneu=find(strcmp(tnames,allneurons{i}(1:end-1)));
        end
        if (isempty(indneu))
            failurecasetallyn(i,t-tstart+1)=failurecasetallyn(i,t-tstart+1)+1;
        end
        if t==tstart % set pha list
            cphalist={};
            contactmade=0;
            for k=1:length (allpharynx)
                indpha=find(strcmp(tnames,allpharynx{k}));
                %if not found try for its parent which should be treated as
                %equivalent
                if isempty(indpha)
                    indpha=find(strcmp(tnames,allpharynx{k}(1:end-1)));
                end
                if ~isempty(indneu)&&~isempty(indpha)&&~isempty(find(graphs{t}{indpha}==indneu,1))
                    cphalist={cphalist{:},allpharynx{k}};%contactmade=contactmade+1;
                end 
            end
        else %on others do actualtest
            contactmade=0;
            for k=1:length (cphalist)
                indpha=find(strcmp(tnames,cphalist{k}));
                %if not found try for its parent which should be treated as
                %equivalent
                if isempty(indpha)
                    indpha=find(strcmp(tnames,cphalist{k}(1:end-1)));
                end
                if ~isempty(indneu)&&~isempty(indpha)&&~isempty(find(graphs{t}{indpha}==indneu,1))
                   contactmade=contactmade+1;
                else
                    if (isempty(indpha))
                        failurecasetallyp(i,t-tstart+1)=failurecasetallyp(i,t-tstart+1)+1;
                    end
                end
                
            end
              perneuroncontact(i,t)=contactmade;
        end
    end
end
realcontact=perneuroncontact(:,tstart:tend);
figure; imagesc(realcontact)
return

phaterm={};
for i=1:length(allpharynx)
    phaterm{i}=SulstontoTerminal(allpharynx{i},partlist);
end


neuterm={};
for i=1:length(allneurons)
    neuterm{i}=SulstontoTerminal(allneurons{i},partlist);
end

neuterm=SulstontoTerminal(allneurons,partlist);

phaterm=SulstontoTerminal(allpharynx,partlist);
%test finding pharynx
pfoundtally=[];
for t=tstart:tend
    tnames=embdat(t).names;
    
    for k=1:length (allpharynx)
        indpha=find(strcmp(tnames,allpharynx{k}));
        %if not found try for its parent which should be treated as
        %equivalent
        if isempty(indpha)
            indpha=find(strcmp(tnames,allpharynx{k}(1:end-1)));
        end
        if ~isempty(indpha)
            pfoundtally(t,k)=1;
        end
    end
end
