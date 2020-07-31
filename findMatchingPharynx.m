function [percellmotion,closeststart]=findMatchingPharynx(pos,tstart,endtime,pharynxlist,emb)
    
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

percellmotion={};
%percellmotion{1}=pos;

pharmotion=[];
posmatch=[];
%start of cell or start of analysis window whichever comes first
closeststart='';
for i=max(tstart,endtime-size(pos,1)):endtime-1
    index=i-(endtime-size(pos,1))+1;
    
    if(isempty(closeststart))
        closeststart=findclosestpharynxname(pos(index,:),pharynxlist,emb(i));
    end
    curpharpos=findmatch(closeststart,emb(i));
  
    %if didnt find our original match, rematch and try again
    if isempty(curpharpos)
            closeststart=findclosestpharynxname(pos(index,:),pharynxlist,emb(i));
             curpharpos=findmatch(closeststart,emb(i));
    end
    
    if ~isempty(curpharpos)
        pharmotion=[pharmotion;curpharpos];
   
    %pharmotion=[pharmotion;findclosestpharynx(pos(index,:),pharynxlist,emb(i),15)];
        posmatch=[posmatch;pos(index,:)];
    end
end
%percellmotion{2}=pharmotion;
percellmotion=[posmatch,pharmotion];
end

function pos=findmatch(name,embi)
ind=find(strcmp(embi.cellnames,name));
pos=embi.finalpoints(ind,:);
end

%note want to use multi point average
function name=findclosestpharynxname(cellpos,pharynxlist,embi)
    [pharynxcells,names]=findcells(embi,pharynxlist);
    distances=distance(cellpos',pharynxcells');
        [v,ind]=min(distances);
    name=names{ind};
end

%note want to use multi point average
function position=findclosestpharynx(cellpos,pharynxlist,embi,num)
    pharynxcells=findcells(embi,pharynxlist);
    distances=distance(cellpos',pharynxcells');
    position=[];
    for i=1:num
        [v,ind]=min(distances);
        position=[position;pharynxcells(ind,:)];
        distances(ind,:)=inf;
    end
    
        if(size(position,1)>1)
            position=mean(position);
        end
    end

function [positions,names]=findcells(embi,celllist)
    positions=[];
    names={};
for i=1:length(embi.cellnames)
    if ~isempty(find(strcmp(celllist,embi.cellnames{i}), 1))
        positions=[positions;embi.finalpoints(i,:)];
        names{end+1}=embi.cellnames{i};
    end
end
end