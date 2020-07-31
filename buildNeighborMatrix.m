function [constraintnames,allc]=buildNeighborMatrix(allnamesa,allposa)
%takes a cell array of name cell arrays one per embryo and a positions
%arraay of matricies
%all embryos are aligned already  and returned is a constraint matrix where
%entry i,j asks is i anterior of j in all data sets etc

%extract list of unique cell names for all embryos
constraintnames=allnamesa{1};
for i=2:length(allnamesa)
    constraintnames=unique({allnamesa{i}{:},constraintnames{:}});
end



allc={};
for j=1:length(allposa)
    c=zeros(length(constraintnames),length(constraintnames));

     dataDist=distFast(allposa{j},allposa{j});
    [dall,indall]=sort(dataDist,'ascend');
   % grapha=GabrielGraph(allposa{j},indall);
    grapha=DelaunayGraph(allposa{j});
   
    for h=1:length(constraintnames) %iterate over list of all names
        for i=1:length(constraintnames) %iterate over list of all names
            
            
            %find cell h if in this emb
            cellindh=find(strcmp(constraintnames(h),allnamesa{j}));
            
            %find cell i if in this emb
            cellindi=find(strcmp(constraintnames(i),allnamesa{j}));
            
            if ~isempty(cellindh)&&~isempty(cellindi)
                neighbors=grapha{cellindh};
                
                if(~isempty(find(neighbors==cellindi, 1)))
                 %previously this checked if it was empty and set to 0 if
                 %was
                    c(h,i)=c(h,i)+1;
                end
            end
        end
    end
    allc{j}=c;
    j
end



