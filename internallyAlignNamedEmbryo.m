function [emb_e1stab]=internallyAlignNamedEmbryo(emb_e1,temb1,t2emb1,anisotropy)
%I wrote this but have not debugged it because I realize that to use it to
%visualize cell tracks in the way I wanted I'd have two write a parser to
%pull out cell data structure like version from emb data structure and I
%dont have the patience right now 
emb_e1stab=emb_e1;
%stabilize embryo affinely over time period provided
for time=temb1:t2emb1-1
    
names_1=emb_e1(time).names;
names_2=emb_e1(time+1).names;

%pos_1=emb_e1(time).finalpoints;
%aligning t+1 to t so t needs to be laready stabilized
pos_1=emb_e1stab(time).finalpoints;

pos_2=emb_e1(time+1).finalpoints;

c=1;
    matchpoint1=[];
    matchpoint2=[];
for i=1:length(names_1)

    for j=1:length(names_2)
        if strcmp(names_1{i},names_2{j})&(isempty(strfind(names_1{i},'Nuc')))
            indmatch1=j;
            matchpoint1=[matchpoint1;pos_1(i,:)];
            matchpoint2=[matchpoint2;pos_2(j,:)];
        end
    end
   %{
     if (~isempty(matchpoint1)&~isempty(matchpoint2))
           lmpositions1=[lmpositions1;matchpoint1];
        lmpositions2=[lmpositions2;matchpoint2];
    end
    %}
end
lmpositions1=matchpoint1;
lmpositions2=matchpoint2;
%compute transformation

lmpositions1(:,3)=lmpositions1(:,3);
lmpositions2(:,3)=lmpositions2(:,3);

transform2to1=[lmpositions1,ones(length(lmpositions1),1)]'/[lmpositions2,ones(length(lmpositions2),1)]';
%alltransforms{e}=transform2to1;

emb_e1stab(time+1).finalpoints=(transform2to1*[emb_e1(time+1).finalpoints,ones(size(emb_e1(time+1).finalpoints,1),1)]')';
emb_e1stab(time+1).finalpoints=emb_e1stab(time+1).finalpoints(:,1:3);
end