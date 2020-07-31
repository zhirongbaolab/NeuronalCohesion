
%this script runs leader analysis based on computing voranoi with
%testcomputation_gabrielgraph then  summing up per tp graph
%and looking at sum

tp=245; 300; 325; 290; 300; %325;%320;%
%34;%good early tp 8 cell
%currnames=constraintnamesunwindowed;
%currconstraints=allcunwindowed{tp};
%currnames=allconstraints(tp).constraintnames;
%currconstraints=allconstraints(tp).neighborconstraints;
   
currnames=constraintnamesunwindowed;
%currconstraints=allcunwindowed{tp-245+1};
win=64;

currconstraints=zeros(size(allcunwindowed{tp-245+1}));
for i=tp:tp+win
   % currconstraints=min(currconstraints,allcunwindowed{i-245+1});
   currconstraints=currconstraints+allcunwindowed{i-245+1};

end


allneuronssuls={};
for i=1:length(neuroncells)
    allneuronssuls{i}=terminalToSulston(neuroncells{i},partlist);
end
%pull out connectivity matrix
pharynxbyneuron=-1*ones(length(pharynxcells),length(neuroncells));
for i=1:length(pharynxcells)
    indi=find(strcmp(currnames, pharynxcells{i}));
    if ~isempty(indi)
        for j=1:length(neuroncells)
            indj=find(strcmp(currnames,allneuronssuls{j}));
            if ~isempty(indj)
                pharynxbyneuron(i,j)=currconstraints(indi,indj);
            end
        end
    end
end

hits=find(pharynxbyneuron>30);
[i,j]=ind2sub(size(pharynxbyneuron),hits);

pharynxhits=pharynxcells(i);
neuronhits=neuroncells(j');

%find missing cells
hitsj=find(max(pharynxbyneuron)==-1);
hitsi=find(max(pharynxbyneuron')==-1);

pharynxmissing=(pharynxcells(hitsi));
neuronmissing=(neuroncells(hitsj'));

return

tp=tp+round(win/2);


currconstraints(currconstraints<(win/2))=0;

alpha1=((tp)-t4emb1)/(temb1-t4emb1);
temb2match=round(t4emb2+(temb2-t4emb2)*alpha1);

currpos=[];
for i=1:length(currnames)
    ind=find(strcmp(currnames{i},emb_e2(temb2match).cellnames));
    offset=0;
    while (isempty(ind)&offset<13)
         ind=find(strcmp(currnames{i},emb_e2(temb2match+offset).cellnames));
         if(isempty(ind))
         offset=offset+1;
         end
    end
    if(isempty(ind))
        offset=0;
    while (isempty(ind)&offset>-13)
         ind=find(strcmp(currnames{i},emb_e2(temb2match+offset).cellnames));
         if(isempty(ind))
         offset=offset-1;
         end
    end
    end
    if(isempty(ind))
        %these are bug cases that dont occur in all embryos and should be 0
        %but are 1
        currconstraints(i,:)=0;
        currconstraints(:,i)=0;
        currpos(i,:)=[-1,-1,-1];
        
    else
    currpos(i,:)=emb_e2(temb2match+offset).finalpoints(ind,:);
    end
end

earlygraph=graph(currconstraints,currnames);

LWidths = 5*(earlygraph.Edges.Weight-min(earlygraph.Edges.Weight)+1)/(max(earlygraph.Edges.Weight)-min(earlygraph.Edges.Weight));

figure
P=plot(earlygraph,'LineWidth',LWidths)

P.MarkerSize=20;
%P.LineWidth=4;
P.XData=currpos(:,2);
P.YData=currpos(:,1);
P.ZData=currpos(:,3);
P.NodeLabel = '';


%{
nl=currnames;
xd = get(P, 'XData');
yd = get(P, 'YData');
zd =get(P, 'ZData');
text(xd, yd,zd, nl, 'FontSize',20, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')
%}

%now color
%color rules
%rules={  'phar',1,'r';'intestin',1,'g';'neuron',1,'b';'body wall',1,'c'};
pharynxcells=importdata('L:\santella\kris_neuronscreen\pharynx_cells.csv',',');
for i=1:length(pharynxcells)
pharynxcells{i}=strtrim(strrep(pharynxcells{i},'"',''));
end
neuroncells=importdata('L:\santella\kris_neuronscreen\neuron_cells.csv',',');

rules={};
for i=1:length (pharynxcells)
    rules{end+1,1}=pharynxcells{i,1};
    rules{end,2}=2;
    rules{end,3}='r';
end
for i=1:length (neuroncells)
    rules{end+1,1}=neuroncells{i};
    rules{end,2}=1;
    rules{end,3}='b';
end
%{
rules{end+1,1}='RMED';
    rules{end,2}=1;
    rules{end,3}='g';
%}
   rules{end+1,1}='ALA';
    rules{end,2}=1;
    rules{end,3}='g'; 
    %{
       rules{end+1,1}='RMEL';
    rules{end,2}=1;
    rules{end,3}='g'; 
           rules{end+1,1}='RMER';
    rules{end,2}=1;
    rules{end,3}='g'; 
           rules{end+1,1}='AVDL';
    rules{end,2}=1;
    rules{end,3}='g'; 
           rules{end+1,1}='AVDR';
    rules{end,2}=1;
    rules{end,3}='g'; 
    %}
%color by rule list
for i=1:size(rules,1)
    nodes=[];
    for j=1:length(currnames)
   %     ind=find(~isempty(strfind(partlist(:,2),currnames{j})));
        ind=find((strcmp(lower(currnames{j}),lower(partlist(:,2)))));
   
        if length(ind)>1
            ind=ind(1);
        end
        if(~isempty(ind))
            matching=false;
            if(rules{i,2}==1)
                matching=~isempty(strfind(lower(partlist{ind,1}),lower(rules{i,1})));
            end
            if (rules{i,2}==3)
                matching=~isempty(strfind(lower(partlist{ind,3}),lower(rules{i,1})));
            end
             if (rules{i,2}==2)
                matching=~isempty(strfind(lower(partlist{ind,2}),lower(rules{i,1})));
            end
            
            if (matching)
                nodes=[nodes,j];
            end
        end
    end
    highlight(P,nodes,'NodeColor',rules{i,3});
end


axis equal

axis off
ax=gca;
ax.Clipping = 'off';    % turn clipping off
axis vis3d
