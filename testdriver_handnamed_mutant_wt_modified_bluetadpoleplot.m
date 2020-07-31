
%load 
tstart=245; %start of pharynx
tend=310; %end of retraction

emb='L:\santella\nih_emb_qc\emb2\Decon_emb1_MGedits.zip'; %e2

anisotropy=1;
%subset to use with neurons in plot 3
limitedindicies=[9,14,15,31,32];
pairings={
'ABplpapaaaa'	'ABalpappppp'
'ABprpapaaaa'	'ABarapapppp'
'ABalppappaa'	'ABalpaaaaaa'
'ABarappppaa'	'ABaraaaaapa'
'ABalppappap'	'ABalpaaaaaa'
'ABarappppap'	'ABarapaaaap'
'ABplpappaaa'	'ABalpappppa'
    };

editedcells={
'AIML'
'AIMR'
'AINL'
'AINR'
'AIYL'
'AIYR'
'ALA'
'AVAL'
'AVAR'
'AVDL'
'AVDR'
'AVHL'
'AVHR'
'CEMVL'
'CEMVR'
'AVL'
'RIPL'
'RIPR'
'RIVL'
'RIVR'
'RMDVL'
'RMDVR'
'RMED'
'RMEL'
'RMER'
'SAAVL'
'SAAVR'
'SIADL'
'SIADR'
'SIBVL'
'SIBVR'
'SMDDL'
'SMDDR'
'SMDVL'
'SMDVR'
'URAVL'
'URAVR'
'URYVL'
'URYVR'
'm4dl'
'm5dl'
'm7d'
'm3dl'
'm4l'
'mc3dl'
'vpi3d'
'm7vl'
'vpi2dl'
'vpi2v'
'vpi3v'
'm3dr'
'm4dr'
'm5dr'
'mc3dr'
'vpi1'
'm4vr'
'm5vr'
'g1ar'
'm7vr'
'm6vl'
'vpi2dr'
'e1d'
'e1vl'
'e3vl'
'e1vr'
'e2dl'
'e2v'
'e2dr'
'e3d'
'm1vl'
'm1dl'
'm1dr'
'm2l'
'm2vl'
'm2dl'
'i1l'
'i2l'
'm4r'
'mc1v'
'mc1dr'
'mc2dr'
'mc2dl'
'm3vl'
'mc3v'
'Posterior_arcade_dl'
'Posterior_arcade_dr'
};



load ('partlist.mat');
editedcellssuls={};
for i=1:length(editedcells)
    editedcellssuls{i}=terminalToSulston(editedcells{i},partlist);
end
%editedcellssuls{8}='MSaaaaaa';
%editedcellssuls{12}='ABalpaapaaa';
%editedcellssuls{13}='ABaraaapaaa';
%editedcellssuls{9}='ABalpaapaaa';
%editedcellssuls{10}='ABaraaapaaa';

editedcellssuls{end-1}='ABalpaapaaa';
editedcellssuls{end}='ABaraaapaaa';

cellcolors={
}
for i=1:length(editedcells)
    if i<=length(editedcells)-47
        cellcolors{i}='b';
    else
        cellcolors{i}='g';
    
    end
end
numcellcolors=[];
for i=1:length(editedcells)
    if i<=length(editedcells)-47
        numcellcolors(i,:)=[0,0,1];
    else
        numcellcolors(i,:)=[0,1,0];
    
    end
end



templocation='temp_unzip\';
%unzip zipfile to temp file
if ~exist(emb,'file')
    errors.zipfilemissing=1;
    return
end
try
    unzip(emb,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end

[ cells,embdat] = loadcells_unnamed(templocation,tend,4,false );
rmdir(templocation,'s');

[embdat_stabilized]=internallyAlignNamedEmbryo(embdat,tstart,tend,anisotropy);
[ cells_stabilized ] = parseCellsFromEmb( embdat_stabilized,tend );

cells_backup=cells;
cells=cells_stabilized;



percellmotion={};
%for all cells

%parse data into target cells 
for i=1:length (cells)
%if in time range
istarget=~isempty(find(strcmp(editedcellssuls,cells{i}.name), 1));
isaliveduringtime=cells{i}.endtime>tstart;

    if(istarget& isaliveduringtime)   
        %note cells can end early by tracking error/death or
        %can start late because of nonterminal period
        %technically could be nonterminal (some are) and then divide but 
        %dont expect this, 
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

%{

%ice cream cone vis
motions=[]
colors=[];
for i=1:length(percellmotion)
if ~isempty(percellmotion{1,i})
        motion1=percellmotion{1,i}; %skip first frames
        motion1(:,3)=motion1(:,3).*anisotropy;
      %  motion1=smoothdata(motion1);
        motions=[motions;motion1(1,:),motion1(end,:)];
        
        target=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        colors=[colors;numcellcolors(target,:)];;
%            plot3( motion1(:,1), motion1(:,2), motion1(:,3),'Color',cellcolors{target},'LineWidth',4);
%            scatter3(motion1(end,1),motion1(end,2),motion1(end,3),cellcolors{target});
        %text(motion1(end,1),motion1(end,2),motion1(end,3),editedcells{target},'FontSize',17);

end
end

figure
hold on

for i=1:400
    scatter3(cells{1,i}.pos(1,1),cells{1,i}.pos(1,2),cells{1,i}.pos(1,3)*anisotropy,'.k');
end

iceCreamConePlot(motions(:,1:3),motions(:,4:6),3,7,50,1,colors);

axis equal
axis vis3d
%}


%alternate vis to assess sorting phenotype
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
        if(strcmp(cellcolors{target},'b'))  
        %plot3( motion1(:,1), motion1(:,2), motion1(:,3),'Color',cellcolors{target},'LineWidth',4);
        %    scatter3(motion1(end,1),motion1(end,2),motion1(end,3),1500,cellcolors{target},'.');
        %text(motion1(end,1),motion1(end,2),motion1(end,3),editedcells{target},'FontSize',17);
        
        maxwidth=16;
        for j=1:length(motion1)-1
            
            widthc=maxwidth*j/(length(motion1)-1);
            if j>1
                plot3( motion1(j-1:j+1,1), motion1(j-1:j+1,2), motion1(j-1:j+1,3),'Color',cellcolors{target},'LineWidth',widthc)
                
            else
                plot3( motion1(j:j+1,1), motion1(j:j+1,2), motion1(j:j+1,3),'Color',cellcolors{target},'LineWidth',widthc)
            end
        end
        
        scatter3(motion1(end,1),motion1(end,2),motion1(end,3),3000,cellcolors{target},'.');
        
        end
end
end
axis equal
axis vis3d


%{
%vis with just neurons
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
        if target<=length(editedcells)-47
            plot3( motion1(:,1), motion1(:,2), motion1(:,3),'Color',cellcolors{target},'LineWidth',4);
            scatter3(motion1(end,1),motion1(end,2),motion1(end,3),cellcolors{target});
        end
      %text(motion1(end,1),motion1(end,2),motion1(end,3),editedcells{target},'FontSize',17);

end
end
axis equal
axis vis3d
%}
%{

%vis with just pha and small neuron listneurons
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
        if (target>length(editedcells)-47)||(~isempty(find(target==limitedindicies)))
            plot3( motion1(:,1), motion1(:,2), motion1(:,3),'Color',cellcolors{target},'LineWidth',4);
            scatter3(motion1(end,1),motion1(end,2),motion1(end,3),cellcolors{target});
        end
      %text(motion1(end,1),motion1(end,2),motion1(end,3),editedcells{target},'FontSize',17);

end
end
axis equal
axis vis3d

%}


%correlation vs graph based comparison
%truncates here which is ok if are careful above
corrv=[];
for i=1:length(percellmotion)
    if ~isempty(percellmotion{1,i})
        for j=1:length(percellmotion)
            if ~isempty(percellmotion{1,j})&&j~=i
                motion1=percellmotion{1,i}; 
                motion2=percellmotion{1,j};
                motion1=smoothdata(motion1);
                motion2=smoothdata(motion2); 
                v1=motion1(1:end-1,:)-motion1(2:end,:);
                v2=motion2(1:end-1,:)-motion2(2:end,:);
                 ind1=(find(strcmp(editedcellssuls,cells{i}.name), 1));
        ind2=(find(strcmp(editedcellssuls,cells{j}.name), 1));
        minlength=min(length(v1),length(v2));
        if (minlength~=length(v1)||minlength~=length(v2))
            ['warning truncating during calculation missing timepoints ',num2str(i),' ',num2str(length(v1)),' ',num2str(j),' ',num2str(length(v2))]
        end
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
figure
imagesc(corrv);
%xticklabels(editedcells)
yticks(linspace(1,length(editedcells),length(editedcells)));

yticklabels(editedcells)
title('Motion correlation between all pairs of cells WT')




%gabriel graph model fo entire embryo
%nothing wrong here that I know of
graphs={};
for i=tstart:tend      
         dataDist=distance_anisotropic(embdat(i).finalpoints',embdat(i).finalpoints',[1,1,anisotropy]);
        [dall,indall]=sort(dataDist,'ascend');
        %graphs{i}=GabrielGraph(embdat(i).finalpoints,indall);
        graphs{i}= DelaunayGraph(embdat(i).finalpoints);
end

% copy just the player connectivity information into a smaller graph
playergraphs={};
for i=tstart:tend
    playergraphs{i}=zeros(length(editedcells),length(editedcells));
    currbiggraph=graphs{i};
    for j=1:length(editedcellssuls)
        %find this target in cell list
        ind=find(strcmp(embdat(i).cellnames,editedcellssuls{j}));
        
        %if found find all the rest in its row and copy their connectivity
        %into mini matrix
        if ~isempty(ind)
            for k=1:length(editedcellssuls)
                if k~=j
                    [ind2]=find(strcmp(embdat(i).cellnames,editedcellssuls{k}));
                    if ~isempty(ind2)
                        %big graph is a neighbor list not  matrix
                        if~isempty(find(currbiggraph{ind}==ind2))
                            playergraphs{i}(j,k)=1;
                        end
                    end
                end
            end
        end
    end
end
%at this point have copied out all the contacts
%{
%vis total matrix
figure
c=1;
num2disp=length(editedcells);
for i=1:num2disp%length(editedcells)
    for j=1:num2disp%length(editedcells)
        if(i>j)
        %subplot(length(editedcells),length(editedcells),c);
        subplot(num2disp,num2disp,c);
        
        currcontact=[];
        for t=tstart:tend
            currcontact=[currcontact;playergraphs{t}(i,j)];
          
        end
          plot(currcontact);
          axis off
        end
        c=c+1;
        
    end
end
%}

totalgraph=zeros(length(editedcells),length(editedcells));
totalchanges=zeros(length(editedcells),length(editedcells));
changes=0;
for i=tstart:tend
    totalgraph=totalgraph+playergraphs{i};
    if (i~=tend)
        changes=changes+sum(sum(abs(playergraphs{i}-playergraphs{i+1})));
        totalchanges=totalchanges+abs(playergraphs{i}-playergraphs{i+1});
    end
    
end
return

%this block is fine and correct, but we're phasing out this figure
%{
figure;imagesc(totalgraph);title('total contacts WT')
%xticklabels(editedcells)
yticks(linspace(1,length(editedcells),length(editedcells)));
yticklabels(editedcells)
figure;imagesc(totalchanges);title('total contact changes WT');
%xticklabels(editedcells)
yticklabels(editedcells)
yticks(linspace(1,length(editedcells),length(editedcells)));


subtotaltotalchanges=sum(totalchanges(1:length(editedcells)-7,:));
barh(subtotaltotalchanges)
title('Total changed in neuron contacts for neurons and pharynx');
yticks(linspace(1,length(editedcells),length(editedcells)));
yticklabels(editedcells)
%}

%new analysis showing positions of paired pharynx neuron cells over time
%number of cells and plotting of wt next to mutuant is hard coded here
%pairings is not loaded from disk
%trimming here does not account for parent-childe quivalence which might be
%desirable 

%now uses clipping and assembly used in correlation 
%{
for i=1:length(pairings)
    pairings{i,2}=terminalToSulston(pairings{i,2},partlist);
end
%}
allpositions={};
for i=1:length(pairings)

    for j=1:length(cells)
        istarget=~isempty(find(strcmp(pairings{i,1},cells{j}.name), 1));
        %position=find(strcmp(cells.n
        if(istarget)
%            startpos=max(1,(cells{j}.endtime-length(cells{j}.pos)-tstart+1));
%            endpos=length(cells{j}.pos);
%            allpositions{i,1}=cells{j}.pos(startpos:endpos,:);
            allpositions{i,1}=percellmotion{j};
        end
        istarget2=~isempty(find(strcmp(pairings{i,2},cells{j}.name), 1));
        if (istarget2)
%              startpos=max(1,(cells{j}.endtime-length(cells{j}.pos)-tstart+1));
%              endpos=length(cells{j}.pos);
%             allpositions{i,2}=cells{j}.pos(startpos:endpos,:);
              allpositions{i,2}=percellmotion{j};
    
        end
    end
end

distances={}
for i=1:length(pairings)
minlen=min(length(allpositions{i,1}),length(allpositions{i,2}));
disvals=[];
for j=1:minlen
    if (minlen>3)
     disvals(j)=distance_anisotropic(allpositions{i,1}(j,:)',allpositions{i,2}(j,:)',[1,1,anisotropy]);
    end
end
distances{i}=disvals;
end

figure
for i=1:length(pairings)
    subplot(7,7,i);
   
    plot(distances{i});
    %hold on
    %plot(distanceswt{i},'g');
end

endstarts=[];
for i=1:length(pairings)
    if(length(distances{i})>0)
    endstarts(i)=distances{i}(end)-distances{i}(1);
    end
end
endstarts=endstarts';


%just want to pull out neuron neuron and neu pharynx
return


totaltotalchanges=sum(totalchanges);
figure;imagesc(totaltotalchanges')
yticks(linspace(1,length(editedcells),length(editedcells)));
yticklabels(editedcells);
title 'total changes per cell wt';

figure
barh(totaltotalchanges)
yticks(linspace(1,length(editedcells),length(editedcells)));
yticklabels(editedcells)
axis([0,120,0,length(editedcells)+1])
title 'total changes per cell wt';



return
%old stuff not doing it now
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
        istarget=~isempty(find(strcmp(pharynxlist,cells{i}.name), 1));
        if istarget
            plot3( motion1(:,1), motion1(:,2), motion1(:,3),'Color','g','LineWidth',4)
            scatter3(motion1(end,1),motion1(end,2),motion1(end,3),'g');
        else
             plot3( motion1(:,1), motion1(:,2), motion1(:,3),'Color','b','LineWidth',4)
            scatter3(motion1(end,1),motion1(end,2),motion1(end,3),'b');
        end
end
end
axis equal

%quantitation
retractiontarget=editedcellssuls(8:10);%'e3d'
involutiontarget=editedcellssuls(1:7);%'ala'
retraction=[];
involution=[];
for i=1:length (cells)
    if(~isempty(find(strcmp(retractiontarget,cells{i}.name))))
        retraction=[retraction;distance_anisotropic(cells{i}.pos(1,:)',cells{i}.pos(end,:)',[1,1,anisotropy])];
    end
        if(~isempty(find(strcmp(involutiontarget,cells{i}.name))))
        involution=[involution;distance_anisotropic(cells{i}.pos(1,:)',cells{i}.pos(end,:)',[1,1,anisotropy])];
    end
end

