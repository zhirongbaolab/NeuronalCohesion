
%load 
%embryo 1 with lots of edited cells
%{
tstart= 250;%285;%correct one %250;start of pharynx
tend= 350;%315;%%; %end of retraction
emb='L:\KB\Head_Formation\Mutant+RNAi_Strains\hmr-1\BV731\10232019_lineaging\10232019_newestedits\KB_BV731_10232020_w2iSIM-TxRed-600-50_s2_emb1_edited.zip'; %e2
anisotropy=1/.16;
%}
%lisf for fig 3 revision


%emb='L:\KB\Head_Formation\Mutant+RNAi_Strains\hmr-1\BV731\10232019_lineaging\10232019_lineaging\KB_BV731_10232020_w2iSIM-TxRed-600-50_s1_emb1_edited.zip';
%{
emb='D:\kris_neuronscreen\figures_3_31_2020\KB_BV731_10232020_w2iSIM-TxRed-600-50_s1_emb1_edited.zip';

tstart=300;
tend=365;
anisotropy=1/.16;
%}
%}
%{
bad
bad 
emb='L:\KB\Head_Formation\Mutant+RNAi_Strains\hmr-1\BV731\11152019_lineaging\11152019_lineaging\KB_BV731_11152019_w2iSIM-TxRed-600-50_s2_emb1_edited.zip'
tstart=235; 
tend= 300; 
bad
%}

%new old 3d one
emb='L:\KB\Head_Formation\Mutant+RNAi_Strains\hmr-1\BV727\09082019\09082019\KB_BV727_hmr1_histones_09082019_w1iSIM-TxRed-600-50_s1_emb1_edited.zip';
tstart=190; %start of pharynx
tend=223; %end of retraction
anisotropy=1/.16;
%}

%emb='L:\santella\nih_emb_qc\new_new_version\Decon_emb1_edited_v13_11202018namingfix.zip'; %e1
%tstart=250;tend=330;
%anisotropy=1;
%emb='L:\santella\nih_emb_qc\emb3\Decon_emb1_updated.zip'; %e3
%tstart=225;tend=305;%
%emb='L:\santella\nih_emb_qc\emb2\Decon_emb1_MGedits.zip'; %e2
%tstart=245;tend=310;
%anisotropy=1;

%limitedindicies=[9,14,15,31,32];
templocation='temp_unzip2\';
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

%}

%'ABplpaapppa'
%'ABprpaapppa'
%{
'ABalpappppp'
'ABarapapppp'
'ABalpaaaaaa'
'ABaraaaaapa'
'ABarapaaaap'
'ABalpappppa'
%}

%chunk list of cells
editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABplpapaapp'
'ABprpapaapp'
'ABplpapaapa'
'ABprpapaapa'
'ABplpapaaap'
'ABprpapaaap'
'ABplpaapppp'
'ABprpaapppp'
};
%'ABplpappaaa'
editedcells={
'SMDDL'
'SMDDR'
'SIBVL'
'SIBVR'
'SIADL'
'SIADR'
'AIYL'
'AIYR'
'ABplpaapppp'
'ABprpaapppp'
};


%chunk list of cells without smddl for fig 2
editedcellssuls={
'ABplpapaapp'
'ABprpapaapp'
'ABplpapaapa'
'ABprpapaapa'
'ABplpapaaap'
'ABprpapaaap'
'ABplpaapppp'
'ABprpaapppp'
};
%'ABplpappaaa'
editedcells={
'SIBVL'
'SIBVR'
'SIADL'
'SIADR'
'AIYL'
'AIYR'
'ABplpaapppp'
'ABprpaapppp'
};

cellcolors={'g','g','b','b','r','r','m','m','c','c','k','k',[.5,.5,.5],[.5,.5,.5]};


%{
%limited list for fig 3 subset (mutant version that is no lonter used
editedcells={
'SMDDL'
'SMDDR'
'RIS'
'RMEV'
'AVG'

'e2v'
'mc1v'
'mc3v'
'posterior_arcade V'
'm5vl'
'm5vr'
};
editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABprpappapa'
'ABplpappaaa'
'ABprpapppap'

'ABalpappapa'
'ABalpappppa'
'ABalpappapp'
'ABarapapapa'
'MSaapaaap'
'MSpapaaap'
};
cellcolors={[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],'g','g','g','g','g','g'};
%}
%variant for wt using original adherence screen based pharynx list
%limited list for fig 2 subset
editedcells={
'SMDDL'
'SMDDR'
'RIS'
'RMEV'
%'AVG'

'p'
'p'
'p'
'p'
'p'

};
editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABprpappapa'
'ABplpappaaa'
%'ABprpapppap'
     
'ABarapapppa'
'ABarapapppp'    
'ABarapapapp'
'ABalpappppa'    
'ABalpappppp' 
};
cellcolors={[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],'g','g','g','g','g','g'};

%cellcolors={[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],'g','g','g','g','g','g'};


%new new for mutant revised figure 3
%{
%new version 

editedcells={
'SMDDL'
'SMDDR'
'RIS'
'RMEV'


'mc1dr'
'm6vl'
'm7vl'
'mc2dl'
'mc2dr'
};

editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABprpappapa'
'ABplpappaaa'


'ABaraaapapa'
'MSaapappa'
'MSaapaapp'
'ABaraapaapp'
'ABaraappapp'
};
cellcolors={[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],'g','g','g','g','g','g','g','g'};




%old version
editedcells={
'SMDDL'
'SMDDR'
'RIS'
'RMEV'
'AVG'

'mc3dl'
'mc3dr'
'mc1dr'
'm6vl'
'm7vl'
'mc2dl'
'mc2dr'
};

editedcellssuls={
'ABplpapaaaa'
'ABprpapaaaa'
'ABprpappapa'
'ABplpappaaa'
'ABprpapppap'

'MSaaapapa'
'MSpaapapa'
'ABaraaapapa'
'MSaapappa'
'MSaapaapp'
'ABaraapaapp'
'ABaraappapp'
};
cellcolors={[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],[1,.5,0],'g','g','g','g','g','g','g','g'};



%}



   %{
   editedcellssuls={
   'ABplpapaaaa'
   'ABprpapaaaa'
   'ABplpappaaa'
   
   'ABalpappppa'
   
   'ABalpappppp'
   'ABarapapppa'
   'ABarapapppp'
   'ABarapapppa'
   'ABarapapapp'
   'ABalpappppa'
   'ABalpappppp'
   'ABarapapapp'
   };
   

editedcells={
'SMDDL'
'SMDDR'
'RMED'
'pha'
'pha'
'pha'
'pha'
'pha'
'pha'
'pha'
'pha'
};
%'RMEV'

%manual list for mutant
 %{
   editedcellssuls={
   'ABplpapaaaa'
   'ABprpapaaaa'
  'ABaraappapp'    
'ABalpappapa'    
'ABaraapaaap'
 'MSpapaaap'
    'ABarapapppp'
   'MSpapaaaa'
   };
   
cellcolors={'g','g','b','b','b','b','b','b','b','b','b','b'};

  
editedcells={
   'ABplpapaaaa'
   'ABprpapaaaa'
  'ABaraappapp'    
'ABalpappapa'    
'ABaraapaaap'
 'MSpapaaap'
    'ABarapapppp'
   'MSpapaaaa'
};

manual list for wt distanceovertime
editedcellssuls={'ABplpapaaaa'
 'ABprpapaaaa'
'ABalpappppa'
    'ABalpappppp'    
'ABarapapppa'
   'ABarapapppp'    
'ABarapapppa'    
'ABarapapapp'
}
%}
%}



%cellcolors={'g','g','b','b','r','r','m','m','c','c','k','k',[.5,.5,.5],[.5,.5,.5]};

%cellcolors={'g','g','g','b','b','b','b','b','b','b','b','b'};

  
    
    
%{
%short list of leaders with wide pharynx
editedcells={
'AVG'
'SMDDL'
'SMDDR'
'RIR'
'RMEV'

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
editedcellssuls{end-1}='ABalpaapaaa';
editedcellssuls{end}='ABaraaapaaa';

}%
%{
editedcells={
'SMDDL'
'SMDDR'
'SIBVL'
'SIBVR'
'AIML'
'AIMR'
'SIADL'
'SIADR'
'AIYL'
'AIYR'
'ABplpaapppp'
'ABprpaapppp'
'RMEV'

'm3vl'
'm3vr'
'e3vl'
'e1vr'
'e3vl'
'm2vr_of_pm2'
'mc1v'
};
%}
pairings={
'ABplpapaaaa'	'ABalpappppp'
'ABprpapaaaa'	'ABarapapppp'
'ABalppappaa'	'ABalpaaaaaa'
'ABarappppaa'	'ABaraaaaapa'
'ABalppappap'	'ABalpaaaaaa'
'ABarappppap'	'ABarapaaaap'
'ABplpappaaa'	'ABalpappppa'
    };


load ('partlist.mat');
%{
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
%}
cellcolors={
}
for i=1:length(editedcellssuls)
    if i<=length(editedcellssuls)-7
        cellcolors{i}='b';
    else
        cellcolors{i}='g';
    
    end
end

%}


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

%{
diststarts=[];
distends=[];
%new analysis of convergence of l/r pairs over window

for i=1:2:length(editedcells)
    ind1=[];ind2=[];
    for j=1:length(cells)
        if(strcmp(editedcellssuls{i},cells{j}.name))
            ind1=j;
        end
        if(strcmp(editedcellssuls{i+1},cells{j}.name))
            ind2=j;
        end
      
    end
      if ~isempty(ind1)&&~isempty(ind2)
            diststarts=[diststarts;distance_anisotropic(percellmotion{1,ind1}(1,:)',percellmotion{1,ind2}(1,:)',[1,1,anisotropy])];
            distends=[distends;distance_anisotropic(percellmotion{1,ind1}(end,:)',percellmotion{1,ind2}(end,:)',[1,1,anisotropy])];
                
        end
end

return 
%}
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
corrv=[];
for i=1:length(percellmotion)
    if ~isempty(percellmotion{1,i})
        for j=1:length(percellmotion)
            if ~isempty(percellmotion{1,j})&&j~=i
                motion1=percellmotion{1,i}; 
                motion2=percellmotion{1,j};
                %note this was somehow missing in the original version 
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
figure
imagesc(corrv);
%xticklabels(editedcells)
yticks(linspace(1,length(editedcells),length(editedcells)));

yticklabels(editedcells)
title('Motion correlation between all pairs of cells mutant')
return

%{
%gabriel graph model fo entire embryo
graphs={};
for i=tstart:tend
    dataDist=distance_anisotropic(embdat(i).finalpoints',embdat(i).finalpoints',[1,1,anisotropy]);
    [dall,indall]=sort(dataDist,'ascend');
    tempdat=embdat(i).finalpoints;
    %hack to permute repeated points not clear why this wasnt a problem
    %earlier
    for k=1:length(tempdat)
        for j=k+1:length(tempdat)
            if(k~=j)&&tempdat(k,1)==tempdat(j,1)&&tempdat(k,2)==tempdat(j,2)&&tempdat(k,3)==tempdat(j,3)
                tempdat(k,3)=tempdat(k,3)+.1*rand(1,1);
            end
        end
    end
    %graphs{i}=GabrielGraph(embdat(i).finalpoints,indall);
 %   graphs{i}= DelaunayGraph(embdat(i).finalpoints);
    graphs{i}= DelaunayGraph(tempdat);
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
figure;imagesc(totalgraph);title('total contacts mutant')
%xticklabels(editedcells)
yticks(linspace(1,length(editedcells),length(editedcells)));
yticklabels(editedcells)
figure;imagesc(totalchanges);title('total contact changes mutant');
%xticklabels(editedcells)
yticklabels(editedcells)
yticks(linspace(1,length(editedcells),length(editedcells)));

subtotaltotalchanges=sum(totalchanges(1:length(editedcells)-7,:));
barh(subtotaltotalchanges)
title('Total changed in neuron contacts for neurons and pharynx');
yticks(linspace(1,length(editedcells),length(editedcells)));
yticklabels(editedcells)
%}

%new analysis compute distance over time of pairs
%{
for i=1:length(pairings)
    pairings{i,2}=terminalToSulston(pairings{i,2},partlist);
end
%}
%{
for i=1:length(pairings)
    pairings{i,1}=terminalToSulston(pairings{i,1},partlist);
end
%}




return

%code for hacked all neuron start end distances
distanceovertime=[];
for i=1:length(pairings)
    if(~isempty(allpositions{i,1}))
     distanceovertime(i)=distance_anisotropic(allpositions{i,1}(1,:)',allpositions{i,1}(end,:)',[1,1,anisotropy]);
    end
end

%try to do 
return


totaltotalchanges=sum(totalchanges);
figure;imagesc(totaltotalchanges')
yticks(linspace(1,length(editedcells),length(editedcells)));
yticklabels(editedcells);
title 'total changes per cell mutant';

figure
barh(totaltotalchanges)
yticks(linspace(1,length(editedcells),length(editedcells)));
yticklabels(editedcells)
axis([0,120,0,length(editedcells)+1])
title 'total changes per cell mutant';



return

return

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
