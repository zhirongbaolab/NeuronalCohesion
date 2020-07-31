


%%new script to pick most contacting 4 to use in distance calculation

editedcellssuls={'ABprpapppap'
    'ABplpapaaaa'
    'ABprpapaaaa'
    'ABprpapppaa'
    'ABplpappaaa'
    
'ABalpaapaap'
'ABalpappapa'
'ABaraaapaap'
'ABaraaaaaaa'
'ABaraaaaapa'
'ABaraaaapap'
'ABalpaaaaaa'
'ABaraapaaaa'
'ABarapaaaaa'
'ABalpaapapa'
'ABalpappppa'
'ABaraaapapa'
'ABaraapaapp'
'ABaraappapp'
'ABarapapppa'
'ABalpappapp'
'MSaaapapa'
'MSpaapapa'
'MSaapaapaa'
'MSpapaapaa'
'MSaaaaapap'
'MSpapapaa'
'MSaapapaa'
'ABalpaaappp'
'ABarapaappp'
'ABalpapppaa'
'ABarapappaa'
'ABalpappaapa'
'ABarapapaapa'
'ABaraapapaav'
'ABaraapppaav'
'MSaaaaapaa'
'MSaaaapaa'
'ABarapapapp'
'MSpaaapaa'
'ABaraappaaa'
'ABalpaapapp'
'ABalpappppp'
'ABaraapappp'
'ABaraappppp'
'ABarapaappa'
'ABarapapppp'
'MSaaapaaa'
'MSpaaaapa'
'ABalpaaaaap'
'ABalpaaapaa'
'ABaraapaapa'
'ABaraapappa'
'ABaraappapa'
'ABaraappppa'
'ABarapaaaap'
'ABarapaapaa'
'ABalpaaaapa'
'ABaraaaaaap'
'ABaraaaaapp'
'ABaraapaaap'
'ABaraappaap'
'ABarapaaapa'
'MSpaapaaa'
'ABaraaapapp'
'MSaaaaapp'
'MSaaapaap'
'MSaapaaaa'
'MSpaaaaaa'
'MSpaaaapp'
'MSpapaaaa'
'ABaraapapap'
'ABaraapppap'
'MSaaaapap'
'MSaapaaap'
'MSpaaapap'
'MSpaaappa'
'MSpapaaap'
'MSaapappa'
'MSpaaappp'
'MSpapappa'
'MSaaaappp'
'MSaapaapp'
'MSpapaapp'
'MSaaapapp'
'MSaaappp'
'MSaapappp'
'MSaappaa'
'MSaappap'
'MSpaapapp'
'MSpapappp'
'ABaraappaa'
'ABaraapppp'
'ABarapaapa'
'ABarapaapap'
'ABaraapapp'
'ABaraapppp'
'ABalpaaaaa'
'ABarapaaap'
'ABaraaaaap'
'ABaraaaaaa'
'ABalpaaapa'
'ABalpaaapap'
'ABalpaaapp'
'ABalpaaappa'
'ABalpaaaap'
'MSpaaaaa'
'ABarapappa'
'MSpaaapa'
'MSpaapaa'
'ABalpapppa'
'ABalpapppap'
'ABaraaaapa'
'MSaaaaapa'
'MSpapapa'
'MSaapapa'
'ABalpappaap'
'MSaapaapa'
'ABarapapaap'
'ABaraapppaa'
'MSpapaapa'
'ABaraapapaa'
'ABaraaapaaa'
'ABalpaapaaa'

    }
%'ABalpaaaapp'
%'ABarapaaapp'
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
    playergraphs{i}=zeros(length(editedcellssuls),length(editedcellssuls));
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

totalgraph=zeros(length(editedcellssuls),length(editedcellssuls));
changes=0;
for i=tstart:tend
    totalgraph=totalgraph+playergraphs{i};    
end

neuphar=totalgraph(1:5,6:end);
mostcontacting={}
for i=1:5
    for j=1:4
        [v,ind]=max(neuphar(i,:));
        mostcontacting{i,j}=editedcellssuls{ind+5};
        neuphar(i,ind)=0;
    end
end