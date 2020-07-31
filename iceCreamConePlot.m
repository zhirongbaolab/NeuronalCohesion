function [allhandles,lighthandle] = iceCreamConePlot( points1, points2,s1,s2 ,sampling,anisotropy,colors)
%3d plot of correspndence between 2 sets of points drawn as tapered cone
%figure
%hold off
%steps=15;
steps=1;
[x,y,z] = sphere(sampling);
c=1;
allhandles=[];
for i=1:size(points1,1)
    for j=0:steps
        alpha=j/steps;
        newcenter=(points1(i,:).*alpha)+((1-alpha).*points2(i,:));
        newcenter(3)=newcenter(3).*anisotropy;
        sizecurr= alpha*s1+(1-alpha)*s2;
        
       
        if(exist('colors'))
            h=surf((x*sizecurr+newcenter(1)),(y*sizecurr+newcenter(2)),(z*sizecurr+newcenter(3)),'FaceColor',colors(i,:),'FaceLighting','gouraud','AmbientStrength',1,'LineStyle','none');
        else
            h=surf((x*sizecurr+newcenter(1)),(y*sizecurr+newcenter(2)),(z*sizecurr+newcenter(3)));
        end
        hold on
        allhandles(c)=h;
    c=c+1;
    end
    tempp1=points1(i,:);
    tempp2=points2(i,:);
    tempp1(3)=tempp1(3).*anisotropy;
    tempp2(3)=tempp2(3).*anisotropy;
    if(exist('colors'))
        
    Cone(tempp1,tempp2,[s1,s2],sampling,colors(i,:),0,0);
  else
    Cone(tempp1,tempp2,[s1,s2],sampling,'b',0,0);

    end
end
axis equal
%shading flat
lighthandle=camlight('headlight');

end

