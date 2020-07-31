function [totaldistance] = cumulativeMotion(motion)
% start with scaled motion path integrate to get total distance travelled
totaldistance=0;
for i=1:length(motion)-1
    totaldistance=totaldistance+distance(motion(i,:)',motion(i+1,:)');
end
end

