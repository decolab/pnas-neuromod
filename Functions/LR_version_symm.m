function [new]=LR_version_symm(TC)
% returns a symmetrical LR version of AAL 90x90 matrix

    odd=[1:2:90];
    even=sort([2:2:90],'descend');     
    new(1:45,:)=TC(odd,:);
    new(46:90,:)=TC(even,:);
    
 end