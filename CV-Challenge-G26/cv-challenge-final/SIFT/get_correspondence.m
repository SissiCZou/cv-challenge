function [corr1,corr2] = get_correspondence(Im1,Im2)
% edited from swSift.m
% Load images
% Im1, Im2 can be both gray image or color image
[des1,loc1] = getFeatures(Im1);
[des2,loc2] = getFeatures(Im2);
matched = match(des1,des2);
nn = find(matched>0);
% drawFeatures(img1,loc1);
% drawFeatures(img2,loc2);
% drawMatched(matched,img1,img2,loc1,loc2);
corr1 = [loc1(nn,2) loc1(nn,1)];
corr2 = [loc2(matched(nn),2) loc2(matched(nn),1)];
% show matches correspondence points

end
