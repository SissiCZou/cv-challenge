function [smooth_view] = free_viewpoint(rec_im1,rec_im2,spdmap_lr,spdmap_rl,p)
%%  read color images and disparity maps for left and right view and convert them to double values.
%   im0: color image of left view
%   im1: color image of right view
%   disp0: disparity image of left view
%   disp1: disparity image of right view
%   theta: parameter for new view location, theta can be any number from 0 to 1
im0 = im2double(rec_im1);
im1 = im2double(rec_im2);
disp0 = double(spdmap_lr);
disp1 = double(spdmap_rl);
%% synthesizing
[im0,im1] = removeGhostContour(im0,im1,disp0,disp1); % remove ghost contour artifacts
[int_view,int_dist] = initialSynthesize(im0,im1,disp0,disp1,p); % synthesis new view and deptg map
%[final_view,final_dist] = depthHoleFill(int_view,int_dist); % fill holese in depth map;
[fl_view,fl_dist] = colorHoleFill(int_view,int_dist);   % fill holes in color image
[smooth_view] = edgeSmooth(fl_view,fl_dist);    % smooth edge for a more realistic look
smooth_view = smooth_view(101:2100,326:3325, :);
end