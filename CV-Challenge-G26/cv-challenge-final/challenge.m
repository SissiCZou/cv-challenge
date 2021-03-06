%% Computer Vision Challenge
clear all;
close all;
% Groupnumber:
time1 = cputime;

group_number = 26;

%Groupmembers:
members = {'Verena Feierle','Ghazal Safian','Sarah Minwegen','Chenxi Zou','Weiyi Zhang'};

%Email-Adress(from Moodle):
mail = {' ga96soz@tum.de', 'ga96fob@tum.de', 'ga@tum.de', 'ga78woy@tum.de', 'ga38zaf@tum.de'};
%% load images
%TODO change RGB2GRAY to my function
%%%% load imageset1 or imageset2 
imageset = 1;  %for imageset2: imageset = 2
time1 = cputime;

imageset = input('Please give 1 for the imageset1, 2 for imageset2 (Default is 1):');
if imageset == 1
     Im1=imread('L1.JPG');GrayI1 = rgb2gray(Im1);
     Im2=imread('R1.JPG');GrayI2 = rgb2gray(Im2);
elseif imageset == 2
     Im1=imread('L2.JPG');GrayI1 = rgb2gray(Im1);
     Im2=imread('R2.JPG');GrayI2 = rgb2gray(Im2);   
end

%%%load corresponding points
if imageset == 1
        load('Data/corr1_1.mat')
        corr1 = corr1_1;
        load('Data/corr2_1.mat')
        corr2 = corr2_1;
elseif imageset == 2
        load('Data/corr1_2.mat')
        corr1 = corr1_1;
        load('Data/corr2_2.mat')
        corr2 = corr2_1;
end
%start executtion timer --> tic
%stop execution timer --> toc
%% get correspondence points from I1, and I2, using SIFT algorithms,very long time
% % computation time is about several minutes, better store the corr
% addpath('SIFT')
% [corr1,corr2] = get_correspondence(GrayI1,GrayI2);
% fprintf('************************************************** \n');
% % fprintf('The computation time using SIFT algorithm is %.2fs \n',(time2-time1));%670s
% % number of matching points
% fprintf('************************************************** \n');
% fprintf('The number of matched correspondence points using SIFT algorithm is %.2fs \n',no_matches);%705s
% rmpath('SIFT');
%% get homogene correspondence points in the center
% operations for note (1) above
no_matches = size(corr1,1);
x1 = [corr1 ones(size(corr1,1), 1)]';
x2 = [corr2 ones(size(corr2,1), 1)]';
% operations for note (2) above
siz = size(GrayI1);
origin = [siz(2); siz(1)]/2;
axis_x = -origin(1) : (origin(1)-1);
axis_y = (origin(2)-1) : -1 : -origin(2);
% T is the 3-by-3 transformation matrix required for operation (2) above
T = [1 0 -origin(1); 0 -1 origin(2); 0 0 1];
x1 = T*x1;
x2 = T*x2;
%% get Fundementalmatrix and epipoles
opt = lmeds_options('func', 'fundmatrix_nonlin', 'prop_outliers', 0.2, 'inlier_noise_level', 1);
[F,inl,outl,errs,avgerr] = lmeds([x1;x2], opt);

% [F,errs] = fundmatrix_ls([x1; x2], [], []);
%% Rotation and translation using E
K = 1.0e+03 * [2.4994         0    1.5018;
        0    2.5127    0.9642;
        0         0    0.0010];
E = K'*F*K;
[T1,R1,T2,R2]=TR_aus_E(E);
[Tans,Rota, lambda, P1] = rekonstruktion(T1,T2,R1,R2, [corr1 corr2]', K);
%% get rectified images and hormography matrics
[rec_im1, rec_im2, box, H1, H2] = rectify_images(Im1, Im2, x1, x2, F);% 14.9s
%% points corresponding to x1 and x2 in the new images are H1*x1 and H2*x2.
minx = box(1); miny = box(2);
maxx = box(3); maxy = box(4);
newCorrx1 = pflat(H1*x1);%new correspondence points 
newCorrx2 = pflat(H2*x2);

% %%
% % plot the input images
% figure;
% subplot(1,2,1)
% imagesc(axis_x, axis_y, Im1), axis xy, axis on, hold on
% 
% plot(x1(1,:), x1(2,:), 'g*')
% text(x1(1,:), x1(2,:), num2str( (1:no_matches)' ));
% title('First original image');
% % if imagetype == 'g', colormap gray; end
% subplot(1,2,2)
% imagesc(axis_x, axis_y, Im2), axis xy, axis on, hold on
% plot(x2(1,:), x2(2,:), 'g*')
% text(x2(1,:), x2(2,:), num2str( (1:no_matches)' ));
% title('Second original image');
% % if imagetype == 'g', colormap gray; end

% plot the outputs of rectified images
figure(1)
subplot(1,2,1)
imagesc(minx:maxx, maxy:-1:miny, rec_im1), axis xy, axis on, hold on
line([minx; maxx]*ones(1,no_matches), [newCorrx1(2,:); newCorrx1(2,:)]);
plot(newCorrx1(1,:), newCorrx1(2,:), 'g*')
%axis equal
title('First rectified image');
% if imagetype == 'g', colormap gray; end

subplot(1,2,2)
imagesc(minx:maxx, maxy:-1:miny, rec_im2), axis xy, axis on, hold on
line([minx; maxx]*ones(1,no_matches), [newCorrx2(2,:); newCorrx2(2,:)]);
plot(newCorrx2(1,:), newCorrx2(2,:), 'g*')
%axis equal
title('Second rectified image');
% if imagetype == 'g', colormap gray; end

%% get the disparity map using SAD algorithms
% time1 = cputime;

%resize the image to speed up zhe computation time
[spdmap_lr, ~, ~, ~] = stereomatch(rgb2gray(rec_im1), rgb2gray(rec_im2), 7, 254, 0);
[spdmap_rl, ~, ~, ~] = stereomatch(rgb2gray(rec_im2), rgb2gray(rec_im1), 7, 254, 0);
% resample the images to the half size,
% [spdmap_lr, ~, ~, ~] = stereomatch(rgb2gray(imresize(rec_im1,0.5)), rgb2gray(imresize(rec_im2,0.5)), 7, 128, 0);
% [spdmap_rl, ~, ~, ~] = stereomatch(rgb2gray(imresize(rec_im2,0.5)), rgb2gray(imresize(rec_im1,0.5)), 7, 128, 0);
% time2 = cputime;
% fprintf('The computation time for disparity map using stereo matching is %.2fs \n',(time2-time1));
figure(2);
subplot(211);imagesc(spdmap_lr);colorbar;title('Disparity Map LR');
subplot(212);imagesc(spdmap_rl);colorbar;title('Disparity Map RL');



%% get the synthesized image from the 3d coordinaten of objects for every pixels
%[out_view dmap rmap] = genIntView(interp, im1, im2, d1, d2);
% tic
%p = input('give p between 0 and 1:\n');
[smooth_view_p] = free_viewpoint(rec_im1,rec_im2,spdmap_lr,spdmap_rl,0.7);
% [smooth_view_p] = free_viewpoint(imresize(rec_im1,0.5),imresize(rec_im2,0.5),spdmap_lr,spdmap_rl,p);
% smooth_view_p = imresize(smooth_view_p,2);
% smooth_view_p = smooth_view_p(101:2100, 326:3325,:);

% [smooth_view_20] = free_viewpoint(rec_im1,rec_im2,spdmap_lr,spdmap_rl,0.2);
% [smooth_view_45] = free_viewpoint(rec_im1,rec_im2,spdmap_lr,spdmap_rl,0.45);
% [smooth_view_70] = free_viewpoint(rec_im1,rec_im2,spdmap_lr,spdmap_rl,0.70);
% [smooth_view_100] = free_viewpoint(rec_im1,rec_im2,spdmap_lr,spdmap_rl,1);
% toc 
%%
time2 = cputime;
elapsed_time = time2-time1;
fprintf('The computation time is %.2fs\n',elapsed_time);

%% plotting the endresult
% figure(3)
% 
% imshow(smooth_view_20);
% title('synthesized view, p = 0.2');
% 
% figure(4)
% 
% imshow(smooth_view_45);
% title('synthesized view, p = 0.45');
% 
% figure(5)
% imshow(smooth_view_70);
% title('synthesized view, p = 0.70');
% 
% figure(6)
% imshow(smooth_view_100);
% title('synthesized view, p = 1.00');