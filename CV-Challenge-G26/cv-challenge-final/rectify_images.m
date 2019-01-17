function [newim1, newim2, b, H1, H2, newe] = rectify_images(im1, im2, x1, x2, F12)

%RECTIFY_IMAGES rectifies two images to achieve correspondences of scanlines.
%
%   [newim1, newim2] = rectify_images(im1, im2, x1, x2, F12)
%   [newim1, newim2, b, H1, H2] = rectify_images(im1, im2, x1, x2, F12)
%   [newim1, newim2, b, H1, H2, newe] = rectify_images(im1, im2, x1, x2, F12)
%   rectifies two images, im1 and im2, using the fundamental matrix, F12, which
%   must satisfy the epipolar equation:
%                          (x1)
%      (x2, y2, 1) * F12 * (y1) = 0
%                          ( 1)
% 
%   The input arguments are:
%   - im1 and im2 should both be an m-by-n array of doubles (or uint8) for some
%     values m and n
%   - x1 and x2 should both be 3-by-n matrix, where each column of the matrix
%     is an image point in homogeneous coordinates.  Corresponding columns in
%     the two matrices contain corresponding points in the two images.
%   - F12 must be a 3-by-3 rank-2 matrix.  Fundamental matrices can be computed
%     using one of Torr's routines (available for download on his Microsoft home
%     page) or Zhang's home page.
%
%     ***Note: the fundamental matrix is assumed to be computed with the following
%              image coordinate systems in both images being adopted: the origins
%              of the image coordinate systems are at the centres of the images;
%              the x-axes point to the right, y-axes up.
%
%   The output arguments are:
%   - the two new rectified images, newim1 and newim2.
%   - (optional) the bounding box, b, of the form [minx,miny,maxx,maxy] which bound
%     the new images newim1 and newim2.
%   - (optional) H1 and H2 are the computed rectification transformations.
%   - (optional) newe is the new epipole in the second image after the rectification.
%     On return, newe is always set [1;0;0] if the horizontal scanlines correspond
%     or [0;1;0] if the vertical scanlines correspond.
%
%   The implementation here is based on that described in the paper:
%   Richard I. Hartley, "Theory and Practice of Projective Rectification"
%   International Journal of Computer Vision, vol 35, no 2, pages 115-127, 1999.
%
%Created July 2002
%
%Copyright Du Huynh
%The University of Western Australia
%School of Computer Science and Software Engineering

% the epipole e1 is the projection of optical centre C2 onto image 1
% the epipole e2 is the projection of optical centre C1 onto image 2
% compute the two epipoles
e1 = null(F12); 
e2 = null(F12');

%[F, e1, e2] = fundmatrix(x1_new, x2_new);
%e1 = [-8039.69129604376;-630.136257971496;1];%% our e1_fund aus fumdenmentalmatrix.m
%e2 = [-1189.48183651498;715.145969174594;1];
% Check that both epipoles are outside the image boundary (a condition for
% Hartley's method)
pe1 = pflat(e1);  pe2 = pflat(e2);
nx = size(im1,2)/2;  ny = size(im1,1)/2;
if (in_range(pe1(1), -nx, nx) & in_range(pe1(2), -ny, ny))| ...
      (in_range(pe2(1), -nx, nx) & in_range(pe2(2), -ny, ny))
   error('rectify_images: rectification not possible if the epipole(s) is(are) inside the image(s)');
end

% compute the two 3-by-4 projection matrices
P1 = [eye(3) zeros(3,1)];
% note that for projective reconstruction we can choose the 3-by-4
% projection matrix P2 to have the form [ skew(e2)*F12+e2*alpha'  beta*e2 ], where
% alpha is an arbitrary 3-column vector and beta is an arbitrary scalar.  Here, for the rectification
% to work, the 3-by-3 submatrix of P2 (ie. the 1st three columns of P2) must be
% non-singular.  So, we adjust the vector alpha to make this submatrix non-singular.
% Since beta is also arbitrary, we can set beta to 1.
P2 = [skew(e2)*F12+e2*rand(1,3)  e2]; 

if nargout == 2
   [H1,H2] = rectification_transf(P1, P2, e1, e2, x1, x2);
elseif nargout == 5
   [H1,H2,newe,G2,R2] = rectification_transf(P1, P2, e1, e2, x1, x2);
end
   

nx = size(im1,2)/2;  ny = size(im1,1)/2;
% look for the smallest image size that encloses all the mapped corners
corners1 = pflat(H1*[-nx -ny 1; nx -ny 1; -nx ny 1; nx ny 1]');

nx = size(im2,2)/2;  ny = size(im2,1)/2;
corners2 = pflat(H2*[-nx -ny 1; nx -ny 1; -nx ny 1; nx ny 1]');

corners = [corners1 corners2];
minx = floor(min(corners(1,:))-1);
miny = floor(min(corners(2,:))-1);
maxx = ceil(max(corners(1,:))+1);
maxy = ceil(max(corners(2,:))+1);

b = [minx,miny,maxx,maxy];
newim1 = rectify_im(im1, H1, b);
newim2 = rectify_im(im2, H2, b);

end


function newim = rectify_im(im, H, box)

%RECTIFY_IM applies the rectification transformation to a given image.
%
%   newim = rectify_im(im, H, box) applies the rectification transformation
%   H (a 3-by-3 matrix) to the given image, im.  The bounding box
%   which determines the size of the output image newim, should be of the
%   format: [minx,miny,maxx,maxy], where minx and miny can be negative
%   numbers as the origin of the image coordinate system is set to the
%   centre of the image when the transformation H is computed.

minx = box(1); miny = box(2); maxx = box(3); maxy = box(4);
% the two matrices xx and yy returned by meshgrid are of the same dimension
[xx,yy] = meshgrid(minx:maxx, maxy:-1:miny);

% dimensions of the new (rectified) image
new_nrows = size(xx,1);  new_ncols = size(xx,2);

x = reshape(xx, 1, new_nrows*new_ncols); clear('xx');
y = reshape(yy, 1, new_nrows*new_ncols); clear('yy');

invH = inv(H);
len = length(x);
% We will encounter the "out of memory" problem if the image is large.
% So, to get around the problem, we do the following operation in
% several steps
mm = 50000;
idx=1:mm;
while (1)
   if idx(1) > len
      break;
   elseif idx(end) > len
      idx = idx(1:len-idx(1)+1);
   end
   newxy(:,idx) = invH*([x(idx); y(idx); ones(1,length(idx))]);
   newxy(:,idx) = pflat(newxy(:,idx));
   idx = idx + mm;
end
clear('idx');

% convert the x-y image coordinate system (origin at the image centre) to
% row-column coordinate system (origin at the top-left corner) before
% calling bilinear interpolation
nrows = size(im,1);  ncols = size(im,2);
newrc = [0 -1 nrows/2; 1 0 ncols/2; 0 0 1]*newxy;
clear('newxy');

% can't interpolate those points that fall outside the image im.  So discard them.
idx = find(newrc(1,:) >= 1 & newrc(1,:) <= nrows & ...
   newrc(2,:) >= 1 & newrc(2,:) <= ncols);
newrc = newrc(1:2,idx);
x = x(idx);
y = y(idx);

val = [];
len = size(newrc,2);
idx = 1:mm;
while (1)
   if idx(1) > len
      break;
   elseif idx(end) > len
      idx = idx(1:len-idx(1)+1);
   end
   val(idx,:) = bilinear_interpolate(im, newrc(:,idx));
   idx = idx + mm;
end

% compose the new image, newim, and the 2 arrays, x and y,
% into 1D array (note that x and y are still defined relative
% to the image origin at the centre of the image buffer.
% newim = zeros(new_ncols*new_nrows, size(im,3));
% covert into row-column coordinate system also
newim = zeros(new_nrows*new_ncols, 1, size(im,3));
rc = (x-minx)*new_nrows + (maxy-y) + 1;
if ~isempty(val)
   newim(rc,:) = val;
end

newim = reshape(uint8(newim), new_nrows, new_ncols, size(im,3));

end


function [ang,newe,R] = rectification_angle(e)

%RECTIFICATION_ANGLE computers the rectificaton angle, ang, required
%   to bring the epipole e to infinity.
%
%   [ang,newe,R] = rectificationAngle(e) also returns the following
%   output arguments:
%   - the new epipole, newe, which would either be [1;0;0] or [0;1;0]
%     depending on the position of the epipole e.
%   - the image plane (2D) rotation matrix, R, required for rectification_H.m
%
%Created July 2002.
%
%Copyright Du Huynh
%The University of Western Australia
%School of Computer Science and Software Engineering

e = e / norm(e,2);
% if e(3) is very small then the epipole must be at infinity already
if abs(e(3)) < 1E-20
   ang = 0;
   if abs(e(1)) > abs(e(2))
      newe = [1;0;0];
   else
      newe = [0;1;0];
   end
else
   theta = atan2(e(2)/e(3), e(1)/e(3))*ones(5,1);
   % below are 5 different angles to be considered.  If theta is
   % close to any of the first three angles then the rectification
   % should align the horizontal scanlines (new epipole should be
   % [1;0;0].  If theta is close to the last two angles then
   % the epipolar lines should be vertical (new epipole should be
   % [0;1;0].
   angles = [0; -pi; pi; -pi/2; pi/2];
   diff = angles - theta;
   [minval,minidx] = min(abs(diff));
   if minidx <= 3
      newe = [1;0;0];
   else
      newe = [0;1;0];
   end
   
   % compute the rotation angle require to bring the epipole e
   % to the x- or y-axis of the image
   ang = diff(minidx);
end

if nargout > 2
   R = [cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1];
end

end


function [H,newe,G,R] = rectification_H(e, orgin)

%RECTIFICATION_H computes the rectification transformation for the image
%   using its epipole e.
%
%   [H,newe] = rectification_H(e)
%   [H,newe,G,R] = rectification_H(e)
%   [H,newe,G,R] = rectification_H(e,origin) computes the rectification
%   transformation, H, for the image which has the epipole given by
%   the 3-by-1 vector, e.  The output argument, H, maps the epipole
%   e to a point at infinity -- typically, [1; 0; 0] for corresponding
%   horizontal scanlines or [0; 1; 0] for corresponding vertical scanlines.
%
%   The two optional output arguments, G and R, satisfy the equation
%      H = G*R*T,
%   with
%   - R being the 3-by-3 2D rotation matrix that brings the epipole e to
%     align with the x- or y-axis;
%   - G being the transformation that brings the epipole e, after rotation
%     R is applied, to infinity.
%   - T is the transformation that would transfer the origin of the image
%     coordinate system to the last input argument, origin.
%
%   The implementation here follows that described in "Theory and Practice
%   of Projective Rectification" by R. I. Hartley, IJCV'99.
if nargin == 1
   % assume that the origin of the image coordinate system has been
   % set at the centre of the image buffer
   T = eye(3);
else
   % assume that the origin was at the top-left corner.  We now move
   % it to the centre of the image, whose size is given by imsize.
   T = [1 0 -origin(1); 0 1 -origin(2); 0 0 1];
end

[ang,newe,R] = rectification_angle(e);
Re = R*e;
if sum(abs(newe-[1;0;0])) == 0
   % bring epipole e to parallel to [1;0;0]
   G = [1 0 0; 0 1 0; -Re(3)/Re(1) 0 1];
else
   % bring epipole e to parallel to [0;1;0]
   G = [1 0 0; 0 1 0; 0 -Re(3)/Re(2) 1];
end
H = G*R*T;

end


function [H1,H2,newe,G2,R2] = rectification_transf(P1, P2, e1, e2, x1, x2)

%RECTIFICATION_TRANSF computes the rectification transformation matrices.
%
%   [H1,H2] = rectification_transf(P1, P2, e1, e2, x1, x2)
%   [H1,H2,newe,G2,R2] = rectification_transf(P1, P2, e1, e2, x1, x2)
%
%   Input arguments are:
%   - P1 and P2 should both be 3-by-4 projection matrices, with P1 = [I 0].
%   - e1 and e2 should both be 3-by-1 vectors containing the two epipoles, with e1
%     being the epipole in image 1 (the projection of the optical centre of the
%     second camera onto the first image), and e2 being the epipole in image 2
%     (projection of the optical centre of the first camera onto the second image).
%   - x1 and x2 should both be 3-by-n matrices with each column of the matrix being
%     an image point in homogeneous coordinates and n being is the number of matching
%     points.
%
%   Output arguments are:
%   - H1 and H2 would both be 3-by-3 rectification matrices.
%   - G2 and R2 are all optional output arguments.  If specified, they would
%     all be 3-by-3 matrices, satisfying the condition H2 = G2*R2.
%
%   SEE ALSO rectification_H.m
% check input argument P1 -- it must be a [I 0] matrix
if max(max(abs(P1-[eye(3) zeros(3,1)]))) ~= 0
   error('rectification_transf: matrix P1 must be of the form [I 0]');
end

% compute rectification matrix H2 for the second image
if nargout >= 4 & nargout <= 5
   [H2,newe,G2,R2] = rectification_H(e2);
else
   [H2,newe] = rectification_H(e2);
end

H0 = H2*P2(:,1:3);
x2hat = pflat(H2*x2);
x1hat = pflat(H0*x1);
B = [
   sum( (ones(3,1)*x1hat(1,:)) .* x1hat, 2 )';
   sum( (ones(3,1)*x1hat(2,:)) .* x1hat, 2 )';
   sum( (ones(3,1)*x1hat(3,:)) .* x1hat, 2 )'
];

if sum(abs(newe-[1;0;0])) == 0
   b = [
      sum( x1hat(1,:) .* (x2hat(1,:)-x1hat(1,:)) );
      sum( x1hat(2,:) .* (x2hat(1,:)-x1hat(1,:)) );
      sum( x1hat(3,:) .* (x2hat(1,:)-x1hat(1,:)) )
];
else
   b = [
      sum( x1hat(1,:) .* (x2hat(2,:)-x1hat(2,:)) );
      sum( x1hat(2,:) .* (x2hat(2,:)-x1hat(2,:)) );
      sum( x1hat(3,:) .* (x2hat(2,:)-x1hat(2,:)) )
   ];
end

abc = B \ b;

if sum(abs(newe-[1;0;0])) == 0
   A = [[1 0 0]+abc'; 0 1 0; 0 0 1];
else
   A = [1 0 0; [0 1 0]+abc'; 0 0 1];
end
H1 = A*H0;

end

function bool = in_range(v, minv, maxv)
%IN_RANGE checks if a given value is within a specified range.
%bool = in_range(v, minv, maxv) sets bool to 1 if minv <= v <= maxv,
%0 otherwise.
if v >= minv & v <= maxv
   bool = 1;
else
   bool = 0;
end
end

function M = skew(t)

M = [0 t(3) -t(2);
   -t(3) 0 t(1);
   t(2) -t(1) 0];

end


function val = bilinear_interpolate(im, rc)

%BILINEAR_INTERPOLATE performs bilinear interpolation.
%
%   val = bilinear_interpolate(im, rc) returns the bilinearly interpolated
%   pixel values at the given list of positions, rc.
%
%   Input arguments:
%   - im should be a m-by-n-by-3 (colour) image or a m-by-n (grey scale) image.
%   - rc should be a 2-by-k matrix where k is the number of positions
%     whose values are to be computed via bilinear interpolation.  The first
%     row of the matrix should contain the row-coordinates and the second row the
%     column-coordinates of these positions.  It is important that all the
%     (row,column)-positions stored in the matrix, rc, are within the boundary
%     of the image im.
%
%   Output argument:
%   - val is the output k-by-3 (if image im has 3 colour bands) or k-by-1
%     (if image im is a grey scale image) matrix containing the interpolated pixel
%     values.

nrows = size(im,1);  ncols = size(im,2);
% number of bands (3 for coloured images; 1 for gray scale images)
nobands = size(im,3);

% check that all the entries in matrix rc are within the boundary of image im.
if sum(rc(1,:) < 1 | rc(1,:) > nrows | rc(2,:) < 1 | rc(2,:) > ncols)
   error('bilinear_interpolate: elements of the rc matrix must be within the image boundary');
end


% The four corner points used in the bilinear interpolation:
%  c4      c3
%  +--------+
%  |        |
%  | o      |     (a point o which is stored in a column vector of rc
%  |        |     and the two corner points used for its bilinear interpoation)
%  +--------+
%  c1      c2
% The row-column coordinate system is used.  All the corner points c1, c2, c3,
% and c4 are two vectors, whose 1st components contains the row coordinates
% and 2nd components contains the column coordinates.

% note that we should have given a small margin for variable rc
% so that the four corner pixels surrounding rc are within the
% image boundary of im.  This condition should be enforced in the
% caller of this function.
c4 = (floor(rc));
c2 = (ceil(rc));
c3 = ([c4(1,:); c2(2,:)]);
c1 = ([c2(1,:); c4(2,:)]);

% d(diffRC_idx) = (rc(diffRC_idx) - c1) ./ (c4 - c1 + eps);
%
% the interpolation procedure above fails for those points in
% rc whose x- or y- component is a whole number (in which case,
% the respective components of these points in the c1 and c4 matrices
% would be the same.  the formula for d below would cause a division
% by zero problem.
sameC_idx = find(c2(2,:) == c4(2,:));
sameR_idx = find(c2(1,:) == c4(1,:));
diffRC_idx = find(c2(1,:) ~= c4(1,:) & c2(2,:) ~= c4(2,:));
% now the formula for d can be safely applied...
d(:,diffRC_idx) = (rc(:,diffRC_idx) - c4(:,diffRC_idx)) ./ ...
   (c2(:,diffRC_idx) - c4(:,diffRC_idx) + eps);
d(2,sameC_idx) = 0;
d(1,sameC_idx) = (rc(1,sameC_idx) - c4(1,sameC_idx)) ./ ...
   (c2(1,sameC_idx) - c4(1,sameC_idx) + eps);
d(1,sameR_idx) = 0;
d(2,sameR_idx) = (rc(2,sameR_idx) - c4(2,sameR_idx)) ./ ...
   (c2(2,sameR_idx) - c4(2,sameR_idx) + eps);

% convert c1, c2, c3, c4 into 1D array for fast retrieval of image
% intensity from im
c1 = (c1(2,:)-1)*nrows + c1(1,:);
c2 = (c2(2,:)-1)*nrows + c2(1,:);
c3 = (c3(2,:)-1)*nrows + c3(1,:);
c4 = (c4(2,:)-1)*nrows + c4(1,:);

im = reshape(im, nrows*ncols, nobands);
c1val = im(c1,:);
c2val = im(c2,:);
c3val = im(c3,:);
c4val = im(c4,:);

for i=1:nobands
   val(:,i) = (1-d(1,:)').*(1-d(2,:)').*double(c4val(:,i)) + ...
      d(1,:)'.*(1-d(2,:)').*double(c1val(:,i)) + ...
      (1-d(1,:)').*d(2,:)'.*double(c3val(:,i)) + ...
      d(1,:)'.*d(2,:)'.*double(c2val(:,i));
   val(:,i) = uint8(val(:,i));
end

end