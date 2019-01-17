function [F,errs,avgerr] = fundmatrix_ls(xy1xy2, idx, ws, scale)

%FUNDMATRIX_LS computes the fundamental matrix using the linear method.
%
%   F = fundmatrix_ls(xy1xy2, idx, ws)
%   F = fundmatrix_ls(xy1xy2, idx, ws, scale)
%   [F,errs] = fundmatrix_ls(xy1xy2, idx, ws)
%   [F,errs] = fundmatrix_ls(xy1xy2, idx, ws, scale)
%   [F,errs,avgerr] = fundmatrix_ls(xy1xy2, idx, ws)
%   [F,errs,avgerr] = fundmatrix_ls(xy1xy2, idx, ws, scale)
%
%   estimates the fundamental matrix using the linear method and assuming that
%   all corresponding points are inliers.  The corresponding point coordinates
%   should be stored in the large 6-by-n matrix of the following format:
%      xy1xy2 = [x1_1  x1_2 ... x1_n;
%                y1_1  y1_2 ... y1_n;
%                z1_1  z1_2 ... z1_n;
%                x2_1  x2_2 ... x2_n;
%                y2_1  y2_2 ... y2_n;
%                z2_1  z2_2 ... z2_n ]
%   where [x1_i, y1_i, z1_i] is the i-th feature point (in homogeneous coordinates)
%         in the first image;
%         [x2_i, y2_i, z2_i] is the i-th feature point (in homogeneous coordinates)
%         in the second image;
%         each column of the matrix stores a pair of corresponding points in
%         homogeneous coordinates.
%
%   The fundmental matrix, F, is a rank-2 3-by-3 matrix satisfying the epipolar
%   constraint:
%                            (x1_i)
%     (x2_i y2_i z2_i) * F * (y1_i) = 0
%                            (z1_i)
%
%   Input arguments:
%   - xy1xy2 should be a 6-by-n matrix where n >= 8.  It is recommended that
%     the coordinates of points in this matrix be relative to an image coordinate
%     system whose origin is at the centre (or principal point, if known) of
%     the image.
%   - idx must be an array of integers in the range [1,n] where n=size(xy1xy2,2).
%     Alternatively, set idx to [] if all the points in xy1xy2 are intended
%     to be used for the estimation of F.
%   - ws should be an array of weighting factors, one for each pair of corresponding
%     points.  If a normal linear method is desired (ie. no weighting factors to
%     be applied) then ws can be set to [].
%   - scale (optional) is the scale factor to be used in least squares to keep
%     things well behaved.  If not specified, a scale factor would be automatically
%     computed (Referencee: see Hartley's ICCV'95 paper).  If specified, a sensible
%     scale factor is half of the average of the image dimensions, eg. if the images
%     are of dimensions 400-by-600 then scale should be set to 500/2 = 250.
%
%   Output arguments:
%   - F is the rank-2 3-by-3 fundamental matrix.
%   - errs is the list of reprojection errors (orthogonal distance errors of
%     the feature points from the epipolar lines) of all the feature point in xy1xy2.
%   - avgerr is the average squared reprojection error, i.e. the mean value of errs.
%
%   Examples:
%   (a)   [F,errs,avgerr] = fundmatrix_ls(xy1xy2, [], []);
%   (b)   [F,errs,avgerr] = fundmatrix_ls(xy1xy2, ones(size(xy1xy2,2), 1), []);
%   (c)   [F,errs,avgerr] = fundmatrix_ls(xy1xy2, [], 1:10);
%   In examples (a) and (b), all the corresponding points in xy1xy2 would be used
%   to estimate F and they all have the weighting factors of 1.  Both examples would
%   produce the same result.
%   In example (c), no weighting factors are specified and only the first 10 pairs of
%   corresponding points would be used for estimating F.  The third argument,
%   idx, thus allows the function to be used in the RANSAC/LMedS framework.   
%
%Created 1996.
%Revised September 2003.
%
%Copyright Du Huynh
%The University of Western Australia
%School of Computer Science and Software Engineering

% check dimensions of input argument
if isempty(idx)
   % the user doesn't care about the list idx, so will use all
   % data points
   idx = 1:size(xy1xy2,2);
elseif any(idx <= 0 | idx > size(xy1xy2,2))
   error('fundmatrix_ls: idx must be an array of integers in the range [1,n], where n=size(xy1xy2,2).');
end

xy = xy1xy2(:,idx);

if nargin < 4
   % compute scale factor and shifting
   xy1 = pflat(xy(1:3,:));
   xy2 = pflat(xy(4:6,:));
   % this scaling is not quite the same as that described in the Hartley's
   % iccv95 paper (and others' papers), but the difference to the estimated
   % fundamental matrix would be insignificant for other ways of scaling.
   % It is more important that the data points are well scattered in the
   % images and are not degenerate (e.g. the image points are projections
   % of coplanar scene points)
   scale = sqrt(2) / ( max(max([xy1;xy2]))-min(min([xy1;xy2])) );
   T1 = [scale 0 -scale*mean(xy1(1,:));
      0 scale -scale*mean(xy1(2,:));
      0 0 1];
   T2 = [scale 0 -scale*mean(xy2(1,:));
      0 scale -scale*mean(xy2(2,:));
      0 0 1];
elseif scale == 0
   T1 = eye(3);  T2 = eye(3);
else
   T1 = diag([1/scale, 1/scale, 1]);
   T2 = T1;
end

% transform the corresponding points by T1 and T2
xy1new = T1*xy1;
xy2new = T2*xy2;

A(:,1) = xy2new(1,:)' .* xy1new(1,:)';
A(:,2) = xy2new(1,:)' .* xy1new(2,:)';
A(:,3) = xy2new(1,:)' .* xy1new(3,:)';

A(:,4) = xy2new(2,:)' .* xy1new(1,:)';
A(:,5) = xy2new(2,:)' .* xy1new(2,:)';
A(:,6) = xy2new(2,:)' .* xy1new(3,:)';

A(:,7) = xy2new(3,:)' .* xy1new(1,:)';
A(:,8) = xy2new(3,:)' .* xy1new(2,:)';
A(:,9) = xy2new(3,:)' .* xy1new(3,:)';

if ~isempty(ws)
   A = diag(ws)*A;
end

% Note that f is a 9-vector
f = lsq(A'*A);
F = reshape(f,3,3)';
% impose the rank-2 constraint
F = imposeRank(F,2);
% Undo the scaling
F = T2'*F*T1;
% scale F so that ||f|| = 1
F = F / norm(F,2);

if nargout >= 2
   xy1 = xy1xy2(1:3,:);
   xy2 = xy1xy2(4:6,:);
   % points in the first image cast epipolar lines in the second image under
   % the mapping of the fundamental matrix
   line2 = F*xy1;
   % similarly...
   line1 = F'*xy2;
   
   % squared reprojection errors (see Hartley & Zisserman's book, Eq (10.10), p.271),
   % (Note: the errors are chi-square variables of 2 degrees of freedom)
   % Here, I keep the array of errors rather than summing them together as in Eq (10.10)
   errs = (sum(line2.*xy2,1)).^2 .* ...
      (1./(line2(1,:).^2 + line2(2,:).^2) + 1./(line1(1,:).^2 + line1(2,:).^2));
end

if nargout >= 3
   % this is the average value of the reprojection errors, i.e. Eq (10.9) divided
   % by the number of corresponding points (see Hartley & Zisserman's book)
   avgerr = mean(errs);
end

end 

% function imposeRank()
function M = imposeRank(M,r)
%IMPOSERANK imposes the rank constraint to a given matrix.
%   M = imposeRank(M,r)
%   imposes the rank constraint to the given matrix M and returns
%   a new matrix which has rank = r.
[n,m] = size(M);
nm = min([n,m]);
if r < nm
   [U,S,V] = svd(M,0);
   S(r+1:nm,r+1:nm) = zeros(nm-r,nm-r);
   M = U*S*V';
end
end

function [x,sval,sval_ratio] = lsq(A)

%LSQ returns the least squares solution of the given matrix.
%
%   x = lsq(A)
%   [x,sval] = lsq(A)
%   [x,sval,sval_ratio] - lsq(A)
%   returns the vector, x, that corresponds to the smallest singular value of
%   matrix A.  That is, the function returns the vector, x, that minimises
%           || A*x ||^2
%   Input argument:
%   - A the matrix of interest
%   Output arguments:
%   - vector x is the least squares solution
%   - sval (optional) is the smallest singular value
%   - sval_ratio is the ratio between the smallest singular value and the
%     next smaller singular value.  The magnitude of the values sval and
%     sval_ratio gives an indication of how accurate x is.
[U,S,V] = svd(A);
last = size(V,2);
x = V(:,last);
if nargout >= 2
   s = diag(S);
   if length(s) < last
      sval = 0;
   else
      sval = s(last);
   end
end
if nargout >= 3
   sval_ratio = sval / s(last-1);
end

end 
