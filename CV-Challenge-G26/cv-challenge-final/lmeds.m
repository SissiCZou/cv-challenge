function [estx, inlier_idx, outlier_idx, errs, avgerr] = ...
   lmeds(data, options, condfunc, varargin)

%LMEDS is a framework for parameter estimation using the LMedS method.
%
%   [estx,...] = ...
%         lmeds(data,options,condfunc,varargin)
%   [estx,inlier_idx,outlier_idx,errs,avgerr] = ...
%         lmeds(data,options,condfunc,varargin)
%
%   lmeds is the LMedS framework for sampling data for parameter estimation.
%   The function takes in a variable number of input and produces a variable
%   number of arguments.
%
%   Input arguments:
%   - data is the observed data for the parameter estimation.  The rule that I
%     impose here is data must be a matrix in which each column must
%     contain one observed datum. For instance,
%     (1) for the estimation of the fundamental matrix, each column of the
%         argument, data, must contain a pair of corresponding points, as
%         a 6-vector (in homogeneous coordinates).
%     (2) for 2D line fitting, each column of data must contain one 2D point,
%         as a 3-vector (in homogeneous coordinates).
%   - options is the struct produced by lmeds_options.  It contains important
%     information about the function to be used for parameter estimation and the
%     various parameters required by the LMedS paradigm, etc.
%   - condfunc is a conditional function for checking that the randomly selected
%     sample is well-behaved (i.e. not degenerate).  You can set condfunc=[] if
%     no conditional function is available.
%   - varargin is a list of optional arguments that you want to pass to function
%     func.  See the examples given below for details.
%
%   Output arguments:
%   - estx stores the estimated parameter, in whatever form produced by your
%     supplied function, func.
%   - inlier_idx (optional): the indices of all the detected inliers.
%   - outlier_idx (optional): the indices of all the detected outliers.
%     The union of inlier_idx and outlier_idx gives 1:n, where n is the
%     number of columns in the argument, data.
%   - errs (optional argument) stores the errors of fitting for the estimated
%     parameter, estx.  If the input argument, data, contains n columns then
%     errs would be a 1-by-n array of error value, one for each observed datum.
%   - avgerr (optional argument) is simply the average value of errs.
%
%   Example 1:
%   options = lmeds_options('func', 'fundmatrix_ls');
%   [F,inlier_idx,outlier_idx] = lmeds(xy, options, []);
%
%   Example 2:
%   options = lmeds_options('func', 'fundmatrix_nonlin', 'prop_outliers', 0.4);
%   [F,inlier_idx,outlier_idx] = lmeds(xy, options, [], 0);
%
%   Example 3:
%   options = lmeds_options('func', 'fitline_ls', 'inlier_noise_level', 1);
%   [abc,inlier_idx,outlier_idx,errs,avgerr] = lmeds(x, options, [], 0);
%
%   In all the examples above, no condfunc is supplied (as shown by the third
%   argument being []).  In examples 1 and 2, the argument xy must be a
%   6-by-n matrix for some n >= 7 or n >= 8.  In example 3, x must be a 3-by-m matrix
%   for some m >= 2, and an additional argument, 0, may be supplied.
%   
%   SEE ALSO fundmatrix_ls, fundmatrix_nonlin, fundmatrix_wls,
%            fitline_ls, fitplane_ls, lmeds_options
%
%Created 2000.
%Last modified September 2003.
%
%Copyright Du Huynh
%The University of Western Australia
%School of Computer Science and Software Engineering

if nargin == 2
   condfunc = [];
end

% number of observations in the data
no_pts = size(data,2);
sample_size = options.sample_size;

med_err = [];
index = 0;
% it is more convenient to use a cell array for the variable x below.
% However, it seems to slow things down a bit...
if options.no_samples == 0, options.no_samples = 1; end;
for i = 1 : options.no_samples
   % pick random data points to form a sample.  Ensure that the data points
   % are well-behaved (as defined by the condfunc function)
   while 1
      idx = randperm(no_pts);
      idx = idx(1:sample_size);
      if isempty(condfunc)
         break;
      elseif feval(condfunc, data(:,idx), varargin{:})
         break;
      end
   end     
   
   % now call the function func to estimate the parameter x.  Note
   % that there may be more than one solution returned to soln
   [soln,error] = feval(options.func, data, idx, varargin{:});
   if iscell(soln)
      % there are more than one solution output by function func.
      % Inspect them all.
      for k=1:length(soln)
         index = index + 1;
         err{index} = error{k};
         med_err(index) = median(error{k});
         x{index} = soln{k};
      end
   else
      % there is only one solution returned from function func
      index = index + 1;
      err{index} = error;
      med_err(index) = median(error);
      x{index} = soln;
   end
end

% find the sample that gives the smallest median squares error
[lmed_err,best_index] = min(med_err);

% x{best_index} is the best solution!
estx = x{best_index};

% now we need to classify the data points into inliers and outliers,
% i.e. we need to compute the indices of the inliers and outliers.
if options.inlier_noise_level > 0
   % the user has specified a known inlier noise level (i.e. the standard
   % deviation of noise)
   threshold = options.inlier_noise_level^2 * chi2inv(options.prob, options.dof);
else
   sigmahat = 1.4826*( 1 + 5/(no_pts - options.sample_size) )*sqrt(lmed_err);
   threshold = 2.0*sigmahat;
end

inlier_idx = find(err{best_index} < threshold);
outlier_idx = find(err{best_index} >= threshold);

if nargout > 3
   errs = err{best_index};
end

if nargout > 4
   avgerr = mean(errs);
end

return
