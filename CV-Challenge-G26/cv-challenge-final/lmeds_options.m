function options = lmeds_options(varargin)

%LMEDS_OPTIONS sets options for the LMEDS function.
%
%   options = lmeds_options('func', funcName, ...)
%   options = lmeds_options('func', funcName,'prob', probValue, ...)
%   options = lmeds_options('func', funcName,'prop_outliers', proportionOfOutlier, ...)
%   options = lmeds_options('func', funcName,'sample_size', sampleSizeValue, ...)
%   options = lmeds_options('func', funcName,'inlier_noise_level', inlierNoiseLevel, ...)
%
%   lmeds_options is to be used in conjunction with lmeds, which implements
%   the Least Median Squares framework.  This function should be used for
%   setting various parameters that are needed by lmeds.  The input arguments
%   must be in the form of attribute-value pairs.  Below are the recognized
%   attributes and their value types:
%
%   Attribute             Value type
%   'func'                string, to specify the name of the function (default is
%                              'fundmatrix_nonlin') that carries out the
%                              implementation for estimating the entity
%                              that is of interest.
%   'prob'                double (between 0 and 1) to specify the probability
%                              that at least one of the samples drawn is good.
%   'inlier_noise_level'  double.  This is the level of inlier noise, e.g.
%                              the standard deviation of Gaussian noise
%                              inherent in the data.
%   'prop_outliers'       double (between 0 and 1) that stores the
%                              proportion of outliers in the data
%   'dof'                 positive integer that stores the degrees of freedom.  This
%                              value will be used for computing the
%                              threshold value required for classifying
%                              inliers and outliers (using the chi-square test)
%   'sample_size'         positive integer to denote the size of each sample in the
%                              random drawing process.  e.g. if 'func' is
%                              'fundmatrix_nonlin' then 'sample_size' should be 7.
%
%   Not all of the above attributes need to be specify, as their default
%   values would be automatically computed if the function name is a
%   recognized one.  For instance, the values for the attributes 'sample_size'
%   and 'dof' can be automatically computed.  These default values can, of
%   course, be overwritten by the specified arguments.
%
%   The output argument, options, is a struct containing the fields:
%   func, prob, prop_outliers, inlier_noise_level, sample_size, dof,
%   no_samples, where the last field, no_samples, denotes the number of
%   samples to be drawn and is automatically computed.
%
%   If you are not sure what values to use for all the fields, it is
%   recommended that you call lmeds_options by providing only the function
%   name and let lmeds_options determine the default values for other
%   fields for you.
%   
%   Examples:
%   (a) use the default parameter setting (see below) for the function fitline_ls
%         options = lmeds_options('func', 'fitline_ls');
%   (b) set the proportion of outliers to 0.3 and use the default values
%       for other parameters for the function fundmatrix_nonlin:
%         options = lmeds_options('func', 'fundmatrix_nonlin', ...
%                                 'prop_outliers', 0.3);
%   (c) set the proportion of outliers to 45 percent and sample size to 5 for
%       fitplane_ls:
%         options = lmeds_options('func', 'fitplane_ls', ...
%                                 'prop_outliers', 0.45, 'sample_size', 5);
%       Note that in this example the sample size is in fact too large for line
%       fitting.
%
%   Below is the list of recognized function names (This list will
%   no doubt be expanded in the future) and the default values for two of the
%   listed attributes:
%                         sample_size   dof
%   fundmatrix_ls              8         2
%   fundmatrix_wls             8         2
%   fundmatrix_nonlin          7         2
%   fitline_ls                 2         1
%   fitplane_ls                3         1
%
%   The following default values are common to all the functions above:
%     prob = 0.99
%     prop_outliers = 0.4 
%
%   *** Note:
%   If the inlier noise level is unknown then by default if would be set to -1.
%   A different formula would then be used in function lmeds to classify
%   inliers and outliers.
%
%   SEE ALSO fundmatrix_ls, fundmatrix_wls, fundmatrix_nonlin,
%            fitline_ls, fitplane_ls, lmeds
%
%   Created September 2003.
%
%   Copyright Du Huynh
%   The University of Western Australia
%   School of Computer Science and Software Engineering

funcName = 'UNDEFINED';
% default values used by all functions
prob = 0.99;
prop_outliers = 0.4;
inlier_noise_level = -1;

if mod(nargin, 2)
   % number of input arguments must be even
   error('lmeds_options: number of input arguments must be even.');
else
   for i=1:2:length(varargin)
      switch varargin{i}
      case 'func'
         funcName = varargin{i+1};
      case 'prob'
         prob = varargin{i+1};
      case 'sample_size'
         sample_size = varargin{i+1};
      case 'prop_outliers'
         prop_outliers = varargin{i+1};
      case 'inlier_noise_level'
         inlier_noise_level = varargin{i+1};
      case 'dof'
         dof = varargin{i+1};
      otherwise
         error(sprintf('lmeds_options: unrecognised field %s.', varargin{i}));
      end
   end
end

switch funcName
case 'UNDEFINED'
   error('lmeds_options: function name not given.');
case {'fundmatrix_ls','fundmatrix_wls'}
   dof = 2;
   samp_size = 8;
case 'fundmatrix_nonlin'
   dof = 2;
   samp_size = 7;
case 'fitline_ls'
   dof = 1;
   samp_size = 2;
case 'fitplane_ls'
   dof = 1;
   samp_size = 3;
case 'trifocal_6pts'
   dof = 3;
   samp_size = 6;
otherwise
   if ~exist('dof')
      warning('lmeds_options: dof field set to -1');
      dof = -1;
   end
end

% use samp_size only if sample_size is not specified by the user.
% i.e. let the user overwrite the default sample size value set up
% above for each function.
if ~exist('sample_size')
   sample_size = samp_size;
end

options = struct(...
   'func', funcName, ...
   'prob', prob, ...
   'prop_outliers', prop_outliers, ...
   'sample_size', sample_size, ...
   'inlier_noise_level', inlier_noise_level, ...
   'dof', dof, ...
   'no_samples', lmeds_comp_no_samples(prob, prop_outliers, sample_size));   

return

% -----------------------------------------------------------

function no_samples = lmeds_comp_no_samples(prob, prop_outliers, sample_size)

no_samples = ceil( log(1-prob) / log(1-(1-prop_outliers)^sample_size) );

return
