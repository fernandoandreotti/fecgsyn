% PRCORR2  - Compute correlation coefficient over all dimentions
%
%    This is a 10 time faster implementation of MATLABs corr2.
%    Implemented as a mex-file generated with the ITERATOR
%    tool that can be downloaded from www.mathworks.com
%    Re compile the source with MEX -O prcorr2.c
%    
%
%    C = PRCORR2(A,B) computes the correlation between A and B.
%    A and B are arrays witn any number of dims but with same 
%    number of elements
% 
%    Class Support
%    -------------
%    A and B can be any numeric type but not complex value.
%    C is a scalar double.
% 
%    See also CORR2, CORRCOEF, STD2, ITERATOR, MEX
%
%
%    (C) 2003 Peter Rydesäter,  http://www.rydesater.com
