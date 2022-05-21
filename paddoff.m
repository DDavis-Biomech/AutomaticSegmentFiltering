function [b] = paddoff( a, nadd)

%  removes padding from a matrix, which has been Butterworth filtered
%  initial padding is required to circumvent endpoint problems
%
%  John H. Challis, The Penn. State University
%  April 26, 1997
%
%  calling - [b] = paddoff( a, nadd)
%  -------
%
%  input
%  -----
%  a    - matrix to have padding removed
%  output
%  ------
%  b    - matrix with padding removed
%  working
%  -------
%  nadd - how many padded values were added
%

%  find size of input matrix
%
n = size(a,1);


%  find start and finish
%
start=nadd+1;
fini=n-nadd;

%  switch data
%
b = a(start:fini,:);

%%% the end %%%
