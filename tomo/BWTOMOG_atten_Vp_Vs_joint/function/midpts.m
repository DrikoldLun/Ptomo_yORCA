function [ B ] = midpts( A )
% [ out_vect ] = midpts( in_vect )
%   This function takes a vector of length N and returns a vector of length
%   N-1 consisting of the linear-interpolated midpoints of the N values in
%   the input vector
%   if A is higher dimension with size MxNx..., then the returned matrix 
%   will be linear-interpolated to midpoints in all dimensions with size
%   (M-1)x(N-1)x...

if isvector(A)
    B = 0.5*(A(1:end-1)+A(2:end));
elseif ismatrix(A)
    B = 0.25*(A(1:end-1,1:end-1) + A(1:end-1,2:end) + A(2:end,1:end-1) + A(2:end,2:end));
elseif ndims(A)==3
    B = (A(1:end-1,1:end-1,1:end-1) +...
         A(2:end  ,1:end-1,1:end-1) + A(1:end-1,2:end  ,1:end-1) + A(1:end-1,1:end-1,2:end  ) +...
         A(1:end-1,2:end  ,2:end  ) + A(2:end  ,1:end-1,2:end  ) + A(2:end  ,2:end  ,1:end-1) +...
         A(2:end  ,2:end  ,2:end  ))/8;
end



end
