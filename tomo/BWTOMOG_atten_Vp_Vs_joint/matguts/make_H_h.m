function [ H,h ] = make_H_h( F,f,nmodel,static )
% fixed estatic and sstatic
% static should be column vector
C = [zeros(length(static),nmodel),eye(length(static))];
H = [2*(F)'*F,C';C,zeros(size(C,1))];
h = [2*(F)'*f;static];
end