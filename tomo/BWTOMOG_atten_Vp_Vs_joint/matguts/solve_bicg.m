function [ dm ] = solve_bicg(F,f,tol,maxit)
%[ dm,bicgfn ] = solve_bicg(F,f,tol,maxit)
% Solve inverse problem using biconjugate gradient method

%% extract model vector

dm = bicg(@bicgfn,F'*f,tol,maxit );

end

function y = bicgfn(v,transp_opt)
global F
y = F'*(F*v);
end