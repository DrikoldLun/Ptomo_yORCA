function logic_if_any_match = regexpany(str,expression,varargin)
%     logic_if_any_match = regexp(str,expression,varargin)
% 
%  Function to find if the string "expression" appears anywhere in a number
%  of strings aranged in a cell array, returning a logical output the same
%  size as the array to be tested. 

if isempty(varargin)
    a = regexp(str,expression);
else
    a = regexp(str,expression,varargin);
end

logic_if_any_match = false(size(str));

for ii = 1:numel(a)
    logic_if_any_match(ii) = any(a{ii});
end

end