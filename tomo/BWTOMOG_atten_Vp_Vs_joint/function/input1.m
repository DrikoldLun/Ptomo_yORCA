function [ out ] = input1(prompt,numopt)
% [ out ] = input1(prompt,numopt)
% function to get one character of user input without pressing enter
%
% 'prompt' = prompt that appears at the command line
% 'stropt' = option specifying  input as number (1) or string (0)

if nargin<2
    numopt=0;
end

fprintf(prompt)

w = waitforbuttonpress;
if w 
    out = get(gcf,'CurrentCharacter');
end

if numopt==1
    out = str2double(out);
end

end

