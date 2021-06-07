function [Y,Mo,D,Jday,H,M,S] = iristime2greg( irisdate, iristime )
% function [Y,Mo,D,Jday,H,M,S] = iristime2greg( irisdate, iristime )

% break down event times dates into Y,Mo,D,JulD
        [Y,remain]=strtok(char(irisdate),'/');  Y=eval(Y);
        [Mo,remain]=strtok(remain,'/');         Mo=eval(Mo);
        [D,remain]=strtok(remain,'/');          D=eval(D);
        %annoying things to get d.o.y. in right format (3-digit string)
        Jday=(datenum(Y,Mo,D)-datenum(Y,0,0)); 
        if Jday >= 100 
            assignin('base','Jday',sprintf('%s',num2str(Jday)));
        elseif Jday >= 10 && Jday < 100
            assignin('base','Jday',sprintf('0%s',num2str(Jday)));
        elseif Jday < 10
            assignin('base','Jday',sprintf('00%s',num2str(Jday))); 
        end
% break down event times dates into H,M,S
        [H,remain]=strtok(char(iristime),':');  H=eval(H);
        [M,remain]=strtok(remain,':');          M=eval(M);
        [S,remain]=strtok(remain,':');          S=eval(S);
end

