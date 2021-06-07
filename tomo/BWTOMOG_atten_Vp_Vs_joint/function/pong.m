%% version of PONG in matlab

% game is played on 10x10 grid
function pong()
%% options
% game mode
opt1 = input('player mode [p] or automatic [a]?  ','s');
opt2 = input('difficulty level (scale from 1-10)?  ');
ex = 0;
d = 0;
dyp = 0;

% paddle dimensions
dpad = 0.3;
hpad = 1.7;

%% starting positions
% ball start
x_ball = 0;
y_ball = random('unif',-10+hpad,10-hpad,1);
u_ball = random('unif',1,3,1); % dx/dt
v_ball = random('unif',0,2,1); % dy/dt


% l_start
y_l = 2*(10-hpad)*(rand-0.5);

% r_start
y_r = 2*(10-hpad)*(rand-0.5);

% plot start
figure(10); clf; hold on
plot_sit(gca,x_ball,y_ball,y_l,y_r,dpad,hpad)

figure('KeyPressFcn',@movepaddle)

function movepaddle(splooge,event) % move paddle function
    switch event.Character
        case 30 % up arrow
            d = 1; % increase amplitudes
        case 31  % down arrow
            d = -1;
        case 29  % right
            d = 0;
        case 28  % left
            d = 0;
        case 'q'
            ex = 1;
    end
end

% start 
dt = 0.3;   
dy = 0.6;
while  ex ~= 1 
    % ball speedup
    [vel_ball,ang_ball] = uv2av(u_ball,v_ball);
    vel_ball = 1.001*vel_ball;
    [u_ball,v_ball] = av2uv(vel_ball,ang_ball);

    
    % move ball one step
    x_ball = x_ball + dt*u_ball;
    y_ball = y_ball + dt*v_ball;
    % reflections
    % top/bottom
    if y_ball > 10
        y_ball = 20 - y_ball;
        v_ball = -v_ball; 
    end
    if y_ball < -10
        y_ball = -20 - y_ball;
        v_ball = -v_ball; 
    end
        
    %% move right paddle 
    y_end = calcend(x_ball,y_ball,u_ball,v_ball);
    if y_r > y_end
        y_r = y_r - min([dy,abs(y_r-y_end)]);
    elseif y_r < y_end
        y_r = y_r + min([dy,abs(y_r-y_end)]);
    end
    % keep in bounds
    if y_r > ( 10-hpad), y_r = ( 10-hpad); end
    if y_r < (-10+hpad), y_r = (-10+hpad); end

    %% move left paddle 
    if strcmp(opt1,'a')
        y_end = calcend(x_ball,y_ball,u_ball,v_ball);
        if y_l > y_end
            y_l = y_l - min([dy,abs(y_l-y_end)]);
        elseif y_l < y_end
            y_l = y_l + min([dy,abs(y_l-y_end)]);
        end
    elseif strcmp(opt1,'p')
        dyp = d*dy + dyp;
        y_l = y_l + dyp;
        d = 0.5*d;
        dyp = 0.5*dyp;
%         if d < 0.1, d=0; end
    end
    
    % keep in bounds
    if y_l > ( 10-hpad), y_l = ( 10-hpad); end
    if y_l < (-10+hpad), y_l = (-10+hpad); end

    plot_sit(gca,x_ball,y_ball,y_l,y_r,dpad,hpad)
    drawnow
    
    %% paddle interaction
    % right
    if x_ball > (10-dpad)
        x_ball = 10;
        % location of ball hit on paddle (0 is in middle)
        f_r = (y_r - y_ball)/hpad;
        if abs(f_r) > 1
            fprintf('\nLEFT PLAYER WINS!\n\n');
            break % game over, r loses
        else
            % reflect at an angle
            [vel_ball,ang_ball] = uv2av(u_ball,v_ball);
            ang_ball = -ang_ball + 50*sind(f_r*90);
            if ang_ball > 0, ang_ball = -20; end
            if ang_ball <-180, ang_ball = -160; end
            [u_ball,v_ball] = av2uv(vel_ball,ang_ball);
        end
    end
    % left
    if x_ball < (-10+dpad)
        x_ball = -10;
        f_l = (y_l - y_ball)/hpad;
        if abs(f_l) > 1
            fprintf('\nRIGHT PLAYER WINS!\n\n');
            break % game over, l loses
        else
            % reflect at an angle
            [vel_ball,ang_ball] = uv2av(u_ball,v_ball);
            ang_ball = -ang_ball + 50*sind(f_l*90);
            if ang_ball < 0
                ang_ball = 20; 
            end
            if ang_ball > 180
                ang_ball = 160; 
            end
            [u_ball,v_ball] = av2uv(vel_ball,ang_ball);
        end
    end
    

pause(0.2./opt2)
    
end
plot_sit(gca,x_ball,y_ball,y_l,y_r,dpad,hpad)
drawnow

end


function plot_sit(ax,x_ball,y_ball,y_l,y_r,dpad,hpad)
llim = -10;
rlim = 10;
cla(ax), hold on
% plot ball
scatter(ax,x_ball,y_ball,100,'k','filled')
% plot left paddle
fill(llim+[0,-dpad,-dpad,0,0],y_l+hpad*[1,1,-1,-1,1],'k')
% plot right paddle
fill(rlim+[0,+dpad,+dpad,0,0],y_r+hpad*[1,1,-1,-1,1],'k')

set(ax,'xlim',[llim-dpad,rlim+dpad],'ylim',[-10,10],'box','on','linewidth',2,...
    'xtick',[],'ytick',[])
end

% compute predicted endpoint (stupidly - no reflection)
function y_end = calcend(x_ball,y_ball,u_ball,v_ball)
% gradient of ball
m_ball = v_ball/u_ball;
    if u_ball > 0 % moving right
        x_end = 10;
    elseif u_ball < 0 % moving left
        x_end = -10;
    end
y_end = y_ball + m_ball*(x_end-x_ball);
end

% ball speed,dir
function [vel_ball,ang_ball] = uv2av(u_ball,v_ball)
    vel_ball = sqrt(u_ball.^2 + v_ball.^2);
    ang_ball = atan2d(u_ball,v_ball);
end

% ball x-vel,y-vel
function [u_ball,v_ball] = av2uv(vel_ball,ang_ball)
     u_ball = vel_ball*sind(ang_ball);
     v_ball = vel_ball*cosd(ang_ball);
end




