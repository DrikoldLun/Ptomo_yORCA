function [out_pos]=fun_mm_subplot_pos(n_row,n_col,sbp_width,sbp_heig);

% subplot pos
% input : -----------------
% n_row     : number of row
% n_col     : number of column
% sbp_width : subplot width  0-1
% sbp_heig  : subplot height 0-1
% output : ----------------
% out_pos
%
%
% sbp_width=0.8;
% sbp_heig=0.85;
% n_row = 4;
% n_col = 6;
% [out_pos]=fun_mm_subplot_pos(n_row,n_col,sbp_width,sbp_heig);
% h=figure(1)
% set(h, 'Position', [100, 100, 800, 600]);
% [x,y,z]=peaks(30);
% caxis_mm=[-5 30];
% for i=1:n_row*n_col
%     ax=axes('position',out_pos(i,:));
%     surf(x,y,z+i,'edgecolor','none')
%     caxis(caxis_mm)
%     view(0,90)
%     hold on
%     text([2],[2],[30],mat2str(i),'fontsize',20)
% end

xmax = 0.85;
for i=1:n_row*n_col

    k1=ceil(i/n_col);
    k2=i-(k1-1)*n_col;
   
    out_pos(i,:)=[((xmax-sbp_width)/n_col/1.5 + (k2-1)/n_col),...
                  ((1-sbp_heig)/n_row/1.5 + (n_row-k1)/n_row), ...
                  xmax*sbp_width/n_col,...
                  sbp_heig/n_row];
end