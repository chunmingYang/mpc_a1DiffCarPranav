clear all
clc
close all
filename = 'data.csv';

T = readtable(filename); %check T.Properties
VariableNames = T.Properties.VariableNames;
Arr = table2array(T);
[m,n] = size(Arr);

%%%% Print Variable names
for i=1:length(VariableNames)
    disp([mat2str(i), ' = ',VariableNames{i}]);
end

for i=1:length(VariableNames)
    tmp_arr = Arr(:,i);
    eval([VariableNames{i},'=tmp_arr;']);
end

index_start = 1;
index_end = length(t);

%code to plot for a given time windows t_start to t_end
%code to plot for a given time windows t_start to t_end
t_start = 5;
t_end = t(end);
[index_start, ~] = get_time_indices(t,t_start,t_end);


index = index_start:index_end;
fig_no = 1;

figure(fig_no); fig_no = fig_no + 1;
plot(WorldPosX,WorldPosY,'k','Linewidth',2);
xlabel('x');
ylabel('y');
axis('equal')

% figure(fig_no); fig_no = fig_no + 1;
% plot(WorldPosX,WorldPosY,'k','Linewidth',2);
% xlabel('x');
% ylabel('y');
% axis('equal')

% figure(fig_no); fig_no = fig_no + 1;
% subplot(2,1,1)
% plot(t(index),error_x(index),'k','Linewidth',2);
% ylabel('error x');
% subplot(2,1,2)
% plot(t(index),error_y(index),'k','Linewidth',2);
% ylabel('error y');
% xlabel('time');


% for i=1:length(t)
%     body_vel = rotation(eulerZ(i,1))'*[WorldVelX(i,1); WorldVelY(i,1)];
%     BodyVelX(i,1) = body_vel(1);
%     BodyVelY(i,1) = body_vel(2);
% end
% figure(fig_no); fig_no = fig_no + 1;
% subplot(3,1,1)
% plot(t(index),BodyVelX(index));
% xlabel('t');
% ylabel('BodyVelX');
% subplot(3,1,2)
% plot(t(index),BodyVelY(index));
% xlabel('t');
% ylabel('BodyVelY');
% subplot(3,1,3)
% plot(t(index),WorldVelZ(index));
% xlabel('t');
% ylabel('BodyVelZ');



function R = rotation(theta)
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
end

function [index_start, index_end] = get_time_indices(t,t_start,t_end) 
    for i=1:length(t)
        if (t(i)< t_start)
            index_start = i;
        else
            break;
        end
    end

    index_end = [];
    if (nargout>2)
        for i=length(t):-1:1
            if (t(i) > t_end)
                index_end = i;
            else
                break;
            end
        end
    end
end

% function index = get_index(VariableNames,index_names)
%     index = [];
%     for i=1:length(index_names)
%         for j=1:length(VariableNames)
%             if (strcmp(index_names{i},VariableNames{j}))
%                 index = [index,j];
%                 break;
%             end
%         end
%     end
% end
