function [Shortest_Length,answer] = ant(Data,flag)
clc ;
%% 初始化参数
N=11;    %% m 蚂蚁个数
alpha=1; %% 信息素重要程度
beta=5;  %% 启发式因子重要程度
rho=0.1; %% 信息素蒸发系数
max_iter=200; %%最大迭代次数
Q=100;        %%信息素增强系数
n=size(Data,1);
road_length=zeros(n,n);%D表示完全图的赋权邻接矩阵

for i=1:n
    for j=1:n
        if i~=j
            road_length(i,j)=Distance(Data(i,1),Data(i,2),Data(j,1),Data(j,2));
        else安全·13   ·32
            road_length(i,j)=eps;% 取倒数时使用
        end
        road_length(j,i)=road_length(i,j);   %对称矩阵
    end
end

%% 初始化路线和信息素矩阵
heu=1./road_length;                 %为启发因子
pheromoneMatrix=ones(n,n);          %信息素矩阵
path=zeros(N,n);   
iter=1;                             %迭代次数
path_best=zeros(max_iter,n);        %迭代中的最佳路线
length_best=inf.*ones(max_iter,1);  %迭代中的最佳路线的长度
length_mean=zeros(max_iter,1);      %迭代中的路线的平均长度

%% 迭代循环；停止条件：达到最大的迭代次数
dots = 0;
while iter<=max_iter        
    % 随机放蚂蚁
    positionInit=[];
    for i=1:(ceil(N/n))
        positionInit = [positionInit,randperm(n)];
    end
    path(:,1)=(positionInit(1,1:N))';  
    
    % 蚂蚁随机选择目的地
    for j=2:n     
        for i=1:N
            visited=path(i,1:(j-1));  
            pos=zeros(1,(n-j+1));
            P=pos;
            pass_cities=1;
            for k=1:n
                if isempty(find(visited==k, 1))
                    pos(pass_cities)=k;
                    pass_cities=pass_cities+1;
                end
            end
            for k=1:length(pos)
                P(k)=(pheromoneMatrix(visited(end),pos(k))^alpha)*(heu(visited(end),pos(k))^beta);
            end
            P=P/(sum(P));
            Psum=cumsum(P);
            Select=find(Psum>=rand);
            to_visit=pos(Select(1));
            path(i,j)=to_visit;
        end
    end
    if iter>=2
        path(1,:)=path_best(iter-1,:);
    end
    
    % 初始化总距离矩阵
    L=zeros(N,1);
    for i=1:N
        P=path(i,:);
        for j=1:(n-1)
            % 距离累和
            L(i)=L(i)+road_length(P(j),P(j+1));
        end
        % 总距离
        L(i)=L(i)+road_length(P(1),P(n));
    end
    length_best(iter)=min(L);
    pos=find(L==length_best(iter));
    path_best(iter,:)=path(pos(1),:);
    length_mean(iter)=mean(L);
    delta_pheromone=zeros(n,n);
    fprintf('.');
    dots = dots + 1;
    for i=1:N
        for j=1:(n-1)
            delta_pheromone(path(i,j),path(i,j+1))=delta_pheromone(path(i,j),path(i,j+1))+Q/L(i);
        end
        
        delta_pheromone(path(i,n),path(i,1))=delta_pheromone(path(i,n),path(i,1))+Q/L(i);
    end
    % 信息素更新
    pheromoneMatrix=(1-rho).*pheromoneMatrix+delta_pheromone;
    path=zeros(N,n);
    
    iter = iter+1;
    if dots>78
        dots = 0;
        fprintf('\n');
    end
end
fprintf('Done!\n');
%% 最终路径
Pos=find(length_best==min(length_best)); 
Shortest_Route=path_best(Pos(1),:);
Shortest_Length=length_best(Pos(1));

%% 数据展示
disp('最优路径为');
Shortest_Route = Shortest_Route(1:end,:)-1;
for i=1:length(Shortest_Route)
    if flag==1 && Shortest_Route(i)==0 
            fprintf('数据中心');
            if i~=length(Shortest_Route)
                fprintf('->');
            else
                fprintf('\n');
            end
    else
        fprintf('站点%d',Shortest_Route(i));
        if i~=length(Shortest_Route)
            fprintf('->');
        else
            fprintf('\n');
        end
    end
end
fprintf('长度为 %d/m\n',Shortest_Length);


%% 可视化
figure(1);
plot(length_best,'k');
xlabel('迭代次数');
ylabel('目标函数值');
title('适应度的进化曲线');


figure(2)
N=length(P);
scatter(Data(:,1),Data(:,2),'r');
if flag==1
    for i = 1:length(Data)
        if i==1 
            text(Data(i,1),Data(i,2),'(数据中心)');
        else
            text(Data(i,1),Data(i,2),['(站点' num2str(i-1) ')']);
        end
    end
else
    for i = 1:length(Data)
        text(Data(i,1),Data(i,2),['(站点' num2str(i) ')']);
    end
end
 hold on
 plot([Data(P(1),1),Data(P(N),1)],[Data(P(1),2),Data(P(N),2)],'k:')
 hold on
 
 answer = path_best(end,:);
 % 画图记得去掉数据中心
for i=1:N
    j=i+1;
    if(i+1>N)
        j=1;
    end
    plot([Data(answer(i),1),Data(answer(j),1)],[Data(answer(i),2),Data(answer(j),2)],'k:')
end
save best_route.mat path_best;
hold on
xlabel('经度');
ylabel('纬度');
title('最优路线规划 ')
grid on
figure(3)
plot(length_best,'b')
hold on                         %保持图形
plot(length_mean,'k')
title('平均距离和最短距离')     %标题
save length_mean.mat length_mean;
save length_best.mat length_best;
