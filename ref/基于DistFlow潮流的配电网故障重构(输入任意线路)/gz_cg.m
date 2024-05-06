%单时段网络重构socp-opf
clear 
clc
tic
warning off
%% 1.设参
mpc = IEEE33BW;
pload = mpc.Pload(:,1);%节点有功负荷
qload = mpc.Qload(:,1);%节点无功负荷
branch = mpc.branch_CG;
branch(:,3) = branch(:,3)*100/(12.66^2);%求阻抗标幺值
r=real(branch(:,3));
x=imag(branch(:,3));

T = 1;%时段数为1小时
nb = 33;%节点数,根节点为33
nl = 37;%支路数
nc = 5;%联络开关数
AA=input('输入故障线路编号（1-32）:');
upstream=zeros(nb,nl);%代表流入节点支路
dnstream=zeros(nb,nl);%代表流出节点支路
for i=1:32
    upstream(i,i)=1;
end
upstream(20,33)=1;%支路33为20-7支路，流入节点20
upstream(14,34)=1;%支路34为14-8支路，流入节点14
upstream(21,35)=1;%支路35为21-11支路，流入节点21
upstream(32,36)=1;%支路36为32-17支路，流入节点32
upstream(28,37)=1;%支路37为28-24支路，流入节点28

for i=[1:16,18:20,22:23,25:31]
    dnstream(i,i+1)=1;
end
dnstream(1,18)=1;
dnstream(2,22)=1;
dnstream(5,25)=1;
dnstream(33,1)=1;
%5条流入，对应5条流出
dnstream(7,33)=1;
dnstream(8,34)=1;
dnstream(11,35)=1;
dnstream(17,36)=1;
dnstream(24,37)=1;

Vmax=[1.06*1.06*ones(32,1);1.06*1.06*ones(1,1)];
Vmin=[0.94*0.94*ones(32,1);1.06*1.06*ones(1,1)];
Pgmax=[zeros(32,1);ones(1,1)];
Qgmax=[zeros(32,1);ones(1,1)];

%% 2.设变量
V = sdpvar(nb,T);%电压的平方
I = sdpvar(nl,T);%电流的平方
P = sdpvar(nl,T);%线路有功
Q = sdpvar(nl,T);%线路无功
Pg = sdpvar(nb,T);%发电机有功
Qg = sdpvar(nb,T);%发电机无功
lamda = sdpvar(nb,T);%弃负荷比例
Zij=binvar(nl,1);%网架结构
Z0=[ones(nl-nc,1);zeros(nc,1)];%初始拓扑
assign(Zij,Z0);

%% 3.设约束
Constraints = [];
%% 网络重构约束
Constraints = [Constraints, sum(Zij) == 32];
Constraints = [Constraints, I(AA) == 0];
Constraints = [Constraints, 0<=lamda<=1];
% P_tree = sdpvar(37,1);%虚拟有功
% Pin_tree = -upstream*P_tree + dnstream*P_tree;%虚拟节点注入有功
% Constraints = [Constraints,-Zij <= P_tree <= Zij];
% Constraints = [Constraints, Pin_tree(1:32) + 0.01==0];
%% 潮流约束
%节点功率约束
Pin = -upstream*P + upstream*(I.*(r*ones(1,T))) + dnstream*P;%节点注入有功
Qin = -upstream*Q + upstream*(I.*(x*ones(1,T))) + dnstream*Q;%节点注入无功
Constraints = [Constraints, Pin + (1-lamda).*pload - Pg ==0];
Constraints = [Constraints, Qin + qload - Qg==0];
%欧姆定律约束
m = 1.06*1.06 - 0.94*0.94;
M = (ones(nl,1) - Zij)*m;
Constraints = [Constraints, V(branch(:,1),:) - V(branch(:,2),:) <= M + 2*(r).*P + 2*(x).*Q - ((r.^2 + x.^2)).*I];
Constraints = [Constraints, V(branch(:,1),:) - V(branch(:,2),:) >= -M + 2*(r).*P + 2*(x).*Q - ((r.^2 + x.^2)).*I];
%二阶锥约束
Constraints = [Constraints, V(branch(:,1),:).*I >= P.^2+Q.^2];
%% 通用约束
%节点电压约束
Constraints = [Constraints, Vmin <= V,V <= Vmax];
%发电机功率约束
Constraints = [Constraints, -Pgmax <= Pg,Pg <= Pgmax,-Qgmax <= Qg,Qg <= Qgmax];
%支路电流约束
Constraints = [Constraints, 0 <= I,I <= 0.1*Zij];
%支路功率约束
Constraints = [Constraints, -0.1*Zij <= P,P <= 0.1*Zij,-0.1*Zij <= Q,Q <= 0.1*Zij,P(AA) == 0,Q(AA) == 0];
%虚拟网络%连续性和辐射性约束
N=1;
Fij=sdpvar(37,N,'full'); 
% Wj=sdpvar(3,N,'full'); 
M=100;
for t=1:N
    for k=2:33
            node_out=find(branch(:,1)==k);
            node_in=find(branch(:,2)==k);
            Constraints=[Constraints,sum(Fij(node_in,t))-sum(Fij(node_out,t))==-1];  
    end

end
Constraints=[Constraints,-M.*Zij<=Fij<=M.*Zij];
Constraints=[Constraints,-M.*(2-Zij)<=Fij<=M.*(2-Zij)];

%% 4.设目标函数
% objective = -mean(V);%
objective = sum(I.*(r*ones(1,T)))+sum(lamda.*pload);
toc%建模时间
%% 5.设求解器
ops=sdpsettings('verbose', 1, 'solver', 'cplex');
sol=optimize(Constraints,objective,ops);
toc%求解时间
objective=100*value(objective)
Zij = -value(Zij)
%% 6.输出AMPL模型
% saveampl(Constraints,objective,'mymodel');

%% 7.分析错误标志
if sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end
%% 8.保存运行结果到XLSX
% Doc_name='OPF_data';
% V = value(V);I = value(I);P = value(P);Q = value(Q);
% Pg = value(Pg);Qg = value(Qg);AA = value(AA);
% xlswrite(Doc_name,V,'V');
% xlswrite(Doc_name,I,'I');
% xlswrite(Doc_name,P,'P');
% xlswrite(Doc_name,Q,'Q');
% xlswrite(Doc_name,Pg,'Pg');
% xlswrite(Doc_name,Qg,'Qg');
% xlswrite(Doc_name,AA,'AA');
A1=[1,3];A2=[2,3];A3=[3,3];A4=[4,3];A5=[5,3];A6=[6,3];
A7=[7,3];A8=[8,3];A9=[9,3];A10=[10,3];A11=[11,3];A12=[12,3];
A13=[13,3];A14=[14,3];A15=[15,3];A16=[16,3];A17=[17,3];A18=[18,3];
A19=[2,1];A20=[3,1];A21=[4,1];A22=[5,1];A23=[3,5];A24=[4,5];
A25=[5,5];A26=[6,5];A27=[7,5];A28=[8,5];A29=[9,5];A30=[10,5];
A31=[11,5];A32=[12,5];A33=[13,5];A34=[9,2];A35=[15,2];A36=[5,6];A37=[9,6];
for p = 1:37  % 放进一个矩阵内方便操作
    c = num2str(p);
    s = ['A(p,:) = A' c ';'];
    eval(s);
end
% scatter(A(:,1),A(:,2));
for p = 1:33
    c = num2str(p);
    plot(A(p,1),A(p,2),'ko','MarkerFaceColor','k');
    hold on
    axis([0 20 0 7]);
    text(A(p,1)+0.05,A(p,2)-0.1,c);
end
set(gca,'xtick',0:1:20,'ytick',0:1:7);

sx = 1:18;
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k');
end

sx = [2 19 20 21 22];
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k');
end

sx = [3 23 24 25];
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k');
end
sx = [6 26 27 28 29 30 31 32 33];
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k');
end
if AA~=18
plot(A(AA,1)+0.5,A(AA,2),'p','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5);
else
    plot(A(19,1),2,'p','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5);
end

sx = [21 8];%33
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k','linestyle','--');
end
sx = [9 34 35 15];%34
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k','linestyle','--');
end
sx = [12 22];%35
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k','linestyle','--');
end
sx = [18 33];%36
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k','linestyle','--');
end
sx = [25 36 37 29];%37
for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','k','linestyle','--');
end
grid on
for p=1:17
    if Zij(p)==0
        sx=1:18;
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);  
     line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
    end
%plot(A(jieru(p),1),A(jieru(p),2),'p','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1.5);
end
for p=1:4
    sx = [2 19 20 21 22];
    if Zij(p+17)==0
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);  
     line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
    end
end
for p=1:3
    sx = [3 23 24 25];
    if Zij(p+21)==0
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);  
     line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
    end
end
for p=1:8
sx = [6 26 27 28 29 30 31 32 33];
    if Zij(p+24)==0
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);  
     line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
    end
end
if Zij(33)==0
 sx = [21 8];%33
  for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
  end
end
if Zij(34)==0
    sx = [9 34 35 15];%34
  for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
  end
end
if Zij(35)==0
    sx = [12 22];%35
  for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
  end
end
if Zij(36)==0
    sx = [18 33];%36
  for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
  end
end
if Zij(37)==0
sx = [25 36 37 29];%37
  for p = 1:length(sx)-1
    P1 = A(sx(p),:);
    P2 = A(sx(p+1),:);
    line([P1(1) P2(1)],[P1(2) P2(2)],'color','r');
  end
end
hold off
title('故障处理网络重构图');