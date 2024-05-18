clear
clc
tic
warning off

%% 1.设参
mpc = data; %读取数据
pload = mpc.Pload; %24h有功负荷
qload = mpc.Qload; %24h无功负荷
branch = mpc.branch; %支路参数
branch(:, 3) = branch(:, 3) * 10 / (12.66 ^ 2);
R = real(branch(:, 3)); %线路电阻
X = imag(branch(:, 3)); %线路电抗
P_WT = mpc.Pwt(:, 2:end); %风电有功出力
Q_WT = mpc.Qwt(:, 2:end); %风电无功出力
P_PV = mpc.Ppv(:, 2:end); %光伏有功出力
T = 24; %时段数为24小时
nb = 69; %节点数
nl = 68; %支路数
nwt = 5; %风机数
npv = 5; %光伏数
ness = 2; %ESS数
%节点连接方式
upstream = zeros(nb, nl);
dnstream = zeros(nb, nl);

for t = 1:nl
    upstream(t, t) = 1;
end

for t = [1:25, 27:33, 35:44, 46:48, 50, 52:63, 65, 67]
    dnstream(t, t + 1) = 1;
end

dnstream(2, 27) = 1;
dnstream(2, 35) = 1;
dnstream(3, 46) = 1;
dnstream(7, 50) = 1;
dnstream(8, 52) = 1;
dnstream(10, 65) = 1;
dnstream(11, 67) = 1;
dnstream(69, 1) = 1;
%电压上下限，节点1是平衡节点
Vmax = [1.1 * 1.1 * ones(nb - 1, T); ones(1, T)];
Vmin = [0.9 * 0.9 * ones(nb - 1, T); ones(1, T)];
%发电机出力上下限
Pgmax = [zeros(nb - 1, T); 1 * ones(1, T)];
Pgmin = [zeros(nb - 1, T); 0 * ones(1, T)];
Qgmax = [zeros(nb - 1, T); 1 * ones(1, T)];
Qgmin = [zeros(nb - 1, T); -1 * ones(1, T)];
%% 2.设变量
V = sdpvar(nb, T); %电压的平方
I = sdpvar(nl, T); %电流的平方
P = sdpvar(nl, T); %线路有功
Q = sdpvar(nl, T); %线路无功
Pg = sdpvar(nb, T); %发电机有功
Qg = sdpvar(nb, T); %发电机无功
p_wt = sdpvar(nwt, T); %风电有功
q_wt = sdpvar(nwt, T); %风电无功
p_pv = sdpvar(npv, T); %光伏有功
p_dch = sdpvar(ness, T); %ESS放电功率
p_ch = sdpvar(ness, T); %ESS充电功率
u_dch = binvar(ness, T); %ESS放电状态
u_ch = binvar(ness, T); %ESS充电状态
E_ess = sdpvar(ness, T); %ESS的电量
S_IL1 = sdpvar(1, T, 'full'); %可削减负荷
S_IL2 = sdpvar(1, T, 'full'); %可削减负荷
%% 3.设约束
C = [];
%% 3.1需求响应约束
C = [C, 0 <= S_IL1 <= 0.8 * pload(10, :)];
C = [C, 0 <= S_IL2 <= 0.8 * pload(26, :)];
IL1 = [zeros(9, T); S_IL1(1, :); zeros(59, T)]; %削减负荷
IL2 = [zeros(25, T); S_IL2(1, :); zeros(43, T)]; %削减负荷
%% 3.2储能装置（ESS）约束
%充放电状态约束
C = [C, u_dch(1, :) + u_ch(1, :) <= 1];
C = [C, u_dch(2, :) + u_ch(2, :) <= 1];
%功率约束
C = [C, 0 <= p_dch <= u_dch * 1];
C = [C, 0 <= p_ch <= u_ch * 1];
%容量约束
for t = 1:23
    C = [C, E_ess(1, t + 1) == E_ess(1, t) + p_ch(1, t) - p_dch(1, t)];
end

for t = 1:23
    C = [C, E_ess(2, t + 1) == E_ess(2, t) + p_ch(2, t) - p_dch(2, t)];
end

C = [C, 0.1 <= E_ess(1, :) <= 0.9];
C = [C, 0.1 <= E_ess(2, :) <= 0.9];
%投入节点选择
P_dch1 = [zeros(22, T); p_dch(1, :); zeros(46, T)];
P_ch1 = [zeros(22, T); p_ch(1, :); zeros(46, T)];
P_dch2 = [zeros(56, T); p_dch(2, :); zeros(12, T)];
P_ch2 = [zeros(56, T); p_ch(2, :); zeros(12, T)];
%% 3.3风机+光伏约束
C = [C, 0 <= p_wt, p_wt <= P_WT];
C = [C, 0 <= q_wt, q_wt <= Q_WT];
P_wt = [zeros(2, T); p_wt(1, :); zeros(15, T); p_wt(2, :); p_wt(3, :); zeros(6, T); p_wt(4, :); zeros(6, T); p_wt(5, :); zeros(35, T)];
Q_wt = [zeros(2, T); q_wt(1, :); zeros(15, T); q_wt(2, :); q_wt(3, :); zeros(6, T); q_wt(4, :); zeros(6, T); q_wt(5, :); zeros(35, T)];
C = [C, 0 <= p_pv, p_pv <= P_PV];
P_pv = [zeros(37, T); p_pv(1, :); zeros(9, T); p_pv(2, :); zeros(2, T); p_pv(3, :); zeros(2, T); p_pv(4, :); zeros(11, T); p_pv(5, :); zeros(3, T)];
%% 3.4潮流约束
Pin = -upstream * P + upstream * (I .* (R * ones(1, T))) + dnstream * P; %节点注入有功
Qin = -upstream * Q + upstream * (I .* (X * ones(1, T))) + dnstream * Q; %节点注入无功
C = [C, Pin + pload - IL1 - IL2 - Pg - P_wt - P_pv - P_dch1 + P_ch1 - P_dch2 + P_ch2 == 0];
C = [C, Qin + qload - Q_wt - Qg == 0];
%欧姆定律约束
C = [C, V(branch(:, 2), :) == V(branch(:, 1), :) - 2 * (R * ones(1, 24)) .* P - 2 * (X * ones(1, 24)) .* Q + ((R .^ 2 + X .^ 2) * ones(1, 24)) .* I];
%二阶锥约束
C = [C, V(branch(:, 1), :) .* I >= P .^ 2 + Q .^ 2];
%% 3.5通用约束
%节点电压约束
C = [C, Vmin <= V, V <= Vmax];
%发电机功率约束
C = [C, Pgmin <= Pg, Pg <= Pgmax, Qgmin <= Qg, Qg <= Qgmax];
%支路电流约束
C = [C, 0 <= I, I <= 11];
%% 4.设目标函数
f1 = sum(sum(I .* (R * ones(1, T)))); %网损
f2 = mean((1.6 - S_IL1 ./ pload(10, :) - S_IL2 ./ pload(26, :)) / 2) %负荷削减量
f3 = mean(mean((3 - p_wt ./ P_WT - q_wt ./ Q_WT - p_pv ./ P_PV) / 3)) %新能源消纳量
f = f1 + f2 + f3;
toc

%% 5.求解
ops = sdpsettings('solver', 'cplex');
sol = optimize(C, f, ops);
%分析错误
if sol.problem == 0
    disp('successful solved求解成功');
else
    disp('error');
    error = P .^ 2 + Q .^ 2 - V(Line(:, 1), :) .* I;
    yalmiperror(sol.problem)
end

%% 6.绘图
p_wt = value(p_wt);
q_wt = value(q_wt);
p_pv = value(p_pv);
p_dch = value(p_dch);
p_ch = value(p_ch);
E_ess = value(E_ess);
S_IL1 = value(S_IL1);
S_IL2 = value(S_IL2);

figure(1)
t = 1:1:24;
plot(t, sum(p_wt), '--k', 'linewidth', 2);
hold on
plot(t, sum(q_wt), '-*k', 'linewidth', 2);
hold on
plot(t, sum(p_pv), '--r', 'linewidth', 2);
xlabel('时刻/h');
ylabel('新能源出力')
legend('风电有功', '风电无功', '光伏有功');

figure(2)
t = 1:1:24;
plot(t, p_dch(1, :) - p_ch(1, :), '--', 'linewidth', 2);
hold on
plot(t, p_dch(2, :) - p_ch(2, :), '-*', 'linewidth', 2);
xlabel('时刻/h');
ylabel('储能出力')
legend('ESS1', 'ESS2');

figure(3)
t = 1:1:24;
plot(t, E_ess(1, :), '--', 'linewidth', 1.8);
hold on
plot(t, E_ess(2, :), '-*', 'linewidth', 1.8);
xlabel('时刻/h');
ylabel('SOC');
legend('ESS1', 'ESS2');

figure(4)
t = 1:1:24;
plot(t, S_IL1, '--', 'linewidth', 2);
hold on
plot(t, S_IL2, '-*', 'linewidth', 2);
xlabel('时刻/h');
ylabel('负荷削减');
legend('可削减负荷1', '可削减负荷2');
