clear
clc
tic
warning off

%% 0.注释
% Pre:预测值（只有负荷和新能源出力上限的数据具有波动性）
% Act:实际值
% Del:修正值Δ
% UP(Unprocessed):未处理的
% P_Line:有功
% Q_Line:无功
% wt:风电
% pv:光伏
% ess:储能
% num:数量
% plan:基于预测数据的计划值
% fluc:基于实际数据（波动）的实际值

%% 1.设参
%读取数据
mpc = data;
%24h有功/无功负荷
Pre_P_Load = mpc.Pre_P_Load;
Pre_Q_Load = mpc.Pre_Q_Load;
%支路参数
branch = mpc.branch;
branch(:, 3) = branch(:, 3) * mpc.baseMVA / (mpc.Vb ^ 2);
R = real(branch(:, 3)); %线路电阻
X = imag(branch(:, 3)); %线路电抗
%RES参数
Pre_P_wt = mpc.Pre_P_wt(:, 2:end); %风电有功
Pre_Q_wt = mpc.Pre_Q_wt(:, 2:end); %风电无功
Pre_P_pv = mpc.Pre_P_pv(:, 2:end); %光伏有功
%时段与机组数
T = 24; %时段数为24小时
num_Bus = length(Pre_P_Load); %节点数
num_Line = length(branch); %支路数
num_wt = size(Pre_P_wt, 1); %风机数
num_pv = size(Pre_P_pv, 1); %光伏数
num_ess = 2; %ESS数
num_load_sd = 2; %可削减负荷数
%储能系统
P_ess_ch_Max = repmat(mpc.ess(:, 2), 1, T);
P_ess_dch_Max = repmat(mpc.ess(:, 3), 1, T);
storage_ess = repmat(mpc.ess(:, 4), 1, T);
%目标函数系数
eta_loss_plan = mpc.eta_loss_plan;
eta_RES_plan = mpc.eta_RES_plan;
%节点连接方式
upstream = zeros(num_Bus, num_Line);
dnstream = zeros(num_Bus, num_Line);
for t = 1:num_Line
    upstream(t, t) = 1;
    dnstream(branch(t, 1), branch(t, 2)) = 1;
end
%电压上下限，节点69是平衡节点
V_Bus_Max = [1.1 * 1.1 * ones(num_Bus - 1, T); ones(1, T)];
V_Bus_Min = [0.9 * 0.9 * ones(num_Bus - 1, T); ones(1, T)];
%发电机出力上下限
P_Gen_Max = [zeros(num_Bus - 1, T); 1 * ones(1, T)];
P_Gen_Min = [zeros(num_Bus - 1, T); 0 * ones(1, T)];
Q_Gen_Max = [zeros(num_Bus - 1, T); 1 * ones(1, T)];
Q_Gen_Min = [zeros(num_Bus - 1, T); -1 * ones(1, T)];
%% 2.设变量
%电压电流（平方）
V = sdpvar(num_Bus, T); %电压的平方
I = sdpvar(num_Line, T); %电流的平方
%线路
P_Line = sdpvar(num_Line, T); %线路有功
Q_Line = sdpvar(num_Line, T); %线路无功
%发电机
P_Gen = sdpvar(num_Bus, T); %发电机有功
Q_Gen = sdpvar(num_Bus, T); %发电机无功
%新能源
UP_P_wt = sdpvar(num_wt, T); %未矩阵处理的风电实际有功出力
UP_Q_wt = sdpvar(num_wt, T); %未矩阵处理的风电实际无功出力
UP_P_pv = sdpvar(num_pv, T); %未矩阵处理的光伏实际有功出力
P_wt = zeros(num_Bus - num_wt, T); %将风电有功出力变量放入矩阵
Q_wt = zeros(num_Bus - num_wt, T); %将风电无功出力变量放入矩阵
P_pv = zeros(num_Bus - num_pv, T); %将光伏有功出力变量放入矩阵
%储能
UP_P_ess_dch = sdpvar(num_ess, T); %ESS放电功率
UP_P_ess_ch = sdpvar(num_ess, T); %ESS充电功率
U_ess_dch = binvar(num_ess, T); %ESS放电状态
U_ess_ch = binvar(num_ess, T); %ESS充电状态
SOC_ess = sdpvar(num_ess, T); %ESS的电量
P_ess_dch = zeros(num_Bus - num_ess, T);
P_ess_ch = zeros(num_Bus - num_ess, T);
%可削减负荷
UP_P_load_sd = sdpvar(num_load_sd, T, 'full');
P_load_sd = zeros(num_Bus - num_load_sd, T);
%节点功率
P_bus_in = -upstream * P_Line + upstream * (I .* (R * ones(1, T))) + dnstream * P_Line; %节点注入有功
Q_bus_in = -upstream * Q_Line + upstream * (I .* (X * ones(1, T))) + dnstream * Q_Line; %节点注入无功
%% 3.设约束
C_plan = [];
%% 3.1储能装置（ESS）约束
%充放电状态约束
C_plan = [C_plan, U_ess_dch(1:2, :) + U_ess_ch(1:2, :) <= 1];
%充放电功率约束
C_plan = [C_plan, 0 <= UP_P_ess_ch <= U_ess_ch .* P_ess_ch_Max];
C_plan = [C_plan, 0 <= UP_P_ess_dch <= U_ess_dch .* P_ess_dch_Max];
%容量约束
for t = 1:23
    C_plan = [C_plan, SOC_ess(1, t + 1) == SOC_ess(1, t) + UP_P_ess_ch(1, t) - UP_P_ess_dch(1, t)];
    C_plan = [C_plan, SOC_ess(2, t + 1) == SOC_ess(2, t) + UP_P_ess_ch(2, t) - UP_P_ess_dch(2, t)];
end
C_plan = [C_plan, 0.1 * storage_ess <= SOC_ess <= 0.9 * storage_ess];
C_plan = [C_plan, SOC_ess(1:2, 1) == SOC_ess(1:2, 24)];
%投入节点选择
for i = 1:num_ess
    P_ess_dch = [P_ess_dch(1:mpc.ess(i, 1) - 1, :); UP_P_ess_dch(i, :); P_ess_dch(mpc.ess(i, 1):end, :)];
    P_ess_ch = [P_ess_ch(1:mpc.ess(i, 1) - 1, :); UP_P_ess_ch(i, :); P_ess_ch(mpc.ess(i, 1):end, :)];
end
%% 3.2风机+光伏约束
%出力功率约束
C_plan = [C_plan, 0 <= UP_P_wt, UP_P_wt <= Pre_P_wt];
C_plan = [C_plan, 0 <= UP_Q_wt, UP_Q_wt <= Pre_Q_wt];
C_plan = [C_plan, 0 <= UP_P_pv, UP_P_pv <= Pre_P_pv];
%投入节点选择
for i = 1:num_wt
    P_wt = [P_wt(1:mpc.Pre_P_wt(i, 1) - 1, :); UP_P_wt(i, :); P_wt(mpc.Pre_P_wt(i, 1):end, :)];
    Q_wt = [Q_wt(1:mpc.Pre_Q_wt(i, 1) - 1, :); UP_Q_wt(i, :); Q_wt(mpc.Pre_Q_wt(i, 1):end, :)];
end
for i = 1:num_pv
    P_pv = [P_pv(1:mpc.Pre_P_pv(i, 1) - 1, :); UP_P_pv(i, :); P_pv(mpc.Pre_P_pv(i, 1):end, :)];
end
%% 3.3潮流约束
%功率平衡约束
C_plan = [C_plan, P_bus_in + Pre_P_Load - P_Gen - P_wt - P_pv - P_ess_dch + P_ess_ch == 0];
C_plan = [C_plan, Q_bus_in + Pre_Q_Load - Q_wt - Q_Gen == 0];
%欧姆定律约束
C_plan = [C_plan, V(branch(:, 2), :) == V(branch(:, 1), :) - 2 * (R * ones(1, 24)) .* P_Line - 2 * (X * ones(1, 24)) .* Q_Line + ((R .^ 2 + X .^ 2) * ones(1, 24)) .* I];
%二阶锥约束
C_plan = [C_plan, V(branch(:, 1), :) .* I >= P_Line .^ 2 + Q_Line .^ 2];
%% 3.4通用约束
%节点电压约束
C_plan = [C_plan, V_Bus_Min <= V <= V_Bus_Max];
%发电机功率约束
C_plan = [C_plan, P_Gen_Min <= P_Gen <= P_Gen_Max, Q_Gen_Min <= Q_Gen <= Q_Gen_Max];
%支路电流约束
C_plan = [C_plan, 0 <= I, I <= 11];
%% 4.设目标函数
f_loss_plan = sum(sum(I .* (R * ones(1, T)))); %网损
f_RES_plan = mean(mean((3 - UP_P_wt ./ Pre_P_wt - UP_Q_wt ./ Pre_Q_wt - UP_P_pv ./ Pre_P_pv) / 3)); %新能源消纳量
f = eta_loss_plan * f_loss_plan + eta_RES_plan * f_RES_plan;
toc
%% 5.求解
%设求解器
ops = sdpsettings('solver', 'cplex');
sol = optimize(C_plan, f, ops);
%分析错误
if sol.problem == 0
    disp('successful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end
%% 6.画图
UP_P_wt_plan = value(UP_P_wt);
UP_Q_wt_plan = value(UP_Q_wt);
UP_P_pv_plan = value(UP_P_pv);
P_wt_plan = value(P_wt);
Q_wt_plan = value(Q_wt);
P_pv_plan = value(P_pv);
UP_P_ess_dch_plan = value(UP_P_ess_dch);
UP_P_ess_ch_plan = value(UP_P_ess_ch);
SOC_ess_plan = value(SOC_ess);
P_Gen_plan = value(P_Gen);
Q_Gen_plan = value(Q_Gen);
P_Line_plan = value(P_Line);
%% 7.波动模型设参
P_Gen_plan = value(P_Gen);
areaA = mpc.areaA;
areaB = mpc.areaB;
areaC = mpc.areaC;
areaD = mpc.areaD;
RES_cspt_all = mean(sum(UP_P_wt_plan + UP_Q_wt_plan + UP_P_pv_plan)); %计划新能源消纳量
RES_cspt_A = mean(sum(P_wt_plan(areaA, :) + Q_wt_plan(areaA, :) + P_pv_plan(areaA, :)));
RES_cspt_B = mean(sum(P_wt_plan(areaB, :) + Q_wt_plan(areaB, :) + P_pv_plan(areaB, :)));
RES_cspt_C = mean(sum(P_wt_plan(areaC, :) + Q_wt_plan(areaC, :) + P_pv_plan(areaC, :)));
RES_cspt_D = mean(sum(P_wt_plan(areaD, :) + Q_wt_plan(areaD, :) + P_pv_plan(areaD, :)));
RES_cspt = repmat([RES_cspt_A; RES_cspt_B; RES_cspt_C; RES_cspt_D], 1, T);
Act_P_Load = mpc.Act_P_Load;
Act_Q_Load = mpc.Act_Q_Load;
Act_P_wt = mpc.Act_P_wt(:, 2:end); %风电有功
Act_Q_wt = mpc.Act_Q_wt(:, 2:end); %风电无功
Act_P_pv = mpc.Act_P_pv(:, 2:end); %光伏有功
%目标函数系数
eta_loss_fluc = mpc.eta_loss_fluc;
eta_RES_fluc = mpc.eta_RES_fluc;
eta_sd_fluc = mpc.eta_sd_fluc;
eta_ess_fluc = mpc.eta_ess_fluc;
%% 8.波动模型设变量
P_wt = zeros(num_Bus - num_wt, T); %将风电有功出力变量放入矩阵
Q_wt = zeros(num_Bus - num_wt, T); %将风电无功出力变量放入矩阵
P_pv = zeros(num_Bus - num_pv, T); %将光伏有功出力变量放入矩阵
P_ess_dch = zeros(num_Bus - num_ess, T);
P_ess_ch = zeros(num_Bus - num_ess, T);
P_FCE = P_Gen - P_Gen_plan;
UP_P_FCE = [abs(P_Gen(69, :) - P_Gen_plan(69, :)) + abs(P_Line(27, :) - P_Line_plan(27, :)) + abs(P_Line(36, :) - P_Line_plan(36, :)); abs(P_Line(27, :) - P_Line_plan(27, :)); abs(P_Line(36, :) - P_Line_plan(36, :)); abs(P_Line(10, :) - P_Line_plan(10, :))]
%% 9.波动模型设约束
C_fluc = [];
%% 9.1需求响应约束
C_fluc = [C_fluc, 0 <= UP_P_load_sd <= 0.8 * Act_P_Load(mpc.load_sd(:, 1), :)];
for i = 1:num_load_sd
    P_load_sd = [P_load_sd(1:mpc.load_sd(i, 1) - 1, :); UP_P_load_sd(i, :); P_load_sd(mpc.load_sd(i, 1):end, :)];
end
%% 9.2储能装置（ESS）约束
%充放电状态约束
C_fluc = [C_fluc, U_ess_dch(1:2, :) + U_ess_ch(1:2, :) <= 1];
%充放电功率约束
C_fluc = [C_fluc, 0 <= UP_P_ess_ch <= U_ess_ch .* P_ess_ch_Max];
C_fluc = [C_fluc, 0 <= UP_P_ess_dch <= U_ess_dch .* P_ess_dch_Max];
%容量约束
for t = 1:23
    C_fluc = [C_fluc, SOC_ess(1, t + 1) == SOC_ess(1, t) + UP_P_ess_ch(1, t) - UP_P_ess_dch(1, t)];
end
for t = 1:23
    C_fluc = [C_fluc, SOC_ess(2, t + 1) == SOC_ess(2, t) + UP_P_ess_ch(2, t) - UP_P_ess_dch(2, t)];
end
C_fluc = [C_fluc, 0.1 * storage_ess <= SOC_ess <= 0.9 * storage_ess];
C_fluc = [C_fluc, SOC_ess(1:2, 1) == SOC_ess(1:2, 24)]
%投入节点选择
for i = 1:num_ess
    P_ess_dch = [P_ess_dch(1:mpc.ess(i, 1) - 1, :); UP_P_ess_dch(i, :); P_ess_dch(mpc.ess(i, 1):end, :)];
    P_ess_ch = [P_ess_ch(1:mpc.ess(i, 1) - 1, :); UP_P_ess_ch(i, :); P_ess_ch(mpc.ess(i, 1):end, :)];
end
%% 9.3风电+光伏约束
C_fluc = [C_fluc, 0 <= UP_P_wt, UP_P_wt <= Act_P_wt];
C_fluc = [C_fluc, 0 <= UP_Q_wt, UP_Q_wt <= Act_Q_wt];
for i = 1:num_wt
    P_wt = [P_wt(1:mpc.Act_P_wt(i, 1) - 1, :); UP_P_wt(i, :); P_wt(mpc.Act_P_wt(i, 1):end, :)];
    Q_wt = [Q_wt(1:mpc.Act_Q_wt(i, 1) - 1, :); UP_Q_wt(i, :); Q_wt(mpc.Act_Q_wt(i, 1):end, :)];
end
C_fluc = [C_fluc, 0 <= UP_P_pv, UP_P_pv <= Act_P_pv];
for i = 1:num_pv
    P_pv = [P_pv(1:mpc.Act_P_pv(i, 1) - 1, :); UP_P_pv(i, :); P_pv(mpc.Act_P_pv(i, 1):end, :)];
end
C_fluc = [C_fluc, mean(sum(UP_P_wt + UP_Q_wt + UP_P_pv)) >= min(mean(sum(Act_P_wt + Act_Q_wt + Act_P_pv)), RES_cspt_all)]
%% 9.4潮流约束
P_bus_in = -upstream * P_Line + upstream * (I .* (R * ones(1, T))) + dnstream * P_Line; %节点注入有功
Q_bus_in = -upstream * Q_Line + upstream * (I .* (X * ones(1, T))) + dnstream * Q_Line; %节点注入无功
C_fluc = [C_fluc, P_bus_in + Act_P_Load - P_load_sd - P_Gen - P_wt - P_pv - P_ess_dch + P_ess_ch == 0];
C_fluc = [C_fluc, Q_bus_in + Act_Q_Load - Q_wt - Q_Gen == 0];
%欧姆定律约束
C_fluc = [C_fluc, V(branch(:, 2), :) == V(branch(:, 1), :) - 2 * (R * ones(1, 24)) .* P_Line - 2 * (X * ones(1, 24)) .* Q_Line + ((R .^ 2 + X .^ 2) * ones(1, 24)) .* I];
%二阶锥约束
C_fluc = [C_fluc, V(branch(:, 1), :) .* I >= P_Line .^ 2 + Q_Line .^ 2];
%% 9.5通用约束
%节点电压约束
C_fluc = [C_fluc, V_Bus_Min <= V, V <= V_Bus_Max];
%发电机功率约束
C_fluc = [C_fluc, P_Gen_Min <= P_Gen, P_Gen <= P_Gen_Max, Q_Gen_Min <= Q_Gen, Q_Gen <= Q_Gen_Max];
%支路电流约束
C_fluc = [C_fluc, 0 <= I, I <= 11];
%区域FCE约束
%C_fluc = [C_fluc, P_FCE(69, :) <= RES_cspt_all / 2]
C_fluc = [C_fluc, UP_P_FCE <= RES_cspt / 4];
%% 10.波动模型设目标函数
f_loss_fluc = sum(sum(I .* (R * ones(1, T)))); %网损
f_RES_fluc = mean(mean((3 - UP_P_wt ./ Act_P_wt - UP_Q_wt ./ Act_Q_wt - UP_P_pv ./ Act_P_pv) / 3)); %新能源消纳量
f_sd_fluc = mean(sum(UP_P_load_sd)) / mean(sum(Act_P_Load(mpc.load_sd(:, 1), :))); %负荷削减量
f_ess_fluc = mean(sum(UP_P_ess_dch + UP_P_ess_ch)); %储能出力
f_fluc = eta_loss_fluc * f_loss_fluc + eta_RES_fluc * f_RES_fluc + eta_sd_fluc * f_sd_fluc + eta_ess_fluc * f_ess_fluc;
toc
%% 11.求解
%设求解器
ops = sdpsettings('solver', 'cplex');
sol = optimize(C_fluc, f_fluc, ops);
%分析错误
if sol.problem == 0
    disp('successful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end
%% 12.绘图
UP_P_wt_fluc = value(UP_P_wt);
UP_Q_wt_fluc = value(UP_Q_wt);
UP_P_pv_fluc = value(UP_P_pv);
UP_P_ess_dch_fluc = value(UP_P_ess_dch);
UP_P_ess_ch_fluc = value(UP_P_ess_ch);
SOC_ess_fluc = value(SOC_ess);
UP_P_load_sd_fluc = value(UP_P_load_sd);
P_Gen_fluc = value(P_Gen);
Q_Gen_fluc = value(Q_Gen);
UP_P_wt = value(UP_P_wt);

figure(1)
t = 1:1:24;
plot(t, UP_P_ess_dch_fluc(1, :) - UP_P_ess_ch_fluc(1, :), '--', 'linewidth', 2);
hold on
plot(t, UP_P_ess_dch_fluc(2, :) - UP_P_ess_ch_fluc(2, :), '-*', 'linewidth', 2);
xlabel('时刻/h');
ylabel('储能出力')
legend('ESS1', 'ESS2');

figure(2)
t = 1:1:24;
plot(t, SOC_ess_fluc(1, :), '--', 'linewidth', 1.8);
hold on
plot(t, SOC_ess_fluc(2, :), '-*', 'linewidth', 1.8);
xlabel('时刻/h');
ylabel('SOC');
legend('ESS1', 'ESS2');

figure(3)
t = 1:1:24;
plot(t, sum(Pre_P_wt), '-ko', 'linewidth', 1.6);
hold on
plot(t, sum(Act_P_wt), '--bo', 'linewidth', 1.6);
hold on
plot(t, sum(Pre_Q_wt), '-k+', 'linewidth', 1.6);
hold on
plot(t, sum(Act_Q_wt), '--b+', 'linewidth', 1.6);
hold on
plot(t, sum(Pre_P_pv), '-k*', 'linewidth', 1.6);
hold on
plot(t, sum(Act_P_pv), '--b*', 'linewidth', 1.6);
xlabel('时刻/h');
ylabel('风力与光照强度')
legend('预测风力P', '实际风力P', '预测风力Q', '实际风力Q', '预测光照P', '实际光照P')

figure(4)
t = 1:1:24;
plot(t, UP_P_wt_plan(1, :), '--k*', 'linewidth', 1.6);
hold on
plot(t, UP_P_wt(1, :), '--r*', 'linewidth', 1.6);
hold on
plot(t, Pre_P_wt(1, :), '-k', 'linewidth', 1.6);
hold on
plot(t, Act_P_wt(1, :), '-r', 'linewidth', 1.6);
xlabel('时刻/h');
ylabel('风电1');
legend('计划有功', '实际有功', '预测风力', '实际风力');

figure(5)
t = 1:1:24;
plot(t, UP_P_wt_plan(2, :), '--k*', 'linewidth', 1.6);
hold on
plot(t, UP_P_wt(2, :), '--r*', 'linewidth', 1.6);
hold on
plot(t, Pre_P_wt(2, :), '-k', 'linewidth', 1.6);
hold on
plot(t, Act_P_wt(2, :), '-r', 'linewidth', 1.6);
xlabel('时刻/h');
ylabel('风电2');
legend('计划有功', '实际有功', '预测风力', '实际风力');

figure(6)
t = 1:1:24;
plot(t, UP_P_wt_plan(3, :), '--k*', 'linewidth', 1.6);
hold on
plot(t, UP_P_wt(3, :), '--r*', 'linewidth', 1.6);
hold on
plot(t, Pre_P_wt(3, :), '-k', 'linewidth', 1.6);
hold on
plot(t, Act_P_wt(3, :), '-r', 'linewidth', 1.6);
xlabel('时刻/h');
ylabel('风电3');
legend('计划有功', '实际有功', '预测风力', '实际风力');

figure(7)
t = 1:1:24;
plot(t, UP_P_wt_plan(4, :), '--k*', 'linewidth', 1.6);
hold on
plot(t, UP_P_wt(4, :), '--r*', 'linewidth', 1.6);
hold on
plot(t, Pre_P_wt(4, :), '-k', 'linewidth', 1.6);
hold on
plot(t, Act_P_wt(4, :), '-r', 'linewidth', 1.6);
xlabel('时刻/h');
ylabel('风电4');
legend('计划有功', '实际有功', '预测风力', '实际风力');

figure(8)
t = 1:1:24;
plot(t, UP_P_wt_plan(5, :), '--k*', 'linewidth', 1.6);
hold on
plot(t, UP_P_wt(5, :), '--r*', 'linewidth', 1.6);
hold on
plot(t, Pre_P_wt(5, :), '-k', 'linewidth', 1.6);
hold on
plot(t, Act_P_wt(5, :), '-r', 'linewidth', 1.6);
xlabel('时刻/h');
ylabel('风电5');
legend('计划有功', '实际有功', '预测风力', '实际风力');
