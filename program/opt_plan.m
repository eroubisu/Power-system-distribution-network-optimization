function data2 = opt_plan(fixed_data, Pre_data)
    mpc0 = fixed_data;
    mpc1 = Pre_data; % 预测数据

    %% 1.设参
    T = 24; %时段数为24小时
    %读取数据
    %24h有功/无功负荷
    P_Load = mpc1.P_Load;
    Q_Load = mpc1.Q_Load;
    %支路参数
    branch = mpc0.branch;
    branch(:, 3) = branch(:, 3) * mpc0.baseMVA / (mpc0.Vb ^ 2); %标幺化
    R = real(branch(:, 3)); %线路电阻
    X = imag(branch(:, 3)); %线路电抗
    %RES参数
    P_WT = mpc1.P_WT; %风电有功
    Q_WT = mpc1.Q_WT; %风电无功
    P_PV = mpc1.P_PV; %光伏有功
    %储能ESS
    Bus_ess = mpc0.ess(:, 1); %储能接入节点
    P_ess_ch_Max = repmat(mpc0.ess(:, 2), 1, T); %最大充电功率
    P_ess_dch_Max = repmat(mpc0.ess(:, 3), 1, T); %最大放电功率
    storage_ess = repmat(mpc0.ess(:, 4), 1, T); %各储能容量
    %机组数
    num_Bus = length(P_Load); %节点数
    num_Line = length(branch); %支路数
    num_wt = size(P_WT, 1); %风机数
    num_pv = size(P_PV, 1); %光伏数
    num_ess = size(Bus_ess, 1); %ESS数
    %目标函数系数
    eta_loss_plan = mpc0.eta_loss_plan;
    eta_RES_plan = mpc0.eta_RES_plan;
    eta_ess_plan = mpc0.eta_ess_plan;
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
    P_Gen_Max = [zeros(num_Bus - 1, T); 0.2 * ones(1, T)];
    P_Gen_Min = [zeros(num_Bus - 1, T); 0.1 * ones(1, T)];
    Q_Gen_Max = [zeros(num_Bus - 1, T); 0.2 * ones(1, T)];
    Q_Gen_Min = [zeros(num_Bus - 1, T); -0.2 * ones(1, T)];

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
    %节点功率
    P_bus_in = -upstream * P_Line + upstream * (I .* (R * ones(1, T))) + dnstream * P_Line; %节点注入有功
    Q_bus_in = -upstream * Q_Line + upstream * (I .* (X * ones(1, T))) + dnstream * Q_Line; %节点注入无功

    %% 3.设约束
    C = [];
    %% 3.1储能装置（ESS）约束
    %充放电状态约束
    C = [C, U_ess_dch(1:2, :) + U_ess_ch(1:2, :) <= 1];
    %充放电功率约束
    C = [C, 0 <= UP_P_ess_ch <= U_ess_ch .* P_ess_ch_Max];
    C = [C, 0 <= UP_P_ess_dch <= U_ess_dch .* P_ess_dch_Max];
    %容量约束
    for t = 1:23
        C = [C, SOC_ess(:, t + 1) == SOC_ess(:, t) + UP_P_ess_ch(:, t) - UP_P_ess_dch(:, t)];
    end
    C = [C, 0.1 * storage_ess <= SOC_ess <= 0.9 * storage_ess];
    C = [C, SOC_ess(:, 1) == SOC_ess(:, 24)];
    %投入节点选择
    for i = 1:num_ess
        P_ess_dch = [P_ess_dch(1:mpc0.ess(i, 1) - 1, :); UP_P_ess_dch(i, :); P_ess_dch(mpc0.ess(i, 1):end, :)];
        P_ess_ch = [P_ess_ch(1:mpc0.ess(i, 1) - 1, :); UP_P_ess_ch(i, :); P_ess_ch(mpc0.ess(i, 1):end, :)];
    end
    %% 3.2风机+光伏约束
    %出力功率约束
    C = [C, 0 <= UP_P_wt, UP_P_wt <= P_WT];
    C = [C, 0 <= UP_Q_wt, UP_Q_wt <= Q_WT];
    C = [C, 0 <= UP_P_pv, UP_P_pv <= P_PV];
    %投入节点选择
    for i = 1:num_wt
        P_wt = [P_wt(1:mpc0.WT(i) - 1, :); UP_P_wt(i, :); P_wt(mpc0.WT(i):end, :)];
        Q_wt = [Q_wt(1:mpc0.WT(i) - 1, :); UP_Q_wt(i, :); Q_wt(mpc0.WT(i):end, :)];
    end
    for i = 1:num_pv
        P_pv = [P_pv(1:mpc0.PV(i) - 1, :); UP_P_pv(i, :); P_pv(mpc0.PV(i):end, :)];
    end
    %% 3.3潮流约束
    %功率平衡约束
    C = [C, P_bus_in + P_Load - P_Gen - P_wt - P_pv - P_ess_dch + P_ess_ch == 0];
    C = [C, Q_bus_in + Q_Load - Q_wt - Q_Gen == 0];
    %欧姆定律约束
    C = [C, V(branch(:, 2), :) == V(branch(:, 1), :) - 2 * (R * ones(1, 24)) .* P_Line - 2 * (X * ones(1, 24)) .* Q_Line + ((R .^ 2 + X .^ 2) * ones(1, 24)) .* I];
    %二阶锥约束
    C = [C, V(branch(:, 1), :) .* I >= P_Line .^ 2 + Q_Line .^ 2];
    %% 3.4通用约束
    %节点电压约束
    C = [C, V_Bus_Min <= V <= V_Bus_Max];
    %发电机功率约束
    C = [C, P_Gen_Min <= P_Gen <= P_Gen_Max, Q_Gen_Min <= Q_Gen <= Q_Gen_Max];
    %支路电流约束
    C = [C, 0 <= I, I <= 11];

    %% 4.设目标函数
    f_loss_plan = sum(sum(I .* (R * ones(1, T)))); %网损
    f_RES_plan = mean(mean((3 - UP_P_wt ./ P_WT - UP_Q_wt ./ Q_WT - UP_P_pv ./ P_PV) / 3)); %新能源消纳量
    f_ess_plan = mean(sum(UP_P_ess_dch + UP_P_ess_ch)); %储能出力
    f = eta_loss_plan * f_loss_plan + eta_RES_plan * f_RES_plan + eta_ess_plan * f_ess_plan;
    %f = eta_RES_plan * f_RES_plan + eta_ess_plan * f_ess_plan;

    %% 5.求解
    %设求解器
    ops = sdpsettings('solver', 'cplex', 'verbose', 0);
    sol = optimize(C, f, ops);
    %分析错误
    if sol.problem
        disp('error');
        yalmiperror(sol.problem)
    end

    %% 6.传出参数
    data2 = struct;
    data2.UP_P_wt = value(UP_P_wt);
    data2.UP_Q_wt = value(UP_Q_wt);
    data2.UP_P_pv = value(UP_P_pv);
    data2.P_wt = value(P_wt);
    data2.Q_wt = value(Q_wt);
    data2.P_pv = value(P_pv);
    data2.UP_P_ess_dch = value(UP_P_ess_dch);
    data2.UP_P_ess_ch = value(UP_P_ess_ch);
    data2.SOC_ess = value(SOC_ess);
    data2.P_Gen = value(P_Gen);
    data2.Q_Gen = value(Q_Gen);
    data2.P_Line = value(P_Line);
    data2.f_RES_plan = value(f_RES_plan);
    data2.f_ess_plan = value(f_ess_plan);
end