function result = opt_fluc(fixed_data, Pre_data, Act_data, plan_data, run)
    mpc0 = fixed_data;
    mpc1 = Pre_data; % 预测数据
    mpc2 = Act_data; % 实际数据
    %% 1.设参
    T = 24;
    %读取数据
    %24h有功/无功负荷
    P_Load = mpc2.P_Load;
    Q_Load = mpc2.Q_Load;
    %支路参数
    branch = mpc0.branch;
    branch(:, 3) = branch(:, 3) * mpc0.baseMVA / (mpc0.Vb ^ 2); %标幺化
    R = real(branch(:, 3)); %线路电阻
    X = imag(branch(:, 3)); %线路电抗
    %RES参数
    P_WT = mpc2.P_WT; %风电有功
    Q_WT = mpc2.Q_WT; %风电无功
    P_PV = mpc2.P_PV; %光伏有功
    %储能ESS
    Bus_ess = mpc0.ess(:, 1); %储能接入节点
    P_ess_ch_Max = repmat(mpc0.ess(:, 2), 1, T); %最大充电功率
    P_ess_dch_Max = repmat(mpc0.ess(:, 3), 1, T); %最大放电功率
    storage_ess = repmat(mpc0.ess(:, 4), 1, T); %各储能容量
    %负荷响应
    Bus_load_sd = mpc0.load_sd;
    %机组数
    num_Bus = length(P_Load); %节点数
    num_Line = length(branch); %支路数
    num_wt = size(P_WT, 1); %风机数
    num_pv = size(P_PV, 1); %光伏数
    num_ess = size(Bus_ess, 1); %ESS数
    num_load_sd = size(Bus_load_sd, 1); %支持负荷削减的节点数
    %目标函数系数
    eta_loss_fluc = mpc0.eta_loss_fluc;
    eta_RES_fluc = mpc0.eta_RES_fluc;
    eta_sd_fluc = mpc0.eta_sd_fluc;
    eta_ess_fluc = mpc0.eta_ess_fluc;
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
    % 新能源消纳
    areaA = mpc0.areaA;
    areaB = mpc0.areaB;
    areaC = mpc0.areaC;
    areaD = mpc0.areaD;
    RES_cspt_all = mean(sum(plan_data.UP_P_wt + plan_data.UP_Q_wt + plan_data.UP_P_pv)); %计划新能源消纳量
    RES_cspt_A = sum(mean(plan_data.P_wt(areaA, :) + plan_data.Q_wt(areaA, :) + plan_data.P_pv(areaA, :)));
    RES_cspt_B = sum(mean(plan_data.P_wt(areaB, :) + plan_data.Q_wt(areaB, :) + plan_data.P_pv(areaB, :)));
    RES_cspt_C = sum(mean(plan_data.P_wt(areaC, :) + plan_data.Q_wt(areaC, :) + plan_data.P_pv(areaC, :)));
    RES_cspt_D = sum(mean(plan_data.P_wt(areaD, :) + plan_data.Q_wt(areaD, :) + plan_data.P_pv(areaD, :)));
    RES_cspt = repmat([RES_cspt_A; RES_cspt_B; RES_cspt_C; RES_cspt_D], 1, T);
    RES_dmdA = [];
    RES_dmdB = [];
    RES_dmdC = [];
    RES_dmdD = [];
    for i = 1:num_wt
        if ismember(mpc0.WT(i), areaA)
            RES_dmdA = [RES_dmdA; mpc1.P_WT(i, :) - mpc2.P_WT(i, :)];
            RES_dmdA = [RES_dmdA; mpc1.Q_WT(i, :) - mpc2.Q_WT(i, :)];
        end
        if ismember(mpc0.WT(i), areaB)
            RES_dmdB = [RES_dmdB; mpc1.P_WT(i, :) - mpc2.P_WT(i, :)];
            RES_dmdB = [RES_dmdA; mpc1.Q_WT(i, :) - mpc2.Q_WT(i, :)];
        end
        if ismember(mpc0.WT(i), areaC)
            RES_dmdC = [RES_dmdC; mpc1.P_WT(i, :) - mpc2.P_WT(i, :)];
            RES_dmdC = [RES_dmdA; mpc1.Q_WT(i, :) - mpc2.Q_WT(i, :)];
        end
        if ismember(mpc0.WT(i), areaD)
            RES_dmdD = [RES_dmdD; mpc1.P_WT(i, :) - mpc2.P_WT(i, :)];
            RES_dmdD = [RES_dmdA; mpc1.Q_WT(i, :) - mpc2.Q_WT(i, :)];
        end
    end
    for i = 1:num_pv
        if ismember(mpc0.PV(i), areaA)
            RES_dmdA = [RES_dmdA; mpc1.P_PV(i, :) - mpc2.P_PV(i, :)];
        end
        if ismember(mpc0.PV(i), areaB)
            RES_dmdB = [RES_dmdB; mpc1.P_PV(i, :) - mpc2.P_PV(i, :)];
        end
        if ismember(mpc0.PV(i), areaC)
            RES_dmdC = [RES_dmdC; mpc1.P_PV(i, :) - mpc2.P_PV(i, :)];
        end
        if ismember(mpc0.PV(i), areaD)
            RES_dmdD = [RES_dmdD; mpc1.P_PV(i, :) - mpc2.P_PV(i, :)];
        end
    end
    RES_dmdA = sum(mean(RES_dmdA) .* mean(RES_dmdA));
    RES_dmdB = sum(mean(RES_dmdB) .* mean(RES_dmdB));
    RES_dmdC = sum(mean(RES_dmdC) .* mean(RES_dmdC));
    RES_dmdD = sum(mean(RES_dmdD) .* mean(RES_dmdD));
    RES_dmd = repmat([RES_dmdA; RES_dmdB; RES_dmdC; RES_dmdD], 1, T);
    RES_eta = RES_cspt + RES_dmd * 12;

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
    %FCE，其中UP_P_FCE表示各区域跟踪关口功率的不平衡量
    P_FCE = P_Gen - plan_data.P_Gen;
    UP_P_FCE = [P_Gen(69, :) - plan_data.P_Gen(69, :) - (P_Line(27, :) - plan_data.P_Line(27, :)) - (P_Line(36, :) - plan_data.P_Line(36, :)); P_Line(27, :) - plan_data.P_Line(27, :); P_Line(36, :) - plan_data.P_Line(36, :); P_Line(10, :) - plan_data.P_Line(10, :)];
    UP_P_FCE = abs(UP_P_FCE);
    UP_P_FCE(:, [12 21]) = 0; %通过运行程序，找到最容易不平衡的两个节点，视为特殊情况，不将其纳入考虑（直接令不平衡量为0）
    %% 9.波动模型设约束
    C = [];
    %% 9.1需求响应约束
    C = [C, 0 <= UP_P_load_sd <= 0.8 * P_Load(mpc0.load_sd(:, 1), :)];
    %C = [C, -0.8 * P_Load(mpc0.load_sd(:, 1), :) <= UP_P_load_sd, UP_P_load_sd <= 0.8 * P_Load(mpc0.load_sd(:, 1), :)];
    for i = 1:num_load_sd
        P_load_sd = [P_load_sd(1:mpc0.load_sd(i, 1) - 1, :); UP_P_load_sd(i, :); P_load_sd(mpc0.load_sd(i, 1):end, :)];
    end
    %% 9.2储能装置（ESS）约束
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
    %% 9.3风电+光伏约束
    C = [C, 0 <= UP_P_wt, UP_P_wt <= P_WT];
    C = [C, 0 <= UP_Q_wt, UP_Q_wt <= Q_WT];
    for i = 1:num_wt
        P_wt = [P_wt(1:mpc0.WT(i) - 1, :); UP_P_wt(i, :); P_wt(mpc0.WT(i):end, :)];
        Q_wt = [Q_wt(1:mpc0.WT(i) - 1, :); UP_Q_wt(i, :); Q_wt(mpc0.WT(i):end, :)];
    end
    C = [C, 0 <= UP_P_pv, UP_P_pv <= P_PV];
    for i = 1:num_pv
        P_pv = [P_pv(1:mpc0.PV(i) - 1, :); UP_P_pv(i, :); P_pv(mpc0.PV(i):end, :)];
    end
    C = [C, mean(sum(UP_P_wt + UP_Q_wt + UP_P_pv)) >= min(mean(sum(P_WT + Q_WT + P_PV)), RES_cspt_all) * 0.6];
    %% 9.4潮流约束
    P_bus_in = -upstream * P_Line + upstream * (I .* (R * ones(1, T))) + dnstream * P_Line; %节点注入有功
    Q_bus_in = -upstream * Q_Line + upstream * (I .* (X * ones(1, T))) + dnstream * Q_Line; %节点注入无功
    C = [C, P_bus_in + P_Load - P_load_sd - P_Gen - P_wt - P_pv - P_ess_dch + P_ess_ch == 0];
    C = [C, Q_bus_in + Q_Load - Q_wt - Q_Gen == 0];
    %欧姆定律约束
    C = [C, V(branch(:, 2), :) == V(branch(:, 1), :) - 2 * (R * ones(1, 24)) .* P_Line - 2 * (X * ones(1, 24)) .* Q_Line + ((R .^ 2 + X .^ 2) * ones(1, 24)) .* I];
    %二阶锥约束
    C = [C, V(branch(:, 1), :) .* I >= P_Line .^ 2 + Q_Line .^ 2];
    %% 9.5通用约束
    %节点电压约束
    C = [C, V_Bus_Min <= V, V <= V_Bus_Max];
    %发电机功率约束
    C = [C, P_Gen_Min <= P_Gen, P_Gen <= P_Gen_Max, Q_Gen_Min <= Q_Gen, Q_Gen <= Q_Gen_Max];
    %支路电流约束
    C = [C, 0 <= I, I <= 11];
    %区域FCE约束
    if run.FCE == 2
        C = C;
    elseif run.FCE == 3
        C = [C, UP_P_FCE == 0];
    elseif run.FCE == 4
        C = [C, UP_P_FCE <= RES_eta / 8];
    elseif run.FCE == 5
        %C = [C, UP_P_FCE <= 0.0106];
        C = [C, UP_P_FCE <= 0.0106];
    end
    %% 10.波动模型设目标函数
    f_loss_fluc = sum(sum(I .* (R * ones(1, T)))); %网损
    f_RES_fluc = mean(mean((3 - UP_P_wt ./ P_WT - UP_Q_wt ./ Q_WT - UP_P_pv ./ P_PV) / 3)); %新能源消纳量
    f_sd_fluc = mean(sum(UP_P_load_sd)) / mean(sum(P_Load(mpc0.load_sd(:, 1), :))); %负荷削减量
    f_ess_fluc = mean(sum(UP_P_ess_dch + UP_P_ess_ch)); %储能出力
    f_fluc = eta_loss_fluc * f_loss_fluc + eta_RES_fluc * f_RES_fluc + eta_sd_fluc * f_sd_fluc + eta_ess_fluc * f_ess_fluc;

    %% 11.求解
    %设求解器
    ops = sdpsettings('solver', 'cplex', 'verbose', 0);
    sol = optimize(C, f_fluc, ops);
    %分析错误
    if sol.problem
        disp('error');
        yalmiperror(sol.problem)
    end
    %% 12.传出参数
    result = struct;
    result.UP_P_wt = value(UP_P_wt);
    result.UP_Q_wt = value(UP_Q_wt);
    result.UP_P_pv = value(UP_P_pv);
    result.P_wt = value(P_wt);
    result.Q_wt = value(Q_wt);
    result.P_pv = value(P_pv);
    result.UP_P_ess_dch = value(UP_P_ess_dch);
    result.UP_P_ess_ch = value(UP_P_ess_ch);
    result.SOC_ess = value(SOC_ess);
    result.P_Gen = value(P_Gen);
    result.Q_Gen = value(Q_Gen);
    result.P_Line = value(P_Line);
    result.UP_P_load_sd = value(UP_P_load_sd);
    result.UP_P_FCE = value(UP_P_FCE);
    result.RES_cspt = value(RES_cspt);
    result.f_RES_fluc = value(f_RES_fluc);
    result.f_sd_fluc = value(f_sd_fluc);
    result.f_ess_fluc = value(f_ess_fluc);
    result.RES_dmdA = value(RES_dmdA);
    result.RES_dmdB = value(RES_dmdB);
    result.RES_dmdC = value(RES_dmdC);
    result.RES_dmdD = value(RES_dmdD);
    result.RES_eta = value(RES_eta);
    result.f_loss_fluc = value(f_loss_fluc);
end
