function draw(plan_data, fluc_data)
    mpc1 = Pre_data;
    mpc2 = Act_data;
    mpc3 = plan_data;
    mpc4 = fluc_data;

    figure(3)
    t = 1:1:24;
    plot(t, sum(mpc1.P_WT(:, 2:end)), '-ko', 'linewidth', 1.6);
    hold on
    plot(t, sum(mpc2.P_WT(:, 2:end)), '--bo', 'linewidth', 1.6);
    hold on
    plot(t, sum(mpc1.Q_WT(:, 2:end)), '-k+', 'linewidth', 1.6);
    hold on
    plot(t, sum(mpc2.Q_WT(:, 2:end)), '--b+', 'linewidth', 1.6);
    hold on
    plot(t, sum(mpc1.P_PV(:, 2:end)), '-k*', 'linewidth', 1.6);
    hold on
    plot(t, sum(mpc2.P_PV(:, 2:end)), '--b*', 'linewidth', 1.6);
    xlabel('时刻/h');
    ylabel('风力与光照强度')
    legend('预测风力P', '实际风力P', '预测风力Q', '实际风力Q', '预测光照P', '实际光照P')

    figure(4)
    t = 1:1:24;
    plot(t, mpc3.UP_P_wt(1, :), '--k*', 'linewidth', 1.6);
    hold on
    plot(t, mpc4.UP_P_wt(1, :), '--r*', 'linewidth', 1.6);
    hold on
    plot(t, mpc1.P_WT(1, 2:end), '-k', 'linewidth', 1.6);
    hold on
    plot(t, mpc2.P_WT(1, 2:end), '-r', 'linewidth', 1.6);
    xlabel('时刻/h');
    ylabel('风电1');
    legend('计划有功', '实际有功', '预测风力', '实际风力');

    figure(5)
    t = 1:1:24;
    plot(t, mpc3.UP_P_wt(2, :), '--k*', 'linewidth', 1.6);
    hold on
    plot(t, mpc4.UP_P_wt(2, :), '--r*', 'linewidth', 1.6);
    hold on
    plot(t, mpc1.P_WT(2, 2:end), '-k', 'linewidth', 1.6);
    hold on
    plot(t, mpc2.P_WT(2, 2:end), '-r', 'linewidth', 1.6);
    xlabel('时刻/h');
    ylabel('风电2');
    legend('计划有功', '实际有功', '预测风力', '实际风力');

    figure(6)
    t = 1:1:24;
    plot(t, mpc3.UP_P_wt(3, :), '--k*', 'linewidth', 1.6);
    hold on
    plot(t, mpc4.UP_P_wt(3, :), '--r*', 'linewidth', 1.6);
    hold on
    plot(t, mpc1.P_WT(3, 2:end), '-k', 'linewidth', 1.6);
    hold on
    plot(t, mpc2.P_WT(3, 2:end), '-r', 'linewidth', 1.6);
    xlabel('时刻/h');
    ylabel('风电3');
    legend('计划有功', '实际有功', '预测风力', '实际风力');

    figure(7)
    t = 1:1:24;
    plot(t, mpc3.UP_P_wt(4, :), '--k*', 'linewidth', 1.6);
    hold on
    plot(t, mpc4.UP_P_wt(4, :), '--r*', 'linewidth', 1.6);
    hold on
    plot(t, mpc1.P_WT(4, 2:end), '-k', 'linewidth', 1.6);
    hold on
    plot(t, mpc2.P_WT(4, 2:end), '-r', 'linewidth', 1.6);
    xlabel('时刻/h');
    ylabel('风电4');
    legend('计划有功', '实际有功', '预测风力', '实际风力');

    figure(8)
    t = 1:1:24;
    plot(t, mpc3.UP_P_wt(5, :), '--k*', 'linewidth', 1.6);
    hold on
    plot(t, mpc4.UP_P_wt(5, :), '--r*', 'linewidth', 1.6);
    hold on
    plot(t, mpc1.P_WT(5, 2:end), '-k', 'linewidth', 1.6);
    hold on
    plot(t, mpc2.P_WT(5, 2:end), '-r', 'linewidth', 1.6);
    xlabel('时刻/h');
    ylabel('风电5');
    legend('计划有功', '实际有功', '预测风力', '实际风力');
end
