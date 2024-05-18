UP_P_wt_plan = value(UP_P_wt);
UP_Q_wt_plan = value(UP_Q_wt);
UP_P_pv_plan = value(UP_P_pv);
UP_P_ess_dch_plan = value(UP_P_ess_dch);
UP_P_ess_ch_plan = value(UP_P_ess_ch);
SOC_ess_plan = value(SOC_ess);
P_Gen_plan = value(P_Gen);
Q_Gen_plan = value(Q_Gen);

figure(1)
t = 1:1:24;
plot(t, sum(UP_P_wt_plan), '--k', 'linewidth', 2);
hold on
plot(t, sum(UP_Q_wt_plan), '-*k', 'linewidth', 2);
hold on
plot(t, sum(UP_P_pv_plan), '--r', 'linewidth', 2);
xlabel('时刻/h');
ylabel('新能源出力')
legend('风电有功', '风电无功', '光伏有功');

figure(2)
t = 1:1:24;
plot(t, UP_P_ess_dch_plan(1, :) - UP_P_ess_ch(1, :), '--', 'linewidth', 2);
hold on
plot(t, UP_P_ess_dch_plan(2, :) - UP_P_ess_ch(2, :), '-*', 'linewidth', 2);
xlabel('时刻/h');
ylabel('储能出力')
legend('ESS1', 'ESS2');

figure(3)
t = 1:1:24;
plot(t, SOC_ess_plan(1, :), '--', 'linewidth', 1.8);
hold on
plot(t, SOC_ess_plan(2, :), '-*', 'linewidth', 1.8);
xlabel('时刻/h');
ylabel('SOC');
legend('ESS1', 'ESS2');

figure(4)
t = 1:1:24;
plot(t, P_Gen_plan(69, :), '--', 'linewidth', 1.6);
hold on
plot(t, Q_Gen_plan(69, :), '-*', 'linewidth', 1.6);
xlabel('时刻/h');
ylabel('机组出力');
legend('Pg', 'Qg');
