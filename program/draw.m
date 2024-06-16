function draw(mpc3, mpc4, run)
    if run.fig6 == 0
        mpc1 = Pre_data;
        mpc2 = Act_data;

        for i = 1:5
            figure('visible', 'off');
            t = 1:1:24;
            plot(t, mpc1.P_WT(i, :), '-r', 'linewidth', 1.6);
            hold on
            plot(t, mpc3.UP_P_wt(i, :), '--k*', 'linewidth', 1.6);
            xlabel('时刻/h');
            ylabel(['风电', int2str(i), '日前预测风力与计划出力']);
            legend('预测风力', '计划出力');
            saveas(gcf, [run.locate, '/风电', int2str(i), '日前预测风力与计划出力.jpg']);

            figure('visible', 'off');
            t = 1:1:24;
            plot(t, mpc2.P_WT(i, :), '-r', 'linewidth', 1.6);
            hold on
            plot(t, mpc4.UP_P_wt(i, :), '--k*', 'linewidth', 1.6);
            xlabel('时刻/h');
            ylabel(['风电', int2str(i), '日前实际风力与实际出力']);
            legend('实际风力', '实际出力');
            saveas(gcf, [run.locate, '/风电', int2str(i), '日内实际风力与实际出力.jpg']);
        end

        for i = 1:5
            figure('visible', 'off');
            t = 1:1:24;
            plot(t, mpc1.P_PV(i, :), '-r', 'linewidth', 1.6);
            hold on
            plot(t, mpc3.UP_P_pv(i, :), '--k*', 'linewidth', 1.6);
            xlabel('时刻/h');
            ylabel(['光伏', int2str(i), '日前预测光照强度与计划光伏发电']);
            legend('预测光照强度', '计划光伏发电');
            saveas(gcf, [run.locate, '/光伏', int2str(i), '日前预测光照强度与计划光伏发电.jpg']);

            figure('visible', 'off');
            t = 1:1:24;
            plot(t, mpc2.P_PV(i, :), '-r', 'linewidth', 1.6);
            hold on
            plot(t, mpc4.UP_P_pv(i, :), '--k*', 'linewidth', 1.6);
            xlabel('时刻/h');
            ylabel(['光伏', int2str(i), '日内实际光照强度与实际光伏发电']);
            legend('实际光照强度', '实际光伏发电');
            saveas(gcf, [run.locate, '/光伏', int2str(i), '日内实际光照强度与实际光伏发电.jpg']);
        end

        figure('visible', 'off');
        for i = 1:5
            if i == 1
                subplot(4, 2, 1);
            elseif i == 2
                t = title('AreaA储能出力');
                t.FontSize = 8;
                subplot(4, 2, 3);
            elseif i == 3
                subplot(4, 2, [2, 4]);
            elseif i == 4
                legend([h{1}; h{2}], '计划出力', '实际出力', 'Orientation', 'horizontal', 'location', 'NorthOutside');
                subplot(4, 2, [5, 7]);
            elseif i == 5
                subplot(4, 2, [6, 8]);
            end
            if i <= 2
                x = 'A';
            elseif i == 3
                x = 'B';
            elseif i == 4
                x = 'C';
            elseif i == 5
                x = 'D';
            end
            t = 1:1:24;
            h{1} = plot(t, mpc3.UP_P_ess_ch(i, :) - mpc3.UP_P_ess_dch(i, :), '--k*', 'linewidth', 1.6);
            hold on
            h{2} = plot(t, mpc4.UP_P_ess_ch(i, :) - mpc4.UP_P_ess_dch(i, :), '--r', 'linewidth', 1.6);
            %xlabel('时刻/h');
            ylabel(['储能', int2str(i), '出力']);
            if i >= 3
                t = title(['Area', x, '储能出力']);
                t.FontSize = 8;
            end
        end
        saveas(gcf, [run.locate, '/储能出力.jpg']);

        figure('visible', 'off');
        for i = 1:5
            if i == 1
                subplot(4, 2, 1);
            elseif i == 2
                t = title('AreaA储能SOC');
                t.FontSize = 8;
                subplot(4, 2, 3);
            elseif i == 3
                subplot(4, 2, [2, 4]);
            elseif i == 4
                legend([h{1}; h{2}], '计划SOC', '实际SOC', 'Orientation', 'horizontal', 'location', 'NorthOutside');
                subplot(4, 2, [5, 7]);
            elseif i == 5
                subplot(4, 2, [6, 8]);
            end
            if i <= 2
                x = 'A';
            elseif i == 3
                x = 'B';
            elseif i == 4
                x = 'C';
            elseif i == 5
                x = 'D';
            end
            t = 1:1:24;
            h{1} = plot(t, mpc3.SOC_ess(i, :), '-r', 'linewidth', 1.6);
            hold on
            h{2} = plot(t, mpc4.SOC_ess(i, :), '--b', 'linewidth', 1.6);
            %xlabel('时刻/h');
            ylabel(['储能', int2str(i), 'SOC']);
            if i >= 3
                t = title(['Area', x, '储能SOC']);
                t.FontSize = 8;
            end
        end
        saveas(gcf, [run.locate, '/储能SOC.jpg']);

        figure('visible', 'off');
        for i = 1:4
            subplot(2, 2, i)
            t = 1:1:24;
            if run.fig5 == 1
                h{1} = plot(t, repmat(0.0106, 1, 24), '-r', 'linewidth', 1.6);
            elseif run.fig5 == 0
                h{1} = plot(t, mpc4.RES_eta(i, :) / 8, '-r', 'linewidth', 1.6);
            end
            hold on
            h{2} = plot(t, mpc4.UP_P_FCE(i, :), '--b', 'linewidth', 1.6);
            xlabel('时刻/h');
            if i == 1
                x = 'A';
            elseif i == 2
                x = 'B';
                legend([h{1}; h{2}], '难度系数', '波动值', 'Orientation', 'horizontal', 'location', 'NorthOutside');
            elseif i == 3
                x = 'C';
            elseif i == 4
                x = 'D';
            end
            ylabel(['区域', x, '关口功率波动值']);
        end
        saveas(gcf, [run.locate, '/区域关口功率波动值.jpg']);

        figure('visible', 'off');
        for i = 1:5
            if i == 1
                subplot(4, 2, 1);
            elseif i == 2
                title('AreaA负荷削减')
                subplot(4, 2, 3);
            elseif i == 3
                subplot(4, 2, [2, 4]);
            elseif i == 4
                legend('负荷削减量');
                subplot(4, 2, [5, 7]);
            elseif i == 5
                subplot(4, 2, [6, 8]);
            end
            if i <= 2
                x = 'A';
            elseif i == 3
                x = 'B';
            elseif i == 4
                x = 'C';
            elseif i == 5
                x = 'D';
            end
            t = 1:1:24;
            plot(t, mpc4.UP_P_load_sd(i, :), '-k*', 'linewidth', 1.6);
            xlabel('时刻/h');
            ylabel(['节点', int2str(i), '负荷削减']);
            if i >= 3
                title(['Area', x, '负荷削减'])
            end
        end
        saveas(gcf, [run.locate, '/负荷削减.jpg']);
    elseif run.fig6 == 1
        s1.AreaA = mpc3.UP_P_wt(1, :) + sum(mpc3.UP_P_pv(2:5, :));
        s2.AreaA = mpc4.UP_P_wt(1, :) + sum(mpc4.UP_P_pv(2:5, :));
        s1.AreaB = mpc3.UP_P_wt(5, :);
        s2.AreaB = mpc4.UP_P_wt(5, :);
        s1.AreaC = mpc3.UP_P_pv(1, :);
        s2.AreaC = mpc4.UP_P_pv(1, :);
        s1.AreaD = sum(mpc3.UP_P_wt(2:4, :));
        s2.AreaD = sum(mpc4.UP_P_wt(2:4, :));
        x = [sum(s1.AreaA) sum(s1.AreaB) sum(s1.AreaC) sum(s1.AreaD);
             sum(s2.AreaA) sum(s2.AreaB) sum(s2.AreaC) sum(s2.AreaD)];
        for i = 1:4
            x(3, i) = min(x(1, i), x(2, i)) / max(x(1, i), x(2, i));
            x(4, i) = 1;
        end
        figure('visible', 'off');
        t = 1:1:4;
        line([-1, 6], [1, 1], 'Color', 'red', 'LineStyle', '--', 'linewidth', 1.6)
        hold on
        bar(t, x(3, :), 'w');
        plot(t, x(3, :), '-*', 'linewidth', 1.6)
        ylim([0, 1.2]);
        legend('情况1', '情况2占比');
        title(['新能源出力比较'])
        saveas(gcf, ['fig1', '/各区域新能源日内总出力比较.jpg']);
    end
end
