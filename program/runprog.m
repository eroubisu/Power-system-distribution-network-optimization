function result = runprog(prog)
    disp(['mode:', int2str(prog.mode)]);
    tic
    run.fig6 = 0;
    run.fig5 = 0;
    if prog.mode == 1
        if prog.sit(1) == 1
            plan_data = opt_plan(fixed_data, Pre_data);
            result.sit1 = plan_data;
            toc
            disp(['optimized solution found', newline, '++++++++++++++++++++++++++++++++']);
        end
        for i = 2:5
            sit = ['sit', int2str(i)];
            if prog.sit(i) == 1
                if i == 5
                    run.fig5 = 1;
                end
                run.FCE = i;
                fluc_data = opt_fluc(fixed_data, Pre_data, Act_data, plan_data, run);
                result.(sit) = fluc_data;
                toc
                disp(['situation ', int2str(i), ' solved']);
                if prog.fig(i) == 1
                    run.locate = ['fig', int2str(i)];
                    draw(plan_data, fluc_data, run); % 画图
                    disp(['jpg saved to file fig', int2str(i), newline, '++++++++++++++++++++++++++++++++']);
                else
                    disp([newline, '++++++++++++++++++++++++++++++++']);
                end
            else
                toc
                disp(['situation ', int2str(i), ' skiped', newline, '++++++++++++++++++++++++++++++++']);
            end
        end
        if prog.fig(6) == 1
            run.fig6 = 1;
            draw(result.sit4, result.sit5, run);
        end
    elseif prog.mode == 0
        result = prog.result;
        for i = 2:5
            sit = ['sit', int2str(i)];
            if prog.fig(i) == 1
                if i == 5
                    run.fig5 = 1;
                end
                run.locate = ['fig', int2str(i)];
                draw(prog.result.sit1, prog.result.(sit), run); % 画图
                disp(['jpg saved to file fig', int2str(i), newline, '++++++++++++++++++++++++++++++++']);
            else
                disp(['situation ', int2str(i), ' skiped', newline, '++++++++++++++++++++++++++++++++']);
            end
            toc
        end
        if prog.fig(6) == 1
            run.fig6 = 1;
            draw(prog.result.sit4, prog.result.sit5, run);
        end
    end
end
