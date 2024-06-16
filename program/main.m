prog.fig = [1 1 1 1 1 1]; % 0：不画图；1：画图
prog.sit = [1 1 1 1 1]; % 日前计划 & 场景1 2 3 & 控制变量对照组
prog.mode = 1; %1：求解+画图；0：只画图，不求解（需要先运行一遍1，别clear）
prog.result = runprog(prog);
