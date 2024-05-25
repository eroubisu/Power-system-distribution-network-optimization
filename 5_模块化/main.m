clear
clc
tic
warning off
plan_data = opt_plan(Pre_data); %计划发电
fluc_data = opt_fluc(Act_data, plan_data); %考虑负荷波动的分区协同平衡
draw(plan_data, fluc_data); %绘图
toc
