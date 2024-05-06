clear
clc
tic
warning off
%% 1.���
mpc = IEEE33BW;
pload = mpc.Pload * 100; %�ڵ��й�����,ԭ����100MVA����Ϊ1MVA
PPload = sum(mpc.Pload * 100);
qload = mpc.Qload * 100; %�ڵ��޹�����
branch = mpc.branch;
branch(:, 3) = branch(:, 3) * 1 / (10 ^ 2); %���迹����ֵ
R = real(branch(:, 3));
X = imag(branch(:, 3));
T = 24; %ʱ����Ϊ24Сʱ
nb = 33; %�ڵ���
nl = 32; %֧·��
nwt = 2; %���
npv = 1; %���
S_pv = 1.5;
RR = xlsread('����ǿ������.xlsx', 'B2:B25');
pp_pv = xlsread('����ǿ������.xlsx', 'D3:AA3');
C_e = xlsread('���.xlsx', 'A1:X1');
C_ae = xlsread('�������.xlsx', 'A1:X1');
C_de = xlsread('������Ӧ�����۸�.xlsx', 'A1:X1');
R_STD = 1000; %��׼̫������
R_C = 150; %�����
PP_pv = zeros(npv, T);
% %%%%�������ģ��
for i = 1:1:npv
    for t = 1:1:T
        if (0 <= RR(t)) & (RR(t) <= R_C)
            PP_pv(i, t) = S_pv * ((RR(t)) ^ 2) / (R_STD * R_C);
        else
            if (R_C <= RR(t)) & (RR(t) <= R_STD)
                PP_pv(i, t) = S_pv * RR(t) / R_STD;
            else
                if (R_STD <= RR(t))
                    PP_pv(i, t) = S_pv;
                end
            end
        end
    end
end
v = xlsread('��������.xlsx', 'B2:B25');
pp_wt = xlsread('��������.xlsx', 'G2:AD2');
v_in = 3; %�������
v_r = 11.3; %�����
v_out = 25; %�г�����
S_wt = [1; 1.5];
PP_wt = zeros(nwt, T);
%%%������ģ��
for i = 1:1:nwt
    for t = 1:1:T
        if (0 <= v(t)) & (v(t) <= v_in) | (v_out <= v(t))
            PP_wt(i, t) = 0;
        else
            if (v_in <= v(t)) & (v(t) < v_r)
                PP_wt(i, t) = [(v(t) - v_in) / (v_r - v_in)] * S_wt(i);
            else
                if (v_r <= v(t)) & (v(t) < v_out)
                    PP_wt(i, t) = S_wt(i);
                end
            end
        end
    end
end
ness = 2; %ESS��
v_in = 3; %�������
v_r = 11.3; %�����
v_out = 25; %�г�����
upstream = zeros(nb, nl);
dnstream = zeros(nb, nl);
for t = 1:nl
    upstream(t, t) = 1;
end
for t = [1:16, 18:20, 22:23, 25:31]
    dnstream(t, t + 1) = 1;
end
dnstream(1, 18) = 1;
dnstream(2, 22) = 1;
dnstream(5, 25) = 1;
dnstream(33, 1) = 1;
Vmax = [1.05 * 1.05 * ones(nb - 1, T)
        1.05 * 1.05 * ones(1, T)];
Vmin = [0.95 * 0.95 * ones(nb - 1, T)
        0.95 * 0.95 * ones(1, T)];
Pgmax = [zeros(nb - 1, T)
         3 * ones(1, T)];
Pgmin = [zeros(nb - 1, T)
         0 * ones(1, T)];
Qgmax = [zeros(nb - 1, T)
         3 * ones(1, T)];
Qgmin = [zeros(nb - 1, T)
         -1 * ones(1, T)];
%% 2.�����
V = sdpvar(nb, T); %��ѹ��ƽ��
I = sdpvar(nl, T); %������ƽ��
P = sdpvar(nl, T); %��·�й�
Q = sdpvar(nl, T); %��·�޹�
Pg = sdpvar(nb, T); %������й�
Qg = sdpvar(nb, T); %������޹�
p_wt = sdpvar(nwt, T); %����й�
p_pv = sdpvar(npv, T); %����й�
s_IL = binvar(1, T);
Temp_shift = binvar(3, T);
Lshift = sdpvar(3, T);
Lshift_old = xlsread('��ת�Ƹ���.xlsx', 'A1:X3');
PLOAD = sdpvar(nb, T);
Pload = sum(PLOAD);
p_dch = sdpvar(ness, T); %ESS�ŵ繦��
p_ch = sdpvar(ness, T); %ESS��繦��
u_dch = binvar(ness, T); %ESS�ŵ�״̬
u_ch = binvar(ness, T); %ESS���״̬
E_ess = sdpvar(ness, T); %ESS�ĵ���
S_IL1 = sdpvar(1, T, 'full');
S_IL2 = sdpvar(1, T, 'full');
%% 3.��Լ��
C = [];
%% ������ӦԼ��
C = [C, 0 <= S_IL1 <= 0.8 * pload(10, :)];
C = [C, 0 <= S_IL2 <= 0.8 * pload(26, :)];
IL1 = [zeros(9, T); S_IL1(1, :); zeros(23, T)]; %��������
IL2 = [zeros(25, T); S_IL2(1, :); zeros(7, T)];
C = [C, pload - IL1 - IL2 == PLOAD];
%% ����װ�ã�ESS��Լ��
%��ŵ�״̬Լ��
C = [C, u_dch(1, :) + u_ch(1, :) <= 1];
C = [C, u_dch(2, :) + u_ch(2, :) <= 1];
C = [C, u_dch(1, :) + u_dch(2, :) <= 1];
C = [C, u_ch(1, :) + u_ch(2, :) <= 1];
%��ŵ�ʱ��Լ��
C = [C, [sum(u_dch(1, 8:15)) == 0, sum(u_dch(1, 18:21)) == 0, sum(u_ch(1, 20:21)) == 1]];
C = [C, [sum(u_dch(2, 8:15)) == 0, sum(u_dch(2, 18:21)) == 0, sum(u_ch(2, 20:21)) == 1]];
%����Լ��
C = [C, 0 <= p_dch <= u_dch * 0.375];
C = [C, 0 <= p_ch <= u_ch * 0.375];
%����Լ��
for t = 1:23
    C = [C, E_ess(1, t + 1) == E_ess(1, t) + 0.9 * p_ch(1, t) - 0.9 * p_dch(1, t)];
end
for t = 1:23
    C = [C, E_ess(2, t + 1) == E_ess(2, t) + 0.9 * p_ch(2, t) - 0.9 * p_dch(2, t)];
end
C = [C, 0.1 <= E_ess(1, :) / 1.03 <= 0.9];
C = [C, 0.1 <= E_ess(2, :) / 1.21 <= 0.9];
%Ͷ��ڵ�ѡ��
P_dch1 = [zeros(16, T); p_dch(1, :); zeros(16, T)];
P_ch1 = [zeros(16, T); p_ch(1, :); zeros(16, T)];
P_dch2 = [zeros(16, T); p_dch(2, :); zeros(16, T)];
P_ch2 = [zeros(16, T); p_ch(2, :); zeros(16, T)];
%% ���+���Լ��
P_WT = xlsread('���ʵ�ʳ���', 'A1:X2');
P_PV = xlsread('���ʵ�ʳ���', 'A1:X1');
C = [C, 0 <= p_wt, p_wt <= 1.2 * P_WT];
P_wt = [zeros(15, T); p_wt(1, :); zeros(10, T); p_wt(2, :); zeros(6, T)];
C = [C, 0 <= p_pv, p_pv <= 1.2 * P_PV];
P_pv = [zeros(15, T); p_pv; zeros(17, T)];
%% ����Լ��
Pin = -upstream * P + upstream * (I .* (R * ones(1, T))) + dnstream * P; %�ڵ�ע���й�
Qin = -upstream * Q + upstream * (I .* (X * ones(1, T))) + dnstream * Q; %�ڵ�ע���޹�
C = [C, Pin + pload - IL1 - IL2 - Pg - P_wt - P_pv - P_dch1 + P_ch1 - P_dch2 + P_ch2 == 0];
C = [C, Qin + qload - Qg == 0];
%ŷķ����Լ��
C = [C, V(branch(:, 2), :) == V(branch(:, 1), :) - 2 * (R * ones(1, 24)) .* P - 2 * (X * ones(1, 24)) .* Q + ((R .^ 2 + X .^ 2) * ones(1, 24)) .* I];
%����׶Լ��
C = [C, V(branch(:, 1), :) .* I >= P .^ 2 + Q .^ 2];
%% ͨ��Լ��
%�ڵ��ѹԼ��
C = [C, Vmin <= V, V <= Vmax];
%���������Լ��
C = [C, Pgmin <= Pg, Pg <= Pgmax, Qgmin <= Qg, Qg <= Qgmax];
%֧·����Լ��
C = [C, 0 <= I, I <= 12];
%% 4.��Ŀ�꺯��
Pload2 = PPload + sum(p_dch - p_ch);
Pload3 = Pload + sum(p_dch - p_ch);
%%%%%%%%%
for t = 1:1:23
    f1_0 = 0.25 * [(max(PPload) - min(PPload)) / mean(PPload(:))] + 0.75 * [(max(abs(PPload(t) - Pload(t + 1)))) / mean(PPload(:))];
    f1_1 = 0.25 * [(max(Pload) - min(Pload)) / mean(Pload(:))] + 0.75 * [(max(abs(Pload(t) - Pload(t + 1)))) / mean(Pload(:))];
    f1_2 = 0.25 * [(max(Pload2) - min(Pload2)) / mean(Pload2)] + 0.75 * [(max(abs(Pload2(t) - Pload2(t + 1)))) / mean(Pload2(:))];
    f1_3 = 0.25 * [(max(Pload3) - min(Pload3)) / mean(Pload3)] + 0.75 * [(max(abs(Pload3(t) - Pload3(t + 1)))) / mean(Pload3(:))];
end
%%%%%%%��ѹ����С%%%%%%%%%
f1 = 0;
for t = 1:1:T
    for n = 1:1:nb - 1
        f1 = f1 + abs(V(n, t) - V(n + 1, t));
    end
end
%%%%%%%%%
B1 = sum(sum(p_dch - p_ch) .* C_ae);
B2 = sum(sum(p_dch)) * 200;
R_ess = B1 + B2; %���ܲ��������� % % % % %
C_LA = sum([0.062, 0.043, 0.048] * Temp_shift .* C_de) + sum((S_IL1 + S_IL2) .* C_de); % % %���ɵ��óɱ� % % % %
C_buy = sum(sum(Pg) .* C_e);
C_loss = sum(sum(I .* (R * ones(1, T)))) * 400; % % %����ɱ� % % % %
C_ess = 250 * sum(sum(p_dch) + sum(p_ch)); % % %���ܵ��óɱ� % % % %
f2 = C_LA + C_ess - (B1 + B2) + C_buy + C_loss;
% +C_loss
%%%%����������ĳɱ�%%%%
n = sum(sum(u_dch + u_ch));
L = 0;
for D = 0:0.1:1
    L = L + n / (71470 * D ^ 4 - 170100 * D ^ 3 + 146400 * D ^ 2 - 56500 * D + 12230);
end
C_ESS1 = 819000 * 2 + 2953000 * 0.375;
C_ESS2 = 0;
for t = 1:1:10
    C_ESS2 = C_ESS2 + 69000 * 2 * ((1 + 0.015) / (1 + 0.09)) ^ t;
end
C_ESS = C_ESS1 + C_ESS2;
C_day = C_ESS * L; % % % % %����������ĳɱ� % % % % % %
%%%%%%%%%%%%%%%
%%%%%%������ۺ�����%%%%%%
peak_0 = PPload(9) + PPload(10) + PPload(11) + PPload(12) + PPload(13) + PPload(14) + PPload(19) + PPload(20);
peak_3 = Pload3(9) + Pload3(10) + Pload3(11) + Pload3(12) + Pload3(13) + Pload3(14) + Pload3(19) + Pload3(20);
lamda3 = (peak_0 - peak_3) / peak_0; % % % %��������������
C_gridup = 1;
delta_n3 = ((log10(1 + lamda3)) / log10(1 + 0.015)); % % % % % %�ӻ��������� % % % %
F1 = C_gridup * [1 - ((1 + 0.015) / (1 + 0.09)) ^ delta_n3]; % % % % % %�������������������� % % % %
T = 9; % % % %��ŵ����� % % % %
e_s = 1000; % % % %���糧�ȵ����� % % % % %
F2 = 0.5 * 0.375 * T * e_s; % % % % %���ٱ��óɱ� % % % % % %
toc %��ģʱ��
%% 5.�������
ops = sdpsettings('verbose', 1, 'solver', 'cplex');
sol = optimize(C, f2);
P = value(P);
p_pv = value(p_pv);
p_wt = value(p_wt);
p_ch = value(p_ch);
p_dch = value(p_dch);
Pg = value(Pg);
E_ess = value(E_ess);
R_ess = value(R_ess);
C_LA = value(C_LA);
C_loss = value(C_loss);
C_ess = value(C_ess);
C_day = value(C_day);
Lshift = value(Lshift);
Temp_shift = value(Temp_shift);
L = value(L);
S_IL1 = value(S_IL1);
S_IL2 = value(S_IL2);
f1_0 = value(f1_0);
f1_1 = value(f1_1);
f1_2 = value(f1_2);
f1_3 = value(f1_3);
F1 = value(F1);
delta_n3 = value(delta_n3);
V = value(V);
f1 = value(f1);
I = value(I);
%% 6.���������־
if sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end
%% ��ͼ
figure(1)
t = 1:1:24;
plot(t, 1.2 * PP_pv, '--k', 'linewidth', 2);
hold on
plot(t, 1.2 * P_PV, '-*r', 'linewidth', 2);
xlabel('ʱ��/h');
ylabel('�������/MW')
legend('Ԥ��ֵ', 'ʵ��ֵ');

figure(2)
t = 1:1:24;
plot(t, 1.2 * sum(PP_wt), '--', 'linewidth', 2);
hold on
plot(t, 1.2 * sum(P_WT), '-*', 'linewidth', 2);
xlabel('ʱ��/h');
ylabel('�������/MW')
legend('Ԥ��ֵ', 'ʵ��ֵ');

figure(3)
t = 1:1:24;
plot(t, sum(p_dch - p_ch), '-^k', 'linewidth', 1.8);
hold on
plot(t, Pg(33, :), '-dm', 'linewidth', 1.8);
xlabel('ʱ��/h');
ylabel('�й�����/MW');
legend('���ܳ�ŵ繦��', '��������');

figure(4)
t = 1:1:24;
plot(t, PPload, '-', 'linewidth', 1.5);
hold on
plot(t, PPload + sum(Lshift_old - Lshift) - S_IL1 - S_IL2, '-.', 'linewidth', 1.5);
hold on
plot(t, PPload + sum(p_dch - p_ch), '-d', 'linewidth', 1.5);
hold on
plot(t, PPload + sum(p_dch - p_ch) + sum(Lshift_old - Lshift) - S_IL1 - S_IL2, '-*', 'linewidth', 1.5);
xlabel('ʱ��/h');
ylabel('�����й�����/MW')
legend('ʵ��ֵ', '����һ', '������', '������');

figure(5)
t = 1:1:24;
plot(t, p_dch(1, :) - p_ch(1, :), '-.', 'linewidth', 1.8);
hold on
plot(t, p_dch(2, :) - p_ch(2, :), '-d', 'linewidth', 1.8);
xlabel('ʱ��/h');
ylabel('�й�����/MW');
legend('ESS1', 'ESS2');

figure(6)
t = 1:1:24;
plot(t, E_ess(1, :) / 1.03, '-.', 'linewidth', 1.8);
hold on
plot(t, E_ess(2, :) / 1.21, '-d', 'linewidth', 1.8);
xlabel('ʱ��/h');
ylabel('SOC');
legend('ESS1', 'ESS2');

figure(7)
t = 1:1:24;
PPPload = [1.77184023	1.760691425	1.985343246	2.164759614	2.525216573	2.836282757	2.717289933	2.912703959	2.938170633	3.200578131	3.93626314	3.680872774	3.53168293	3.215685929	2.903666471	2.626416093	2.33835944	2.582209352	3.063339069	3.465731913	3.01938202	2.591419695	2.254491338	1.784180138];
plot(t, PPPload, '-.d', 'linewidth', 1.5);
hold on
plot(t, PPload, '-*', 'linewidth', 1.5);
xlabel('ʱ��/h');
ylabel('�й�����/MW');
legend('����Ԥ��ֵ', '����ʵ��ֵ');

figure(8)
t = 1:1:24;
S_DR = zeros(5, 24);
S_DR(1:3, :) = Lshift_old - Lshift;
S_DR(4, :) = -S_IL1;
S_DR(5, :) = -S_IL2;
bar(t, (S_DR)', 'stacked');
xlabel('ʱ��/h');
ylabel('���ɹ���/MW');
legend('��ת�Ƹ���1', '��ת�Ƹ���2', '��ת�Ƹ���3', '����������1', '����������2');
toc %���ʱ��
