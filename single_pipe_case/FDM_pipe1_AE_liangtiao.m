clear,clc
%% 0 程序说明
% 单管道算例：水力部分采用首末压力边界和AE压降方程反算质量流量，温度部分沿用FDM递推。
% 如需对比，先运行 FDM_pipe1_liangtiao.m 生成温度和动态水力保存文件。

%% 1 定义物理参数
data_c = 4200;         %工质的比热容
data_Ta = 263.15;      %环境温度
data_lamda_w = 0.67*5;   %管道的导热系数
data_lamda = 0.015;    %管道的粗糙系数
data_rho = 960;        %水的密度
data_L = 9250;         %管道长度(m)
data_D = 1.4;          %管道直径
data_A = 1.54;         %管道横截面积(m2)
k_loss = data_lamda_w/data_c/data_A/data_rho; %温度与NCI方程统一使用的热损系数
data_X = 1;            %热损碳责任由下游负荷承担的比例，取值范围0到1
data_X_compare = 0;    %对比工况：不把热损碳责任分摊到下游
NCI_source_base = 0.30; %入口NCI基准值
NCI_source_min = 0.2;   %入口NCI缩放后的最小值
NCI_source_max = 0.3;   %入口NCI缩放后的最大值
NCI_wave_amp1 = 0.12;   %入口NCI低频波动幅值
NCI_wave_amp2 = 0.07;   %入口NCI中频波动幅值
NCI_wave_amp3 = 0.04;   %入口NCI高频波动幅值
NCI_wave_period1 = 6;   %低频波动周期(h)
NCI_wave_period2 = 2.4; %中频波动周期(h)
NCI_wave_period3 = 1.1; %高频波动周期(h)
pressure_drop_scale = 0.0625; %AE低流量算例：首末压差在上一版1/4基础上继续缩小为1/4，即原始压差的1/16

%% 2 定义迭代信息
delta_t=180;
delta_x=300;
t_end=75600;
time_steps = t_end/delta_t + 1;
w=0:delta_t:t_end;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
time_col = w(:);
time_hour = time_col/3600;
NCI_source_raw = NCI_source_base ...
    + NCI_wave_amp1*sin(2*pi*time_hour/NCI_wave_period1) ...
    + NCI_wave_amp2*sin(2*pi*time_hour/NCI_wave_period2 + pi/5) ...
    + NCI_wave_amp3*sin(2*pi*time_hour/NCI_wave_period3 + pi/3);
NCI_source = NCI_source_min + (NCI_source_raw-min(NCI_source_raw))/(max(NCI_source_raw)-min(NCI_source_raw))*(NCI_source_max-NCI_source_min);

%% 3 设置分段数
length_yu = mod(data_L,delta_x);
if length_yu == 0
    deta_h_zheng=fix(data_L/delta_x);
else
    deta_h_zheng=fix(data_L/delta_x)+1;
end

%% 4 读取边界数据
data_file = fullfile(script_dir,'single pipe.xlsx');
if exist(data_file,'file') ~= 2
    error('未找到入口边界数据文件：%s', data_file);
end

try
    singlepipe = readtable(data_file,'Sheet','Sheet1','VariableNamingRule','preserve');
catch
    singlepipe = readtable(data_file,'Sheet','Sheet1');
end

if height(singlepipe) < time_steps
    error('single pipe.xlsx 数据行数不足：需要 %d 行，实际只有 %d 行。', time_steps, height(singlepipe));
end

p_out_file = fullfile(script_dir, sprintf('P_out_%d.mat', delta_x));
if exist(p_out_file,'file') ~= 2
    error('未找到压力边界数据文件：%s', p_out_file);
end
load(p_out_file,'P_out');
if size(P_out,1) < time_steps || size(P_out,2) < deta_h_zheng + 1
    error('P_out 维度不足：需要至少 %d x %d，实际为 %d x %d。', ...
        time_steps, deta_h_zheng + 1, size(P_out,1), size(P_out,2));
end
p_in_AE = P_out(1:time_steps,1);
p_out_terminal_raw_AE = P_out(1:time_steps,deta_h_zheng+1);
p_out_terminal_AE = p_in_AE - pressure_drop_scale.*(p_in_AE - p_out_terminal_raw_AE);
P_boundary_AE = [p_in_AE,p_out_terminal_AE];

%% 5 水力AE方程：压降与质量流量关系
% 对单管道，首末压力边界给出总压差；
% 管道压降方程 delta_p = K_f * m * abs(m) 反算整根管道的质量流量。
Kf_total = data_lamda*data_L/(2*data_rho*data_A^2*data_D);
delta_p_raw_AE = p_in_AE - p_out_terminal_raw_AE;
delta_p_AE = p_in_AE - p_out_terminal_AE;
mass_flow_scale_estimate_AE = sqrt(pressure_drop_scale);
M_boundary_AE = sign(delta_p_AE).*sqrt(abs(delta_p_AE)./Kf_total);
M_out_AE = repmat(M_boundary_AE,1,deta_h_zheng);
u_out_AE = M_out_AE/data_rho/data_A;
fprintf('AE首末压差缩放系数：%.4f\n', pressure_drop_scale);
fprintf('AE质量流量理论缩放系数：%.4f\n', mass_flow_scale_estimate_AE);

P_AE = zeros(time_steps,deta_h_zheng+1);
P_AE(:,1)=p_in_AE;
for i=1:deta_h_zheng
    if i==deta_h_zheng && length_yu~=0
        delta_x_i = length_yu;
    else
        delta_x_i = delta_x;
    end
    Kf_i = data_lamda*delta_x_i/(2*data_rho*data_A^2*data_D);
    P_AE(:,i+1)=P_AE(:,i)-Kf_i.*M_boundary_AE.*abs(M_boundary_AE);
end
P_terminal_residual_AE = P_AE(:,deta_h_zheng+1) - p_out_terminal_AE;

%% 6 温度初值和矩阵初始化
T_in0=90.1724638 + 273.15;
T_out0=90.10096985 + 273.15;

T_in_AE=zeros(time_steps,deta_h_zheng);
T_out_AE=zeros(time_steps,deta_h_zheng);
T_in_AE(:,1) = table2array(singlepipe(1:time_steps,3)) + 273.15;

T_10=(T_out0-T_in0)/data_L;
for i=1:deta_h_zheng
    if i==1
        T_out_AE(1,i)=T_in0+T_10*delta_x;
        T_in_AE(1,i+1)=T_out_AE(1,i);
    elseif i==deta_h_zheng
        if length_yu==0
            T_out_AE(1,i)=T_out_AE(1,i-1)+T_10*delta_x;
        else
            T_out_AE(1,i)=T_out_AE(1,i-1)+T_10*length_yu;
        end
    else
        T_out_AE(1,i)=T_out_AE(1,i-1)+T_10*delta_x;
        T_in_AE(1,i+1)=T_out_AE(1,i);
    end
end

ConH_T = zeros(deta_h_zheng,2);

%% 7 求解温度方程
for kk=2:time_steps
    T_state=zeros(2*deta_h_zheng,1);

    for i = 1:1:deta_h_zheng
        if i==deta_h_zheng && length_yu~=0
            delta_x_i = length_yu;
        else
            delta_x_i = delta_x;
        end

        ConH_T(i,1) =  k_loss/2 ...
            -  (M_out_AE(kk,i)+M_out_AE(kk-1,i)) / delta_x_i / data_A / data_rho/2;
        ConH_T(i,2) =  k_loss/2 ...
            +  (M_out_AE(kk,i)+M_out_AE(kk-1,i)) / delta_x_i / data_A / data_rho/2;
    end

    T_state(1)=T_in_AE(kk,1);
    for i = 1:1:deta_h_zheng
        T_state(i+deta_h_zheng) = (T_in_AE(kk-1,i) + T_out_AE(kk-1,i) ...
            - (ConH_T(i,1) * T_in_AE(kk-1,i) + ConH_T(i,2) * T_out_AE(kk-1,i) - 2*k_loss * data_Ta) ...
            * delta_t-(1 + ConH_T(i,1) * delta_t)*T_state(i))/(1 + ConH_T(i,2) * delta_t);
        T_state(i+1) = T_state(i+deta_h_zheng);
    end
    T_in_AE(kk,:)=T_state(1 : deta_h_zheng).';
    T_out_AE(kk,:)=T_state(deta_h_zheng+1 : 2*deta_h_zheng).';
end

%% 8 求解碳排放流NCI方程
pipe_node_x_NCI = zeros(1,deta_h_zheng+1);
for i=1:deta_h_zheng
    if i==deta_h_zheng && length_yu~=0
        delta_x_NCI_init = length_yu;
    else
        delta_x_NCI_init = delta_x;
    end
    pipe_node_x_NCI(i+1)=pipe_node_x_NCI(i)+delta_x_NCI_init;
end
u_init_NCI_AE = u_out_AE(1,1);
if u_init_NCI_AE <= 0
    error('NCI初始空间分布要求正向流速，当前AE初始流速为 %.6f m/s。', u_init_NCI_AE);
end

NCI_in_AE=zeros(time_steps,deta_h_zheng);
NCI_out_AE=zeros(time_steps,deta_h_zheng);
NCI_initial_profile_AE = NCI_source(1)*exp(data_X*k_loss*pipe_node_x_NCI/u_init_NCI_AE);
NCI_in_AE(:,1)=NCI_source;
NCI_in_AE(1,:)=NCI_initial_profile_AE(1:deta_h_zheng);
NCI_out_AE(1,:)=NCI_initial_profile_AE(2:deta_h_zheng+1);

for kk=2:time_steps
    for i=1:deta_h_zheng
        if i==deta_h_zheng && length_yu~=0
            delta_x_NCI = length_yu;
        else
            delta_x_NCI = delta_x;
        end
        u_bar_NCI = (u_out_AE(kk,i)+u_out_AE(kk-1,i))/2;
        con_NCI_in = -u_bar_NCI/delta_x_NCI - data_X*k_loss/2;
        con_NCI_out =  u_bar_NCI/delta_x_NCI - data_X*k_loss/2;
        NCI_out_AE(kk,i) = (NCI_in_AE(kk-1,i) + NCI_out_AE(kk-1,i) ...
                        - delta_t*(con_NCI_in*NCI_in_AE(kk-1,i) + con_NCI_out*NCI_out_AE(kk-1,i)) ...
                        - (1 + con_NCI_in*delta_t)*NCI_in_AE(kk,i))/(1 + con_NCI_out*delta_t);
        if i < deta_h_zheng
            NCI_in_AE(kk,i+1)=NCI_out_AE(kk,i);
        end
    end
end

NCI_in_AE_X0=zeros(time_steps,deta_h_zheng);
NCI_out_AE_X0=zeros(time_steps,deta_h_zheng);
NCI_initial_profile_AE_X0 = NCI_source(1)*exp(data_X_compare*k_loss*pipe_node_x_NCI/u_init_NCI_AE);
NCI_in_AE_X0(:,1)=NCI_source;
NCI_in_AE_X0(1,:)=NCI_initial_profile_AE_X0(1:deta_h_zheng);
NCI_out_AE_X0(1,:)=NCI_initial_profile_AE_X0(2:deta_h_zheng+1);

for kk=2:time_steps
    for i=1:deta_h_zheng
        if i==deta_h_zheng && length_yu~=0
            delta_x_NCI = length_yu;
        else
            delta_x_NCI = delta_x;
        end
        u_bar_NCI = (u_out_AE(kk,i)+u_out_AE(kk-1,i))/2;
        con_NCI_in = -u_bar_NCI/delta_x_NCI - data_X_compare*k_loss/2;
        con_NCI_out =  u_bar_NCI/delta_x_NCI - data_X_compare*k_loss/2;
        NCI_out_AE_X0(kk,i) = (NCI_in_AE_X0(kk-1,i) + NCI_out_AE_X0(kk-1,i) ...
                        - delta_t*(con_NCI_in*NCI_in_AE_X0(kk-1,i) + con_NCI_out*NCI_out_AE_X0(kk-1,i)) ...
                        - (1 + con_NCI_in*delta_t)*NCI_in_AE_X0(kk,i))/(1 + con_NCI_out*delta_t);
        if i < deta_h_zheng
            NCI_in_AE_X0(kk,i+1)=NCI_out_AE_X0(kk,i);
        end
    end
end

NCI_mean_in_AE = mean(NCI_source);
NCI_mean_out_AE = mean(NCI_out_AE(:,deta_h_zheng));
NCI_mean_out_AE_X0 = mean(NCI_out_AE_X0(:,deta_h_zheng));
NCI_terminal_diff_X1_X0 = NCI_out_AE(:,deta_h_zheng) - NCI_out_AE_X0(:,deta_h_zheng);
NCI_mean_increase_AE = NCI_mean_out_AE - NCI_mean_in_AE;
NCI_mean_increase_percent_AE = NCI_mean_increase_AE/NCI_mean_in_AE*100;
fprintf('AE入口NCI平均值：%.6f\n', NCI_mean_in_AE);
fprintf('AE管道末端出口NCI平均值：%.6f\n', NCI_mean_out_AE);
fprintf('AE管道末端出口NCI平均值(X=0)：%.6f\n', NCI_mean_out_AE_X0);
fprintf('AE出口相对入口NCI平均增量：%.6f (%.4f%%)\n', NCI_mean_increase_AE, NCI_mean_increase_percent_AE);
fprintf('AE末端NCI差值(X=1减X=0)平均值：%.8f\n', mean(NCI_terminal_diff_X1_X0));
fprintf('AE末端NCI差值(X=1减X=0)最大值：%.8f\n', max(NCI_terminal_diff_X1_X0));

figure(9)
plot(w,NCI_source,'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,NCI_out_AE(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
plot(w,NCI_out_AE_X0(:,deta_h_zheng),'-.','Color','#77AC30','LineWidth',1.5)
xlim([0 t_end])
legend('入口NCI','稳态AE末端NCI(X=1)','稳态AE末端NCI(X=0)')
xlabel('时间/s')
ylabel('NCI')

figure(12)
plot(w,NCI_terminal_diff_X1_X0,'-','Color','#7E2F8E','LineWidth',1.5)
xlim([0 t_end])
legend('末端NCI差值(X=1-X=0)')
xlabel('时间/s')
ylabel('NCI差值')

%% 9 计算热功率流和CEFR
Q_in_AE = data_c .* M_out_AE .* (T_in_AE - data_Ta);
Q_out_AE = data_c .* M_out_AE .* (T_out_AE - data_Ta);
Q_in_AE_kW = Q_in_AE/1000;
Q_out_AE_kW = Q_out_AE/1000;
CEFR_in_AE = NCI_in_AE .* Q_in_AE_kW;
CEFR_out_AE = NCI_out_AE .* Q_out_AE_kW;
CEFR_in_AE_X0 = NCI_in_AE_X0 .* Q_in_AE_kW;
CEFR_out_AE_X0 = NCI_out_AE_X0 .* Q_out_AE_kW;

Q_mean_in_AE = mean(Q_in_AE_kW(:,1));
Q_mean_out_AE = mean(Q_out_AE_kW(:,deta_h_zheng));
Q_mean_loss_AE = Q_mean_in_AE - Q_mean_out_AE;
CEFR_mean_in_AE = mean(CEFR_in_AE(:,1));
CEFR_mean_out_AE = mean(CEFR_out_AE(:,deta_h_zheng));
CEFR_mean_out_AE_X0 = mean(CEFR_out_AE_X0(:,deta_h_zheng));
CEFR_mean_loss_AE = CEFR_mean_in_AE - CEFR_mean_out_AE;
CEFR_integral_in_AE = trapz(w,CEFR_in_AE(:,1))/3600;
CEFR_integral_out_AE = trapz(w,CEFR_out_AE(:,deta_h_zheng))/3600;
CEFR_integral_diff_AE = CEFR_integral_in_AE - CEFR_integral_out_AE;
CEFR_integral_diff_percent_AE = CEFR_integral_diff_AE/CEFR_integral_in_AE*100;
fprintf('AE入口热功率流平均值(kW)：%.6f\n', Q_mean_in_AE);
fprintf('AE管道末端出口热功率流平均值(kW)：%.6f\n', Q_mean_out_AE);
fprintf('AE平均热损功率(kW)：%.6f\n', Q_mean_loss_AE);
fprintf('AE入口CEFR平均值(kgCO2/h)：%.6f\n', CEFR_mean_in_AE);
fprintf('AE管道末端出口CEFR平均值(kgCO2/h)：%.6f\n', CEFR_mean_out_AE);
fprintf('AE管道末端出口CEFR平均值(X=0, kgCO2/h)：%.6f\n', CEFR_mean_out_AE_X0);
fprintf('AE入口与出口CEFR平均差值(kgCO2/h)：%.6f\n', CEFR_mean_loss_AE);
fprintf('AE入口CEFR时间积分(kgCO2)：%.6f\n', CEFR_integral_in_AE);
fprintf('AE末端出口CEFR时间积分(kgCO2)：%.6f\n', CEFR_integral_out_AE);
fprintf('AE首末端CEFR积分差值(kgCO2)：%.6f (%.4f%%)\n', CEFR_integral_diff_AE, CEFR_integral_diff_percent_AE);

[~,middle_node_index_AE] = min(abs(pipe_node_x_NCI - data_L/2));
middle_node_x_AE = pipe_node_x_NCI(middle_node_index_AE);
if middle_node_index_AE <= deta_h_zheng
    NCI_middle_node_AE = NCI_in_AE(:,middle_node_index_AE);
    CEFR_middle_node_AE = CEFR_in_AE(:,middle_node_index_AE);
else
    NCI_middle_node_AE = NCI_out_AE(:,deta_h_zheng);
    CEFR_middle_node_AE = CEFR_out_AE(:,deta_h_zheng);
end
fprintf('AE三点NCI/CEFR对比使用X=%.0f，中间节点位置：%.2f m\n', data_X, middle_node_x_AE);

figure(10)
plot(w,CEFR_in_AE(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,CEFR_out_AE(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
plot(w,CEFR_out_AE_X0(:,deta_h_zheng),'-.','Color','#77AC30','LineWidth',1.5)
xlim([0 t_end])
legend('入口CEFR','稳态AE末端CEFR(X=1)','稳态AE末端CEFR(X=0)')
xlabel('时间/s')
ylabel('CEFR/(kgCO2/h)')

figure(11)
subplot(3,1,1)
plot(w,NCI_out_AE(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
hold on
plot(w,NCI_out_AE_X0(:,deta_h_zheng),'-.','Color','#77AC30','LineWidth',1.5)
xlim([0 t_end])
legend('末端NCI(X=1)','末端NCI(X=0)')
ylabel('NCI')
subplot(3,1,2)
plot(w,Q_out_AE_kW(:,deta_h_zheng),'-','Color','#0072BD','LineWidth',1.5)
xlim([0 t_end])
legend('末端热功率流')
ylabel('Q/kW')
subplot(3,1,3)
plot(w,CEFR_out_AE(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
hold on
plot(w,CEFR_out_AE_X0(:,deta_h_zheng),'-.','Color','#77AC30','LineWidth',1.5)
xlim([0 t_end])
legend('末端CEFR(X=1)','末端CEFR(X=0)')
xlabel('时间/s')
ylabel('CEFR/(kgCO2/h)')

figure(13)
subplot(2,1,1)
plot(w,NCI_in_AE(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,NCI_middle_node_AE,'-.','Color','#77AC30','LineWidth',1.5)
xlim([0 t_end])
legend('入口NCI','中间节点NCI')
ylabel('NCI')
subplot(2,1,2)
plot(w,NCI_in_AE(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,NCI_out_AE(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
xlim([0 t_end])
legend('入口NCI','末端出口NCI')
xlabel('时间/s')
ylabel('NCI')

figure(14)
subplot(2,1,1)
plot(w,CEFR_in_AE(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,CEFR_middle_node_AE,'-.','Color','#77AC30','LineWidth',1.5)
xlim([0 t_end])
legend('入口CEFR','中间节点CEFR')
ylabel('CEFR/(kgCO2/h)')
subplot(2,1,2)
plot(w,CEFR_in_AE(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,CEFR_out_AE(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
xlim([0 t_end])
legend('入口CEFR','末端出口CEFR')
xlabel('时间/s')
ylabel('CEFR/(kgCO2/h)')

save(fullfile(script_dir,'FDM_pipe1_AE_liangtiao_temperature.mat'), ...
    'T_in_AE','T_out_AE','M_out_AE','u_out_AE','P_AE','P_boundary_AE', ...
    'p_in_AE','p_out_terminal_raw_AE','p_out_terminal_AE','delta_p_raw_AE', ...
    'delta_p_AE','pressure_drop_scale','mass_flow_scale_estimate_AE','M_boundary_AE', ...
    'P_terminal_residual_AE','NCI_source','NCI_in_AE','NCI_out_AE', ...
    'NCI_in_AE_X0','NCI_out_AE_X0','Q_in_AE','Q_out_AE','Q_in_AE_kW','Q_out_AE_kW', ...
    'NCI_terminal_diff_X1_X0','CEFR_in_AE','CEFR_out_AE','CEFR_in_AE_X0','CEFR_out_AE_X0', ...
    'CEFR_integral_in_AE','CEFR_integral_out_AE','CEFR_integral_diff_AE','CEFR_integral_diff_percent_AE', ...
    'middle_node_index_AE','middle_node_x_AE','NCI_middle_node_AE','CEFR_middle_node_AE', ...
    'w','deta_h_zheng','Kf_total','k_loss','pipe_node_x_NCI');
save(fullfile(script_dir,'FDM_pipe1_AE_liangtiao_carbon.mat'), ...
    'NCI_source','NCI_in_AE','NCI_out_AE','NCI_in_AE_X0','NCI_out_AE_X0', ...
    'Q_in_AE','Q_out_AE','Q_in_AE_kW','Q_out_AE_kW','NCI_terminal_diff_X1_X0','CEFR_in_AE','CEFR_out_AE', ...
    'CEFR_integral_in_AE','CEFR_integral_out_AE','CEFR_integral_diff_AE','CEFR_integral_diff_percent_AE', ...
    'CEFR_in_AE_X0','CEFR_out_AE_X0','pressure_drop_scale', ...
    'data_X','data_X_compare','k_loss','pipe_node_x_NCI','middle_node_index_AE', ...
    'middle_node_x_AE','NCI_middle_node_AE','CEFR_middle_node_AE','w','deta_h_zheng');

%% 10 与 FDM_pipe1_liangtiao 温度对比
reference_file = fullfile(script_dir,'FDM_pipe1_liangtiao_temperature.mat');
if exist(reference_file,'file') ~= 2
    warning('未找到 %s。请先运行 FDM_pipe1_liangtiao.m 生成温度结果。', reference_file);
else
    ref = load(reference_file,'T_out','w','deta_h_zheng');
    T_out_ref = ref.T_out;
    w_ref = ref.w;
    deta_ref = ref.deta_h_zheng;

    T_terminal_diff = T_out_AE(:,deta_h_zheng) - T_out_ref(:,deta_ref);
    T_terminal_APE = abs(T_terminal_diff ./ T_out_ref(:,deta_ref))*100;
    T_terminal_MAPE = mean(T_terminal_APE);
    T_terminal_maxAPE = max(T_terminal_APE);
    fprintf('末端出口温度MAPE：%.6f%%\n', T_terminal_MAPE);
    fprintf('末端出口温度maxAPE：%.6f%%\n', T_terminal_maxAPE);

    figure(5)
    plot(w_ref,T_out_ref(:,deta_ref),'-','Color','#0072BD','LineWidth',1.5)
    hold on
    plot(w,T_out_AE(:,deta_h_zheng),'--','Color','#D95319','LineWidth',1.5)
    xlim([0 t_end])
    legend('FDM\_pipe1\_liangtiao末端出口温度','AE水力程序末端出口温度')
    xlabel('时间/s')
    ylabel('温度/K')

    figure(6)
    plot(w,T_terminal_diff,'-','Color','#77AC30','LineWidth',1.5)
    xlim([0 t_end])
    xlabel('时间/s')
    ylabel('温度差/K')
end

%% 11 与 FDM_pipe1_liangtiao 水力结果对比
hydraulic_file = fullfile(script_dir,'FDM_pipe1_liangtiao_dynamic_hydraulic.mat');
if exist(hydraulic_file,'file') ~= 2
    warning('未找到 %s。请先运行 FDM_pipe1_liangtiao.m 生成动态水力结果。', hydraulic_file);
else
    hyd = load(hydraulic_file,'P_out','M_out','w','deta_h_zheng');
    P_reference_database = hyd.P_out;
    M_dynamic = hyd.M_out;
    w_dynamic = hyd.w;
    deta_dynamic = hyd.deta_h_zheng;

    pressure_compare_node = ceil((min(deta_h_zheng,deta_dynamic)+2)/2);
    P_compare_reference = P_reference_database(:,pressure_compare_node);
    P_compare_AE = P_AE(:,pressure_compare_node);
    P_compare_diff = P_compare_AE - P_compare_reference;
    P_compare_APE = abs(P_compare_diff ./ P_compare_reference)*100;
    P_compare_MAPE = mean(P_compare_APE);
    P_compare_maxAPE = max(P_compare_APE);

    pipe_node_x = zeros(1,deta_h_zheng+1);
    for i=1:deta_h_zheng
        if i==deta_h_zheng && length_yu~=0
            delta_x_i = length_yu;
        else
            delta_x_i = delta_x;
        end
        pipe_node_x(i+1)=pipe_node_x(i)+delta_x_i;
    end
    P_endpoint_linear = repmat(P_reference_database(:,1),1,deta_h_zheng+1) ...
        - (P_reference_database(:,1)-P_reference_database(:,deta_h_zheng+1))*(pipe_node_x/data_L);
    P_database_linear_diff = P_reference_database(:,1:deta_h_zheng+1) - P_endpoint_linear;
    fprintf('压力数据库由首末压力线性重构的最大误差：%.6e Pa\n', max(abs(P_database_linear_diff(:))));

    M_terminal_dynamic = M_dynamic(:,deta_dynamic);
    M_terminal_AE = M_out_AE(:,deta_h_zheng);
    M_terminal_diff = M_terminal_AE - M_terminal_dynamic;
    M_terminal_APE = abs(M_terminal_diff ./ M_terminal_dynamic)*100;
    M_terminal_MAPE = mean(M_terminal_APE);
    M_terminal_maxAPE = max(M_terminal_APE);

    fprintf('AE末端压力边界残差最大绝对值：%.6f Pa\n', max(abs(P_terminal_residual_AE)));
    fprintf('中间节点压力MAPE：%.6f%%\n', P_compare_MAPE);
    fprintf('中间节点压力maxAPE：%.6f%%\n', P_compare_maxAPE);
    fprintf('末端质量流量MAPE：%.6f%%\n', M_terminal_MAPE);
    fprintf('末端质量流量maxAPE：%.6f%%\n', M_terminal_maxAPE);

    figure(7)
    subplot(2,1,1)
    plot(w_dynamic,P_compare_reference,'-','Color','#0072BD','LineWidth',1.5)
    hold on
    plot(w,P_compare_AE,'--','Color','#D95319','LineWidth',1.5)
    xlim([0 t_end])
    legend(sprintf('原压力数据库节点%d',pressure_compare_node),sprintf('压差1/16稳态AE节点%d',pressure_compare_node))
    ylabel('压力/Pa')
    subplot(2,1,2)
    plot(w,P_compare_diff,'-','Color','#77AC30','LineWidth',1.5)
    xlim([0 t_end])
    xlabel('时间/s')
    ylabel('压力差/Pa')

    figure(8)
    subplot(2,1,1)
    plot(w_dynamic,M_terminal_dynamic,'-','Color','#0072BD','LineWidth',1.5)
    hold on
    plot(w,M_terminal_AE,'--','Color','#D95319','LineWidth',1.5)
    xlim([0 t_end])
    legend('动态FDM末端质量流量','稳态AE末端质量流量')
    ylabel('质量流量/(kg/s)')
    subplot(2,1,2)
    plot(w,M_terminal_diff,'-','Color','#77AC30','LineWidth',1.5)
    xlim([0 t_end])
    xlabel('时间/s')
    ylabel('质量流量差/(kg/s)')
end
