clear,clc
%% 0 程序说明
%热网
%两个节点、一根管道的算例
%FDM计算程序,算例1模
%质量调节和质调节的核心区别，质调节里质量流量假设是一个恒定值，但是在质量调节中，先要用压力算质量流量，然后用质量流量计算温度的迭代


%核心算法流程图

%初始化阶段

%├── 1 定义物理参数（c, ρ, λ, D, L等）
%├── 2 定义迭代信息：时间与空间步长（Δt=180s, Δx=300m）
%├── 3 设置分段数
%├── 4 t=0时刻的温度设置
%├── 5 读取边界数据
%│   ├── 入口流量 M_in(t) [Excel]
%│   ├── 入口温度 T_in(t) [Excel]
%│   └── 压力分布 P_out [.mat文件]
%└── 6 数据结构初始化和常数预计算

%计算阶段

%时间循环 (t = Δt → t_end)
%├── 7 求解动量方程
%│   ├── 7.1 构建非线性方程组 F(M) = 0
%│   ├── 7.2 使用 fsolve 求解
%│
%├── 8 计算热传导系数
%│   ├── 基于新的流量 M
%│   └── 更新 ConH_T 矩阵
%│
%└── 9 求解能量方程
%    ├── 设置入口边界 T(1)
%    ├── 逐段计算温度
%    └── 更新 T_in, T_out 矩阵
%后处理
%└── 绘制出口温度时程曲线
%% 1 定义物理参数

data_c = 4200;         %工质的比热容
data_Ta = 263.15;       %环境温度
data_lamda_w = 0.67*5;   %管道的导热系数，与SAS_Pipe1保持一致
data_lamda = 0.015;    %管道的粗糙系数
data_R = 0.35;         %管段的热阻
data_rho = 960;        %水的密度
data_L = 9250;           %管道b的长度(m)
data_D = 1.4;          %管道b的直径
data_A = 1.54;         %管道b的横截面积(m2)
  
k_loss = data_lamda_w/data_c/data_A/data_rho; %温度与NCI方程统一使用的热损系数
data_X = 1;            %热损碳责任由下游负荷承担的比例，取值范围0到1
NCI_source_base = 0.30; %入口NCI基准值
NCI_source_min = 0.2;   %入口NCI缩放后的最小值
NCI_source_max = 0.3;   %入口NCI缩放后的最大值
NCI_wave_amp1 = 0.12;   %入口NCI低频波动幅值
NCI_wave_amp2 = 0.07;   %入口NCI中频波动幅值
NCI_wave_amp3 = 0.04;   %入口NCI高频波动幅值
NCI_wave_period1 = 6;   %低频波动周期(h)
NCI_wave_period2 = 2.4; %中频波动周期(h)
NCI_wave_period3 = 1.1; %高频波动周期(h)
flow_scale = 1;         %基准工况流量系数

%% 2 定义迭代信息：时间与空间步长

delta_t=180;  %时间步长为180s
delta_x=300;  %空间步长为300m
t_end=75600; %仿真结束时间

%% 3 设置分段数  % 尽可能使用统一的标准段长（300m），只在最后处理余数

length_yu = mod(data_L,delta_x); %  计算余数 mod是计算余数的函数
if length_yu == 0
    deta_h_zheng=fix(data_L/delta_x);%% 能分段整除，段数为长度除空间步长
else
    deta_h_zheng=fix(data_L/delta_x)+1;%% 不能分段整除，段数为长度除空间步长后＋1   
end
% 具体计算示例

%情况1：能够整除
%假设管道长度 = 9000m，空间步长 = 300m
%length_yu = mod(9000, 300) = 0
%deta_h_zheng = fix(9000/300) = 30

%分段示意：
%|--300m--|--300m--|--300m--| ... |--300m--|
%  段1      段2      段3     ...    段30

%情况2：不能整除（实际情况）
%管道长度 = 9250m，空间步长 = 300m
%length_yu = mod(9250, 300) = 250m
%deta_h_zheng = fix(9250/300) + 1 = 30 + 1 = 31

%分段示意：
%|--300m--|--300m--|--300m--| ... |--300m--|--250m--|
%   段1      段2      段3     ...    段30     段31

%% 4 t=0时 温度设置

T_in0=90.1724638 + 273.15;%初始温度设置
T_out0=90.10096985 + 273.15;%
%% 5 读取边界数据
time_steps = t_end/delta_t + 1;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

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

M_in = table2array(singlepipe(1:time_steps,2))*data_rho*flow_scale;
T_in=zeros(time_steps,deta_h_zheng); 
T_in(:,1) = table2array(singlepipe(1:time_steps,3)) + 273.15;

p_out_file = fullfile(script_dir, sprintf('P_out_%d.mat', delta_x));
if exist(p_out_file,'file') ~= 2
    error('未找到压力边界数据文件：%s', p_out_file);
end
load(p_out_file,'P_out');
if size(P_out,1) < time_steps || size(P_out,2) < deta_h_zheng + 1
    error('P_out 维度不足：需要至少 %d x %d，实际为 %d x %d。', ...
        time_steps, deta_h_zheng + 1, size(P_out,1), size(P_out,2));
end


%这是读取数据的第二种方法建立了matlab数据库 当我们需要读取数据时直接调用先前建立好的matlab数据集，这个方法程序运行起来比较快

%这段是压力的数据库，我们要用压力的值去计算质量流量，这里也是质量调节和质调节的核心区别，质调节里质量流量假设是一个恒定值，但是在质量调节中，先要用压力算质量流量

%% 6 数据结构初始化和常数预计算

w=0:delta_t:t_end;%创建时间点数组：[0, 180, 360, ..., 75600]
              %共421个时间点（0到75600秒，步长180秒）
              %用于后续绘图的x轴坐标
time_col = w(:);
time_hour = time_col/3600;
NCI_source_raw = NCI_source_base ...
    + NCI_wave_amp1*sin(2*pi*time_hour/NCI_wave_period1) ...
    + NCI_wave_amp2*sin(2*pi*time_hour/NCI_wave_period2 + pi/5) ...
    + NCI_wave_amp3*sin(2*pi*time_hour/NCI_wave_period3 + pi/3);
NCI_source = NCI_source_min + (NCI_source_raw-min(NCI_source_raw))/(max(NCI_source_raw)-min(NCI_source_raw))*(NCI_source_max-NCI_source_min);

    %每段的末端温度储存矩
T_out=zeros(time_steps,deta_h_zheng); 
%维度：421 × 31（时间步数 × 空间段数）
%T_in(i,j)：第i个时刻、第j段的入口温度
%T_out(i,j)：第i个时刻、第j段的出口温度

M_out=zeros(time_steps,deta_h_zheng);
%维度：421 × 31
%M_out(i,j)：第i个时刻、第j段的出口质量流率
u_out=zeros(time_steps,deta_h_zheng);
%u_out(i,j)：第i个时刻、第j段的出口流速
NCI_in=zeros(time_steps,deta_h_zheng);
NCI_out=zeros(time_steps,deta_h_zheng);
%NCI_in(i,j)、NCI_out(i,j)：第i个时刻、第j段入口和出口NCI

M_out(1,:)=M_in(1,:);
u_out(1,:)=M_out(1,:)/data_rho/data_A;

pipe_node_x_NCI = zeros(1,deta_h_zheng+1);
for i=1:deta_h_zheng
    if i==deta_h_zheng && length_yu~=0
        delta_x_NCI_init = length_yu;
    else
        delta_x_NCI_init = delta_x;
    end
    pipe_node_x_NCI(i+1)=pipe_node_x_NCI(i)+delta_x_NCI_init;
end
u_init_NCI = u_out(1,1);
if u_init_NCI <= 0
    error('NCI初始空间分布要求正向流速，当前初始流速为 %.6f m/s。', u_init_NCI);
end
NCI_initial_profile = NCI_source(1)*exp(data_X*k_loss*pipe_node_x_NCI/u_init_NCI);
NCI_in(:,1)=NCI_source;
NCI_in(1,:)=NCI_initial_profile(1:deta_h_zheng);
NCI_out(1,:)=NCI_initial_profile(2:deta_h_zheng+1);
T_10=(T_out0-T_in0)/data_L;

ConH_T = zeros(deta_h_zheng,2);
con3 =  2*k_loss;%公式30中的环境温度源项系数
con4=data_lamda/2/data_rho/data_A/data_D;     %λ/(2*ρ*A*D)    公式（25）中摩擦阻力的系数

for i=1:deta_h_zheng
    if i==1
        T_out(1,i)=T_in0+T_10*delta_x;
        T_in(1,i+1)=T_out(1,i);
    elseif i==deta_h_zheng
        if length_yu==0%如果余数等于0
            T_out(1,i)=T_out(1,i-1)+T_10*delta_x;%正常分段
        else
            T_out(1,i)=T_out(1,i-1)+T_10*length_yu;%否则余数加一段
        end
    else
        T_out(1,i)=T_out(1,i-1)+T_10*delta_x;
        T_in(1,i+1)=T_out(1,i);%上一阶段的输出等于下一阶段的输入
    end
end
%%  7 求解动量方程（公式24）
    %%  7.1 构建非线性方程组 F(M) = 0
kk=1;
tic%开始计时  
%%计算过程这段代码实现了一个​​动态热力管网（DHN）的数值求解过程​​，主要用于计算管道中的质量流率（M_out）和温度分布（T_in, T_out）随时间的变化。

for t=delta_t:delta_t:t_end%开始求解质量流率从delta_t到t_end，按时间步长逐步求解。
    T=zeros(2*deta_h_zheng,1);
    kk=kk+1;
    x0_H = M_out(kk-1,:); % 初始猜测：上一时间步的质量流率

  
 eqnH = @(x)[% 公式（24）
     (x(1:deta_h_zheng-1)-M_out(kk-1,1:deta_h_zheng-1)-data_A*delta_t.*((P_out(kk-1,1:deta_h_zheng-1)+P_out(kk,1:deta_h_zheng-1))/2-(P_out(kk-1,2:deta_h_zheng)+P_out(kk,2:deta_h_zheng))/2)./(delta_x).'+0.25*con4.*delta_t.*(x(1:deta_h_zheng-1).^2+M_out(kk-1,1:deta_h_zheng-1).^2+2.*M_out(kk-1,1:deta_h_zheng-1).*x(1:deta_h_zheng-1))).';
        (x(deta_h_zheng)-M_out(kk-1,deta_h_zheng)-data_A*delta_t*((P_out(kk-1,deta_h_zheng)+P_out(kk,deta_h_zheng))/2-(P_out(kk-1,deta_h_zheng+1)+P_out(kk,deta_h_zheng+1))/2)/length_yu+0.25*con4*delta_t*(x(deta_h_zheng)^2+M_out(kk-1,deta_h_zheng)^2+2*M_out(kk-1,deta_h_zheng)*x(deta_h_zheng)));
        ]; % 定义非线性方程组（质量守恒方程）
 %fsolve 是求解函数中有平方 立方的一种方法
 % 例子 1
 %求解单个方程：
 %求解方程：x^2 - 4 = 0

 %fun = @(x) x^2 - 4;
 %x0 = 1;  % 初始猜测
 %x = fsolve(fun, x0)
 
 % 结果：x = 2

 %例子2
 %求解方程组：
 
 % x1^2 + x2^2 = 4
 % x1 - x2 = 1
 %fun = @(x) [x(1)^2 + x(2)^2 - 4; 
 %             x(1) - x(2) - 1];
 %x0 = [1; 1];  % 初始猜测
 %x = fsolve(fun, x0)
 % 结果：x = [1.5811; 0.5811]
 
 %% 7.2 使用 fsolve 求解
    x = fsolve(eqnH, x0_H);% % 求解非线性方程组
    M_out(kk,:)=x;            %质量守恒方程离散化为隐式格式，包含：压力梯度项（(P_out(kk-1, :) - P_out(kk, :)) / delta_x）摩擦阻力项（con4 * M_out^2）
    u_out(kk,:)=M_out(kk,:)/data_rho/data_A; %由质量流率计算各管段流速

%%  8   计算热传导系数 （公式31 32）
    for i = 1:1:deta_h_zheng-1%对温度进行求解   ​​目标​​：计算管道在某一时刻（kk）的温度分布，包括入口温度（T_in）和出口温度（T_out）。
        ConH_T(i,1) =  k_loss/2 -  (M_out(kk,i)+M_out(kk-1,i)) / delta_x / data_A / data_rho/2;%热传导系数矩阵，分别表示第i个分段入口和出口的热传导系数  
        ConH_T(i,2) =  k_loss/2 +  (M_out(kk,i)+M_out(kk-1,i)) / delta_x / data_A / data_rho/2;   
    end
    ConH_T(deta_h_zheng,1) =  k_loss/2 -  (M_out(kk,deta_h_zheng)+M_out(kk-1,deta_h_zheng)) / length_yu / data_A / data_rho/2;
    ConH_T(deta_h_zheng,2) =  k_loss/2 +  (M_out(kk,deta_h_zheng)+M_out(kk-1,deta_h_zheng)) / length_yu / data_A / data_rho/2;   
    T(1)=T_in(kk,1);
%% 9  求解能量方程（公式34）
    for i = 1:1:deta_h_zheng%DHN中热水温度沿管道的动态传播由​​一阶对流-热损失方程​​描述，下面这个方程就是对这个方程的分析
        T(i+deta_h_zheng) = (T_in(kk-1,i) + T_out(kk-1,i) - (ConH_T(i,1) * T_in(kk-1,i) + ConH_T(i,2) * T_out(kk-1,i) - con3 * data_Ta) ...
                            * delta_t-(1 + ConH_T(i,1) * delta_t)*T(i))/(1 + ConH_T(i,2) * delta_t);%公式（26）到
        T(i+1) = T(i+deta_h_zheng);%上一阶段的输出等于下一阶段的输入
    end
    T_in(kk,:)=T(1 : deta_h_zheng).';%T  注意T是一个时刻的输入输出矩阵 上面是输入 下面是输出   把T上面矩阵的值 赋值到Tin
    T_out(kk,:)=T(deta_h_zheng+1 : 2*deta_h_zheng).';%把T下面矩阵的值 赋值到Tout

    %% 9.1 求解NCI方程
    for i = 1:1:deta_h_zheng
        if i==deta_h_zheng && length_yu~=0
            delta_x_NCI = length_yu;
        else
            delta_x_NCI = delta_x;
        end
        u_bar_NCI = (u_out(kk,i)+u_out(kk-1,i))/2;
        con_NCI_in = -u_bar_NCI/delta_x_NCI - data_X*k_loss/2;
        con_NCI_out =  u_bar_NCI/delta_x_NCI - data_X*k_loss/2;
        NCI_out(kk,i) = (NCI_in(kk-1,i) + NCI_out(kk-1,i) ...
                        - delta_t*(con_NCI_in*NCI_in(kk-1,i) + con_NCI_out*NCI_out(kk-1,i)) ...
                        - (1 + con_NCI_in*delta_t)*NCI_in(kk,i))/(1 + con_NCI_out*delta_t);
        if i < deta_h_zheng
            NCI_in(kk,i+1)=NCI_out(kk,i);
        end
    end
end
%% 10 画图  
toc

load(fullfile(script_dir,'T_out_300.mat'),'T_out_300')

figure(2)
plot(w,T_out(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
%参数解释：

%w：X轴数据（时间向量）

%w = 0:180:75600
%表示 [0, 180, 360, ..., 75600] 秒


%T_out(:,deta_h_zheng)：Y轴数据（温度）

%T_out(:,31)：所有时刻的第31段出口温度
%这是管道最末端的温度
%维度：421×1的向量


%'-'：线型

%实线（solid line）


%'Color', '#D95319'：颜色设置

%十六进制颜色代码
%#D95319 = 橙色/橘红色
%相当于 RGB(217, 83, 25)


%'LineWidth', 1.5：线宽

%1.5倍的默认线宽
%使曲线更清晰可见
hold on
plot(w,T_out_300(:,31),'--','Color','#0072BD','LineWidth',1.5)
xlim([0 t_end])

figure(3)
plot(w,NCI_source,'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,NCI_out(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
xlim([0 t_end])
legend('入口NCI','管道末端出口NCI')
xlabel('时间/s')
ylabel('NCI')

NCI_mean_in = mean(NCI_source);
NCI_mean_out = mean(NCI_out(:,deta_h_zheng));
NCI_mean_increase = NCI_mean_out - NCI_mean_in;
NCI_mean_increase_percent = NCI_mean_increase/NCI_mean_in*100;
fprintf('入口NCI平均值：%.6f\n', NCI_mean_in);
fprintf('管道末端出口NCI平均值：%.6f\n', NCI_mean_out);
fprintf('出口相对入口平均增量：%.6f (%.4f%%)\n', NCI_mean_increase, NCI_mean_increase_percent);

%% 11 计算热功率流和CEFR
Q_in = data_c .* M_out .* (T_in - data_Ta);
Q_out = data_c .* M_out .* (T_out - data_Ta);
Q_in_kW = Q_in/1000;
Q_out_kW = Q_out/1000;
CEFR_in = NCI_in .* Q_in_kW;
CEFR_out = NCI_out .* Q_out_kW;

Q_mean_in = mean(Q_in_kW(:,1));
Q_mean_out = mean(Q_out_kW(:,deta_h_zheng));
Q_mean_loss = Q_mean_in - Q_mean_out;
CEFR_mean_in = mean(CEFR_in(:,1));
CEFR_mean_out = mean(CEFR_out(:,deta_h_zheng));
CEFR_mean_loss = CEFR_mean_in - CEFR_mean_out;
fprintf('入口热功率流平均值(kW)：%.6f\n', Q_mean_in);
fprintf('管道末端出口热功率流平均值(kW)：%.6f\n', Q_mean_out);
fprintf('平均热损功率(kW)：%.6f\n', Q_mean_loss);
fprintf('入口CEFR平均值(kgCO2/h)：%.6f\n', CEFR_mean_in);
fprintf('管道末端出口CEFR平均值(kgCO2/h)：%.6f\n', CEFR_mean_out);
fprintf('入口与出口CEFR平均差值(kgCO2/h)：%.6f\n', CEFR_mean_loss);

figure(4)
plot(w,CEFR_in(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,CEFR_out(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
xlim([0 t_end])
legend('入口CEFR','管道末端出口CEFR')
xlabel('时间/s')
ylabel('CEFR/(kgCO2/h)')

save(fullfile(script_dir,'FDM_pipe1_liangtiao_temperature.mat'),'T_in','T_out','w','deta_h_zheng','flow_scale','k_loss','pipe_node_x_NCI');
save(fullfile(script_dir,'FDM_pipe1_liangtiao_dynamic_hydraulic.mat'),'P_out','M_out','u_out','w','deta_h_zheng','flow_scale','k_loss','pipe_node_x_NCI');
