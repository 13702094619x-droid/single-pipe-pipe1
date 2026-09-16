clear,clc
script_dir = fileparts(mfilename('fullpath'));
data_file = fullfile(script_dir,'single pipe.xlsx');
%% 程序说明
%热网
%两个节点、一根管道的算例
%SAS二维嵌入程序,算例1模拟
%% 1.1 定义热网参数
data_c = 4200;         %工质的比热容
data_Ta = 263.15;       %环境温度
data_lamda_w = 0.67*5;   %管道的导热系数
data_lamda = 0.015;    %管道的粗糙系数
data_R = 0.35;         %管段的热阻
data_rho = 960;        %水的密度
data_L = 9250;           %管道b的长度(m)
data_D = 1.4;          %管道b的直径
data_A = 1.54;         %管道b的横截面积(m2)
pressure_drop_scale = 0.0625; %与FDM-AE低流量算例一致：首末压差在上一版1/4基础上继续缩小为1/4，即原始压差的1/16
%% 1.2 定义迭代信息：时间与空间步长
delta_t=180;  %时间步长为180s
delta_x=300;  %空间步长为300m
t_end=75600; %仿真结束时间
%% 1.3 t=0时迭代初值
T_in0=90.1724638 + 273.15;
T_out0=90.10096985 + 273.15;
w=0:delta_t:t_end;
ww=0:90:t_end;
[NCI_source_boundary,NCI_source_boundary_name] = load_nci_source_boundary(data_file,w(:));
fprintf('SAS-NCI入口边界来源：%s\n', NCI_source_boundary_name);
mass_reference_file = fullfile(script_dir,'FDM_pipe1_AE_liangtiao_temperature.mat');
M_in_source_time = w(:);
M_in_source_name = 'single pipe.xlsx';
if exist(mass_reference_file,'file') == 2
    mass_reference = load(mass_reference_file,'M_boundary_AE','M_out_AE','w');
    if isfield(mass_reference,'M_boundary_AE')
        M_in1 = mass_reference.M_boundary_AE(:);
    elseif isfield(mass_reference,'M_out_AE')
        M_in1 = mass_reference.M_out_AE(:,1);
    else
        error('历史结果文件中未找到 M_boundary_AE 或 M_out_AE：%s', mass_reference_file);
    end
    if isfield(mass_reference,'w')
        M_in_source_time = mass_reference.w(:);
    else
        M_in_source_time = (0:length(M_in1)-1)'*delta_t;
    end
    M_in_source_name = 'FDM_pipe1_AE_liangtiao_temperature.mat';
else
    M_in1=xlsread(data_file,1,'B2:B422')*960;
end
M_in = interp1(M_in_source_time,M_in1(:),w(:),'linear','extrap');
M_in=M_in(:); 
fprintf('SAS质量流量入口数据来源：%s\n', M_in_source_name);

%%判断是否能分段整除
length_yu = mod(data_L,delta_x); 
if length_yu == 0
    deta_h_zheng=fix(data_L/delta_x);%% 能分段整除，段数为长度除空间步长
else
    deta_h_zheng=fix(data_L/delta_x)+1;%% 不能分段整除，段数为长度除空间步长后＋1   
end
M_in=M_in.*ones(1,deta_h_zheng);
T_10=(T_out0-T_in0)/data_L;
T_in=zeros(t_end/delta_t+1,deta_h_zheng); %每段的末端温度储存矩
T_out=zeros(t_end/delta_t+1,deta_h_zheng); %每段的末端温度储存矩阵
M_out=zeros(t_end/delta_t+1,deta_h_zheng);
M_out(1,:)=M_in(1,:);

for i=1:deta_h_zheng %% 稳态温度分布
    if i==1
        T_out(1,i)=T_in0+T_10*delta_x;
        T_in(1,i+1)=T_out(1,i);
    elseif i==deta_h_zheng
        if length_yu==0
            T_out(1,i)=T_out(1,i-1)+T_10*delta_x;
        else
            T_out(1,i)=T_out(1,i-1)+T_10*length_yu;
        end
    else
        T_out(1,i)=T_out(1,i-1)+T_10*delta_x;
        T_in(1,i+1)=T_out(1,i);
    end
end
T_in1=xlsread(data_file,1,'C2:C422')+273.15;
T_in(:,1) = interp1(w,T_in1,w,'linear')';
con1=data_c*data_rho*data_A;
ConH_T = zeros(deta_h_zheng,2);
con3 =  2*data_lamda_w/data_c/data_A/data_rho;
con4=data_lamda/2/data_rho/data_A/data_D;%动量守恒方程常数
if delta_x==300
    load('P_out_300.mat')
elseif delta_x==500
    load('P_out_500.mat')
elseif delta_x==800
    load('P_out_800.mat')
elseif delta_x==50
    load('P_out_50.mat')    
elseif delta_x==200
    load('P_out_200.mat')    
elseif delta_x==100  
    load('P_out_100.mat')
end
P_out_raw = P_out;
P_in_boundary_sas = P_out_raw(:,1);
P_out = P_in_boundary_sas - pressure_drop_scale.*(P_in_boundary_sas - P_out_raw);
P_terminal_raw_sas = P_out_raw(:,deta_h_zheng+1);
P_terminal_scaled_sas = P_out(:,deta_h_zheng+1);
fprintf('SAS首末压差缩放系数：%.4f\n', pressure_drop_scale);
kk=1;
error_max_M=zeros(t_end/delta_t,1);
P_x=zeros(t_end/delta_t,deta_h_zheng);
P_x1=zeros(t_end/delta_t,deta_h_zheng);
T1=zeros(t_end/delta_t,deta_h_zheng,4,4);
tic
for t=delta_t:delta_t:t_end
    T=zeros(2*deta_h_zheng,1);
    kk=kk+1;
    for i=1:deta_h_zheng-1
        P_x(kk-1,i)=(P_out(kk-1,i+1)-P_out(kk-1,i))/delta_x;
    end
    P_x(kk-1,deta_h_zheng)=(P_out(kk-1,deta_h_zheng+1)-P_out(kk-1,deta_h_zheng))/length_yu;
    for i=1:deta_h_zheng-1
        P_x1(kk-1,i)=(P_out(kk,i+1)-P_out(kk,i)-P_x(kk-1,i)*delta_x)/delta_x/delta_t;
    end
    P_x1(kk-1,deta_h_zheng)=(P_out(kk,deta_h_zheng+1)-P_out(kk,deta_h_zheng)-P_x(kk-1,deta_h_zheng)*length_yu)/length_yu/delta_t;
    H_M=zeros(31,deta_h_zheng);
    for i=1:deta_h_zheng
        H_M(1,i)=M_out(kk-1,i);
        H_M(2,i)=-data_A*P_x(kk-1,i)-con4*H_M(1,i)^2;
    end
    error_M=zeros(length(H_M(2,:)),1);
    for i=1:length(H_M(2,:))
        error_M(i)=con4*(2*H_M(2,i)*H_M(1,i)*delta_t+(H_M(2,i)^2)*delta_t^2);
    end
    error_max_M(kk-1)=max(abs(error_M));
    if error_max_M(kk-1) < 1e-6
        M_out(kk,:)  = M_out(kk-1,:)+H_M(2,1:deta_h_zheng).*delta_t;
    else
        for i=1:deta_h_zheng
            H_M(3,i)=-con4*H_M(2,i)*H_M(1,i)-data_A*P_x1(kk-1,i)/2;
        end
        error_M=zeros(length(H_M(2,:)),1);
        for i=1:length(H_M(1,:))
            error_M(i)=con4*(2*H_M(3,i)*H_M(1,i)*delta_t^2+2*H_M(2,i)*H_M(3,i)*delta_t^3+H_M(2,i)^2*delta_t^2+H_M(3,i)^2*delta_t^4);
        end
        error_max_M(kk-1)=max(abs(error_M));
        if error_max_M(kk-1) < 1e-6
            M_out(kk,:)  = M_out(kk-1,:)+H_M(2,1 : deta_h_zheng).*delta_t+H_M(3,1 : deta_h_zheng).*delta_t^2;
        else
            for i=1:deta_h_zheng
                H_M(4,i)=-con4*(H_M(2,i)^2+2*H_M(3,i)*H_M(1,i))/3;
            end
        end
        error_M=zeros(length(H_M(2,:)),1);
        for i=1:length(H_M(1,:))
            error_M(i)=con4*(H_M(3,i)^2*delta_t^4+2*H_M(4,i)*H_M(1,i)*delta_t^3+2*H_M(2,i)*H_M(3,i)*delta_t^3+2*H_M(2,i)*H_M(4,i)*delta_t^4+2*H_M(3,i)*H_M(4,i)*delta_t^5+H_M(4,i)^2*delta_t^6);
        end
        error_max_M(kk-1)=max(abs(error_M));
        if error_max_M(kk-1) < 1e-6
            M_out(kk,:)  = M_out(kk-1,:)+H_M(2,1 : deta_h_zheng).*delta_t+H_M(3,1 : deta_h_zheng).*delta_t^2+H_M(4,1 : deta_h_zheng).*delta_t^3;
        else
            for i=1:deta_h_zheng
                H_M(5,i)=-con4*(2*H_M(2,i)*H_M(3,i)+2*H_M(4,i)*H_M(1,i))/4;
            end
            error_M=zeros(length(H_M(2,:)),1);
            for i=1:length(H_M(1,:))
                error_M(i)=con4*(H_M(3,i)^2*delta_t^4+2*H_M(5,i)*H_M(1,i)*delta_t^4+2*H_M(2,i)*H_M(4,i)*delta_t^4+2*H_M(3,i)*H_M(4,i)*delta_t^5+2*H_M(4,i)*H_M(5,i)*delta_t^7+H_M(5,i)^2*delta_t^8+H_M(4,i)^2*delta_t^6);
            end
            error_max_M(kk-1)=max(abs(error_M));
        end
        if error_max_M(kk-1) < 1e-6
            M_out(kk,:)  = M_out(kk-1,:)+H_M(2,1 : deta_h_zheng).*delta_t+H_M(3,1 : deta_h_zheng).*delta_t^2+H_M(4,1 : deta_h_zheng).*delta_t^3+H_M(5,1 : deta_h_zheng).*delta_t^4;
        else
            j=6;
            while 1
                for i=1:deta_h_zheng
                    for a=1:j-1
                        H_M(j,i)=H_M(j,i)-con4*(H_M(a,i)*H_M(j-a,i))/(j-1);
                    end
                end
                error_M=zeros(length(H_M(2,:)),1);
                for i=1:length(H_M(1,:))
                    for m=0:j-1
                        for n=0:m
                            error_M(i)=error_M(i)+con4*(H_M(m+1,i)*H_M(j-n,i)*delta_t^(j-n+m-1));
                        end
                    end
                end
                error_max_M(kk-1)=max(abs(error_M));
                if j >= 8
                    break
                else
                    j=j+1;
                end
            end
            for i=1:length(H_M(:,1))
                M_out(kk,:) =M_out(kk,:)+H_M(i,:).*delta_t^(i-1);
            end
%             j=6;
%             for i=1:deta_h_zheng
%                 for a=1:j-1
%                     H_M(j,i)=H_M(j,i)-con4*(H_M(a,i)*H_M(j-a,i))/(j-1);
%                 end
%             end
%             error_M=zeros(length(H_M(2,:)),1);
%             for i=1:length(H_M(1,:))
%                 for m=0:j-1
%                     for n=0:m
%                         error_M(i)=error_M(i)+H_M(m+1,i)*H_M(j-n,i)*delta_t^(j-n+m-1);
%                     end
%                 end
%             end
%             error_max_M(kk-1)=max(abs(error_M));
%             for i=1:length(H_M(:,1))
%                 M_out(kk,:) =M_out(kk,:)+H_M(i,:).*delta_t^(i-1);
%             end
        end
    end
    for i = 1:1:deta_h_zheng-1
        ConH_T(i,1) =  data_lamda_w / data_c / data_A / data_rho / 2 -  (M_out(kk,i)+M_out(kk-1,i)) / delta_x / data_A / data_rho/2;  
        ConH_T(i,2) =  data_lamda_w / data_c / data_A / data_rho / 2 +  (M_out(kk,i)+M_out(kk-1,i)) / delta_x / data_A / data_rho/2;   
    end
    ConH_T(deta_h_zheng,1) =  data_lamda_w / data_c / data_A / data_rho / 2 -  (M_out(kk,deta_h_zheng)+M_out(kk-1,deta_h_zheng)) / length_yu / data_A / data_rho/2;
    ConH_T(deta_h_zheng,2) =  data_lamda_w / data_c / data_A / data_rho / 2 +  (M_out(kk,deta_h_zheng)+M_out(kk-1,deta_h_zheng)) / length_yu / data_A / data_rho/2;   
    T(1)=T_in(kk,1);
    for i = 1:1:deta_h_zheng
        T(i+deta_h_zheng) = (T_in(kk-1,i) + T_out(kk-1,i) - (ConH_T(i,1) * T_in(kk-1,i) + ConH_T(i,2) * T_out(kk-1,i) - con3 * data_Ta) * delta_t-(1 + ConH_T(i,1) * delta_t)*T(i))/(1 + ConH_T(i,2) * delta_t);
        T(i+1) = T(i+deta_h_zheng);
    end
    T_in(kk,:)=T(1 : deta_h_zheng).';
    T_out(kk,:)=T(deta_h_zheng+1 : 2*deta_h_zheng).';
    M1=zeros(deta_h_zheng,1);
    M_ed=zeros(deta_h_zheng,1);
    %% 温度计算
    Mat_T_left=zeros(9,9);
    Mat_T_right=zeros(9,1);
    %% 能量守恒方程
    for i=deta_h_zheng
        T1(kk-1,i,1,1)=T_in(kk-1,i);
        M1(i)=(M_out(kk,i)-M_out(kk-1,i))/delta_t;
        M_ed(i)=M_out(kk-1,i);
    end
    for i=deta_h_zheng
        Mat_T_left(1,1)=con1;
        Mat_T_left(1,4)=data_c*M_ed(i);
        Mat_T_right(1)=-data_lamda_w*(T1(kk-1,i,1,1)-data_Ta);
        Mat_T_left(2,5)=con1;
        Mat_T_left(2,7)=2*data_c*M_ed(i);
        Mat_T_left(2,4)=data_lamda_w;
        Mat_T_left(3,2)=2*con1;
        Mat_T_left(3,4)=data_c*M1(i);    
        Mat_T_left(3,5)=data_c*M_ed(i); 
        Mat_T_left(3,1)=data_lamda_w;
        Mat_T_left(4,8)=con1;
        Mat_T_left(4,9)=3*data_c*M_ed(i); 
        Mat_T_left(4,7)=data_lamda_w;
        Mat_T_left(5,3)=3*con1;
        Mat_T_left(5,6)=data_c*M_ed(i); 
        Mat_T_left(5,5)=data_c*M1(i);  
        Mat_T_left(5,2)=data_lamda_w;  
        Mat_T_left(6,6)=2*con1;
        Mat_T_left(6,8)=2*data_c*M_ed(i); 
        Mat_T_left(6,7)=2*data_c*M1(i); 
        Mat_T_left(6,5)=data_lamda_w; 
        Mat_T_left(7,1)=delta_t; %% 时间边界条件
        Mat_T_left(7,2)=delta_t^2;
        Mat_T_left(7,3)=delta_t^3;
        Mat_T_right(7)=T_in(kk,i)-T_in(kk-1,i);
        if i==deta_h_zheng
            Mat_T_left(8,4)=length_yu; %% 空间边界条件
            Mat_T_left(8,7)=length_yu^2;
            Mat_T_left(8,9)=length_yu^3;
            Mat_T_left(9,1)=delta_t; % 
            Mat_T_left(9,2)=delta_t^2;
            Mat_T_left(9,3)=delta_t^3;
            Mat_T_left(9,4)=length_yu;
            Mat_T_left(9,5)=delta_t*length_yu;
            Mat_T_left(9,6)=length_yu*delta_t^2;
            Mat_T_left(9,7)=length_yu^2;
            Mat_T_left(9,8)=delta_t*length_yu^2;
            Mat_T_left(9,9)=length_yu^3;
        else
            Mat_T_left(8,4)=delta_x; %% 空间边界条件
            Mat_T_left(8,7)=delta_x^2;
            Mat_T_left(8,9)=delta_x^3;
            Mat_T_left(9,1)=delta_t; % 
            Mat_T_left(9,2)=delta_t^2;
            Mat_T_left(9,3)=delta_t^3;
            Mat_T_left(9,4)=delta_x;
            Mat_T_left(9,5)=delta_t*delta_x;
            Mat_T_left(9,6)=delta_x*delta_t^2;
            Mat_T_left(9,7)=delta_x^2;
            Mat_T_left(9,8)=delta_t*delta_x^2;
            Mat_T_left(9,9)=delta_x^3;
        end
        Mat_T_right(8)=T_out(kk-1,i)-T_in(kk-1,i);
        Mat_T_right(9)=T_out(kk,i)-T_in(kk-1,i);
        T_solve=Mat_T_left\Mat_T_right;
        T1(kk-1,i,1,2)=T_solve(1);
        T1(kk-1,i,1,3)=T_solve(2);
        T1(kk-1,i,1,4)=T_solve(3);
        T1(kk-1,i,2,1)=T_solve(4);
        T1(kk-1,i,2,2)=T_solve(5);
        T1(kk-1,i,2,3)=T_solve(6);
        T1(kk-1,i,3,1)=T_solve(7);
        T1(kk-1,i,3,2)=T_solve(8);
        T1(kk-1,i,4,1)=T_solve(9);
    end   
end
toc

M_reference_terminal = [];
M_sas_reference_diff = [];
M_sas_reference_MAE = [];
M_sas_reference_max_abs_diff = [];
if exist('mass_reference','var')
    if isfield(mass_reference,'M_out_AE')
        M_reference_terminal_raw = mass_reference.M_out_AE(:,min(size(mass_reference.M_out_AE,2),deta_h_zheng));
    elseif isfield(mass_reference,'M_boundary_AE')
        M_reference_terminal_raw = mass_reference.M_boundary_AE(:);
    else
        M_reference_terminal_raw = [];
    end
    if ~isempty(M_reference_terminal_raw)
        M_reference_terminal = interp1(M_in_source_time,M_reference_terminal_raw(:),w(:),'linear','extrap');
        M_sas_reference_diff = M_out(:,deta_h_zheng)-M_reference_terminal;
        M_sas_reference_MAE = mean(abs(M_sas_reference_diff));
        M_sas_reference_max_abs_diff = max(abs(M_sas_reference_diff));
    end
end
fprintf('SAS入口质量流量平均值(kg/s)：%.6f\n', mean(M_in(:,1)));
fprintf('SAS末端质量流量平均值(kg/s)：%.6f\n', mean(M_out(:,deta_h_zheng)));
if ~isempty(M_reference_terminal)
    fprintf('AE参考末端质量流量平均值(kg/s)：%.6f\n', mean(M_reference_terminal));
    fprintf('SAS-AE末端质量流量平均绝对差(kg/s)：%.6f\n', M_sas_reference_MAE);
    fprintf('SAS-AE末端质量流量最大绝对差(kg/s)：%.6f\n', M_sas_reference_max_abs_diff);
end

%% 2 SAS-NCI求解
% 方程：dNCI/dt + u(t)*dNCI/dx = X*k_loss*NCI
% 其中 k_loss 与温度方程中的散热系数保持一致。
k_loss = data_lamda_w/data_c/data_A/data_rho;
data_X_NCI = 1;                     % 热损碳责任系数，X=1表示热损责任计入下游
fdm_nci_file = fullfile(script_dir,'FDM_pipe1_AE_liangtiao_carbon.mat');
fdm_nci_reference_loaded = false;
fdm_nci_data = struct();

if exist(fdm_nci_file,'file') == 2
    fdm_nci_data = load(fdm_nci_file, ...
        'NCI_source','NCI_in_AE','NCI_out_AE','k_loss','data_X','w','deta_h_zheng');
    if isfield(fdm_nci_data,'NCI_source')
        fdm_nci_reference_loaded = true;
        if isfield(fdm_nci_data,'k_loss')
            k_loss = fdm_nci_data.k_loss;
        end
        if isfield(fdm_nci_data,'data_X')
            data_X_NCI = fdm_nci_data.data_X;
        end
    end
end

beta_NCI = data_X_NCI*k_loss;

time_steps = length(w);
NCI_source_sas = NCI_source_boundary(:);
pipe_segment_length = delta_x*ones(1,deta_h_zheng);
if length_yu ~= 0
    pipe_segment_length(end) = length_yu;
end
pipe_node_x_NCI = [0,cumsum(pipe_segment_length)];

NCI_in_sas = zeros(time_steps,deta_h_zheng);
NCI_out_sas = zeros(time_steps,deta_h_zheng);
NCI_factor_sas = zeros(time_steps-1,deta_h_zheng,3,3);

u_init_NCI = max(M_out(1,1)/data_rho/data_A,1e-9);
NCI_initial_profile_sas = NCI_source_sas(1) ...
    * exp(data_X_NCI*k_loss*pipe_node_x_NCI/u_init_NCI);
NCI_in_sas(:,1) = NCI_source_sas;
if fdm_nci_reference_loaded ...
        && isfield(fdm_nci_data,'NCI_in_AE') ...
        && isfield(fdm_nci_data,'NCI_out_AE') ...
        && size(fdm_nci_data.NCI_in_AE,2) >= deta_h_zheng ...
        && size(fdm_nci_data.NCI_out_AE,2) >= deta_h_zheng
    NCI_in_sas(1,:) = fdm_nci_data.NCI_in_AE(1,1:deta_h_zheng);
    NCI_out_sas(1,:) = fdm_nci_data.NCI_out_AE(1,1:deta_h_zheng);
    NCI_initial_profile_sas = [NCI_in_sas(1,1),NCI_out_sas(1,:)];
else
    NCI_in_sas(1,:) = NCI_initial_profile_sas(1:deta_h_zheng);
    NCI_out_sas(1,:) = NCI_initial_profile_sas(2:deta_h_zheng+1);
end

for kk = 2:time_steps
    NCI_in_sas(kk,1) = NCI_source_sas(kk);
    for i = 1:deta_h_zheng
        dx_NCI = pipe_segment_length(i);
        rho00 = NCI_in_sas(kk-1,i);
        rhoL_prev = NCI_out_sas(kk-1,i);
        rho_in_now = NCI_in_sas(kk,i);

        u_prev_NCI = max(M_out(kk-1,i)/data_rho/data_A,1e-9);
        u_now_NCI = max(M_out(kk,i)/data_rho/data_A,1e-9);
        u00_NCI = u_prev_NCI;
        u01_NCI = (u_now_NCI-u_prev_NCI)/delta_t;
        u10_NCI = 0;                % 单管道量调节假设同一时刻流速无空间差异

        Mat_NCI_left = zeros(5,5);
        Mat_NCI_right = zeros(5,1);

        % 未知量顺序：[rho01; rho02; rho10; rho11; rho20]
        % 时间边界：rho(dt,0)=当前时刻入口NCI
        Mat_NCI_left(1,1) = delta_t;
        Mat_NCI_left(1,2) = delta_t^2;
        Mat_NCI_right(1) = rho_in_now-rho00;

        % 空间边界：rho(0,dx)=上一时刻管段出口NCI
        Mat_NCI_left(2,3) = dx_NCI;
        Mat_NCI_left(2,5) = dx_NCI^2;
        Mat_NCI_right(2) = rhoL_prev-rho00;

        % PDE常数项
        Mat_NCI_left(3,1) = 1;
        Mat_NCI_left(3,3) = u00_NCI;
        Mat_NCI_right(3) = beta_NCI*rho00;

        % PDE时间一阶项
        Mat_NCI_left(4,1) = -beta_NCI;
        Mat_NCI_left(4,2) = 2;
        Mat_NCI_left(4,3) = u01_NCI;
        Mat_NCI_left(4,4) = u00_NCI;

        % PDE空间一阶项
        Mat_NCI_left(5,3) = u10_NCI-beta_NCI;
        Mat_NCI_left(5,4) = 1;
        Mat_NCI_left(5,5) = 2*u00_NCI;

        if rcond(Mat_NCI_left) < 1e-12
            NCI_solve = pinv(Mat_NCI_left)*Mat_NCI_right;
        else
            NCI_solve = Mat_NCI_left\Mat_NCI_right;
        end

        rho01 = NCI_solve(1);
        rho02 = NCI_solve(2);
        rho10 = NCI_solve(3);
        rho11 = NCI_solve(4);
        rho20 = NCI_solve(5);

        rho_out_now = rho00 ...
            + rho01*delta_t ...
            + rho02*delta_t^2 ...
            + rho10*dx_NCI ...
            + rho11*delta_t*dx_NCI ...
            + rho20*dx_NCI^2;
        rho_out_now = max(rho_out_now,0);

        NCI_factor_sas(kk-1,i,1,1) = rho00;
        NCI_factor_sas(kk-1,i,1,2) = rho01;
        NCI_factor_sas(kk-1,i,1,3) = rho02;
        NCI_factor_sas(kk-1,i,2,1) = rho10;
        NCI_factor_sas(kk-1,i,2,2) = rho11;
        NCI_factor_sas(kk-1,i,3,1) = rho20;

        NCI_out_sas(kk,i) = rho_out_now;
        if i < deta_h_zheng
            NCI_in_sas(kk,i+1) = rho_out_now;
        end
    end
end

Q_in_sas = data_c.*M_out.*(T_in-data_Ta);
Q_out_sas = data_c.*M_out.*(T_out-data_Ta);
Q_in_sas_kW = Q_in_sas/1000;
Q_out_sas_kW = Q_out_sas/1000;
CEFR_in_sas = NCI_in_sas.*Q_in_sas_kW;
CEFR_out_sas = NCI_out_sas.*Q_out_sas_kW;

% 准动态碳流：动态能流 + 稳态碳熵分配模型。
% 单管道无中途取热时，源侧输入碳熵流沿管道分配给实际有损热功率。
Q_node_sas_kW = zeros(time_steps,deta_h_zheng+1);
Q_node_sas_kW(:,1) = Q_in_sas_kW(:,1);
Q_node_sas_kW(:,2:end) = Q_out_sas_kW;
CEFR_quasi_source_sas = NCI_source_sas.*Q_in_sas_kW(:,1);
CEFR_quasi_node_sas = repmat(CEFR_quasi_source_sas,1,deta_h_zheng+1);
NCI_quasi_node_sas = CEFR_quasi_node_sas./max(abs(Q_node_sas_kW),eps);

NCI_mean_in_sas = mean(NCI_in_sas(:,1));
NCI_mean_out_sas = mean(NCI_out_sas(:,deta_h_zheng));
NCI_integral_in_sas = trapz(w,NCI_in_sas(:,1))/3600;
NCI_integral_out_sas = trapz(w,NCI_out_sas(:,deta_h_zheng))/3600;
CEFR_mean_in_sas = mean(CEFR_in_sas(:,1));
CEFR_mean_out_sas = mean(CEFR_out_sas(:,deta_h_zheng));
CEFR_integral_in_sas = trapz(w,CEFR_in_sas(:,1))/3600;
CEFR_integral_out_sas = trapz(w,CEFR_out_sas(:,deta_h_zheng))/3600;
CEFR_integral_diff_sas = CEFR_integral_in_sas-CEFR_integral_out_sas;
CEFR_integral_diff_percent_sas = CEFR_integral_diff_sas/CEFR_integral_in_sas*100;

[~,middle_node_index_sas] = min(abs(pipe_node_x_NCI-data_L/2));
middle_node_x_sas = pipe_node_x_NCI(middle_node_index_sas);
if middle_node_index_sas <= deta_h_zheng
    NCI_middle_node_sas = NCI_in_sas(:,middle_node_index_sas);
    CEFR_middle_node_sas = CEFR_in_sas(:,middle_node_index_sas);
else
    NCI_middle_node_sas = NCI_out_sas(:,deta_h_zheng);
    CEFR_middle_node_sas = CEFR_out_sas(:,deta_h_zheng);
end
NCI_middle_node_quasi_sas = NCI_quasi_node_sas(:,middle_node_index_sas);
NCI_terminal_quasi_sas = NCI_quasi_node_sas(:,end);
CEFR_middle_node_quasi_sas = CEFR_quasi_node_sas(:,middle_node_index_sas);
CEFR_terminal_quasi_sas = CEFR_quasi_node_sas(:,end);
NCI_integral_middle_quasi_sas = trapz(w,NCI_middle_node_quasi_sas)/3600;
NCI_integral_terminal_quasi_sas = trapz(w,NCI_terminal_quasi_sas)/3600;
CEFR_integral_middle_quasi_sas = trapz(w,CEFR_middle_node_quasi_sas)/3600;
CEFR_integral_terminal_quasi_sas = trapz(w,CEFR_terminal_quasi_sas)/3600;

NCI_terminal_fdm_compare = [];
NCI_terminal_sas_fdm_diff = [];
NCI_terminal_sas_fdm_MAE = [];
NCI_terminal_sas_fdm_max_abs_diff = [];
NCI_terminal_sas_fdm_MAPE = [];
NCI_terminal_sas_fdm_maxAPE = [];
NCI_integral_terminal_fdm_compare = [];
NCI_integral_terminal_sas_fdm_diff = [];

if fdm_nci_reference_loaded && isfield(fdm_nci_data,'NCI_out_AE')
    fdm_terminal_col = size(fdm_nci_data.NCI_out_AE,2);
    fdm_compare_time = (0:size(fdm_nci_data.NCI_out_AE,1)-1)'*delta_t;
    if isfield(fdm_nci_data,'w')
        fdm_compare_time = fdm_nci_data.w(:);
    end
    NCI_terminal_fdm_compare = interp1(fdm_compare_time, ...
        fdm_nci_data.NCI_out_AE(:,fdm_terminal_col),w(:),'linear','extrap');
    NCI_terminal_sas_fdm_diff = NCI_out_sas(:,deta_h_zheng) - NCI_terminal_fdm_compare;
    NCI_terminal_sas_fdm_MAE = mean(abs(NCI_terminal_sas_fdm_diff));
    NCI_terminal_sas_fdm_max_abs_diff = max(abs(NCI_terminal_sas_fdm_diff));
    NCI_terminal_sas_fdm_APE = abs(NCI_terminal_sas_fdm_diff./max(abs(NCI_terminal_fdm_compare),eps))*100;
    NCI_terminal_sas_fdm_MAPE = mean(NCI_terminal_sas_fdm_APE);
    NCI_terminal_sas_fdm_maxAPE = max(NCI_terminal_sas_fdm_APE);
    NCI_integral_terminal_fdm_compare = trapz(w,NCI_terminal_fdm_compare)/3600;
    NCI_integral_terminal_sas_fdm_diff = NCI_integral_out_sas - NCI_integral_terminal_fdm_compare;
end

fprintf('SAS-NCI热损系数k_loss：%.10e 1/s\n', k_loss);
if fdm_nci_reference_loaded
    fprintf('SAS-NCI入口与初始条件来源：FDM_pipe1_AE_liangtiao_carbon.mat\n');
else
    fprintf('SAS-NCI入口与初始条件来源：备用正弦边界\n');
end
fprintf('SAS-NCI入口NCI平均值：%.8f\n', NCI_mean_in_sas);
fprintf('SAS-NCI末端出口NCI平均值：%.8f\n', NCI_mean_out_sas);
fprintf('SAS-NCI入口NCI时间积分：%.8f h\n', NCI_integral_in_sas);
fprintf('SAS-NCI末端出口NCI时间积分：%.8f h\n', NCI_integral_out_sas);
fprintf('SAS入口CEFR平均值(kgCO2/h)：%.6f\n', CEFR_mean_in_sas);
fprintf('SAS管道末端出口CEFR平均值(kgCO2/h)：%.6f\n', CEFR_mean_out_sas);
fprintf('SAS入口CEFR时间积分(kgCO2)：%.6f\n', CEFR_integral_in_sas);
fprintf('SAS末端出口CEFR时间积分(kgCO2)：%.6f\n', CEFR_integral_out_sas);
fprintf('SAS首末端CEFR积分差值(kgCO2)：%.6f (%.4f%%)\n', CEFR_integral_diff_sas, CEFR_integral_diff_percent_sas);
fprintf('SAS-NCI中间节点位置：%.2f m\n', middle_node_x_sas);
fprintf('准动态中间节点NCI时间积分：%.8f h\n', NCI_integral_middle_quasi_sas);
fprintf('准动态末端NCI时间积分：%.8f h\n', NCI_integral_terminal_quasi_sas);
fprintf('准动态中间节点CEFR时间积分(kgCO2)：%.6f\n', CEFR_integral_middle_quasi_sas);
fprintf('准动态末端CEFR时间积分(kgCO2)：%.6f\n', CEFR_integral_terminal_quasi_sas);
if fdm_nci_reference_loaded && ~isempty(NCI_terminal_fdm_compare)
    fprintf('SAS-FDM末端NCI平均绝对差：%.8f\n', NCI_terminal_sas_fdm_MAE);
    fprintf('SAS-FDM末端NCI最大绝对差：%.8f\n', NCI_terminal_sas_fdm_max_abs_diff);
    fprintf('SAS-FDM末端NCI MAPE：%.6f%%\n', NCI_terminal_sas_fdm_MAPE);
    fprintf('SAS-FDM末端NCI maxAPE：%.6f%%\n', NCI_terminal_sas_fdm_maxAPE);
    fprintf('FDM末端NCI时间积分：%.8f h\n', NCI_integral_terminal_fdm_compare);
    fprintf('SAS-FDM末端NCI积分差：%.8f h\n', NCI_integral_terminal_sas_fdm_diff);
end

quasi_curve_color_sas = '#7E2F8E';

figure(20)
subplot(2,1,1)
plot(w,NCI_in_sas(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,NCI_middle_node_sas,'-.','Color','#77AC30','LineWidth',1.5)
plot(w,NCI_middle_node_quasi_sas,'--','Color',quasi_curve_color_sas,'LineWidth',1.7)
xlim([0 t_end])
legend('入口NCI','动态中间节点NCI','准动态中间节点NCI')
ylabel('NCI')
subplot(2,1,2)
plot(w,NCI_in_sas(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,NCI_out_sas(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
plot(w,NCI_terminal_quasi_sas,'--','Color',quasi_curve_color_sas,'LineWidth',1.7)
xlim([0 t_end])
legend('入口NCI','动态末端出口NCI','准动态末端出口NCI')
xlabel('时间/s')
ylabel('NCI')

figure(21)
subplot(2,1,1)
plot(w,CEFR_in_sas(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,CEFR_middle_node_sas,'-.','Color','#77AC30','LineWidth',1.5)
plot(w,CEFR_middle_node_quasi_sas,'--','Color',quasi_curve_color_sas,'LineWidth',1.7)
xlim([0 t_end])
legend('入口CEFR','动态中间节点CEFR','准动态中间节点CEFR')
ylabel('CEFR/(kgCO2/h)')
subplot(2,1,2)
plot(w,CEFR_in_sas(:,1),'--','Color','#0072BD','LineWidth',1.5)
hold on
plot(w,CEFR_out_sas(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
plot(w,CEFR_terminal_quasi_sas,'--','Color',quasi_curve_color_sas,'LineWidth',1.7)
xlim([0 t_end])
legend('入口CEFR','动态末端出口CEFR','准动态末端出口CEFR')
xlabel('时间/s')
ylabel('CEFR/(kgCO2/h)')

if fdm_nci_reference_loaded && ~isempty(NCI_terminal_fdm_compare)
    figure(22)
    subplot(2,1,1)
    plot(w,NCI_out_sas(:,deta_h_zheng),'-','Color','#D95319','LineWidth',1.5)
    hold on
    plot(w,NCI_terminal_fdm_compare,'--','Color','#0072BD','LineWidth',1.5)
    xlim([0 t_end])
    legend('SAS末端NCI','FDM-AE末端NCI')
    ylabel('NCI')
    subplot(2,1,2)
    plot(w,NCI_terminal_sas_fdm_diff,'-','Color','#77AC30','LineWidth',1.5)
    xlim([0 t_end])
    xlabel('时间/s')
    ylabel('SAS-FDM')
end

save(fullfile(script_dir,'SAS_Pipe1_NCI_results.mat'), ...
    'NCI_source_sas','NCI_in_sas','NCI_out_sas','NCI_factor_sas', ...
    'NCI_mean_in_sas','NCI_mean_out_sas','NCI_integral_in_sas','NCI_integral_out_sas', ...
    'CEFR_in_sas','CEFR_out_sas','CEFR_mean_in_sas','CEFR_mean_out_sas', ...
    'CEFR_integral_in_sas','CEFR_integral_out_sas','CEFR_integral_diff_sas', ...
    'CEFR_integral_diff_percent_sas','Q_in_sas','Q_out_sas','Q_in_sas_kW','Q_out_sas_kW', ...
    'NCI_middle_node_sas','CEFR_middle_node_sas','middle_node_index_sas','middle_node_x_sas', ...
    'NCI_quasi_node_sas','CEFR_quasi_node_sas','CEFR_quasi_source_sas','Q_node_sas_kW', ...
    'NCI_middle_node_quasi_sas','NCI_terminal_quasi_sas', ...
    'CEFR_middle_node_quasi_sas','CEFR_terminal_quasi_sas', ...
    'NCI_integral_middle_quasi_sas','NCI_integral_terminal_quasi_sas', ...
    'CEFR_integral_middle_quasi_sas','CEFR_integral_terminal_quasi_sas', ...
    'NCI_terminal_fdm_compare','NCI_terminal_sas_fdm_diff','NCI_terminal_sas_fdm_MAE', ...
    'NCI_terminal_sas_fdm_max_abs_diff','NCI_terminal_sas_fdm_MAPE','NCI_terminal_sas_fdm_maxAPE', ...
    'NCI_integral_terminal_fdm_compare','NCI_integral_terminal_sas_fdm_diff', ...
    'fdm_nci_reference_loaded','fdm_nci_file','pipe_node_x_NCI','pipe_segment_length', ...
    'M_in1','M_in','M_out','mass_reference_file','M_in_source_name', ...
    'M_reference_terminal','M_sas_reference_diff','M_sas_reference_MAE','M_sas_reference_max_abs_diff', ...
    'P_out','P_out_raw','P_in_boundary_sas','P_terminal_raw_sas','P_terminal_scaled_sas','pressure_drop_scale', ...
    'k_loss','data_X_NCI','w','deta_h_zheng','data_L');
