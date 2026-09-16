function [NCI_source, source_name] = load_nci_source_boundary(data_file, w)
%LOAD_NCI_SOURCE_BOUNDARY Read inlet NCI boundary from the case workbook.
% If the workbook does not include NCI_source, use the historical default waveform.

time_steps = numel(w);
NCI_source = [];
source_name = 'default NCI waveform';

if exist(data_file,'file') == 2
    try
        boundary_table = readtable(data_file,'Sheet','Sheet1','VariableNamingRule','preserve');
    catch
        boundary_table = readtable(data_file,'Sheet','Sheet1');
    end

    var_names = boundary_table.Properties.VariableNames;
    var_names_lower = lower(string(var_names));
    nci_col = find(strcmpi(var_names,'NCI_source') | contains(var_names_lower,'nci'),1);
    if ~isempty(nci_col)
        if height(boundary_table) < time_steps
            error('NCI_source 数据行数不足：需要 %d 行，实际只有 %d 行。', ...
                time_steps, height(boundary_table));
        end

        values = table2array(boundary_table(1:time_steps,nci_col));
        if iscell(values) || isstring(values)
            values = str2double(values);
        end
        values = double(values(:));
        if any(~isfinite(values))
            error('NCI_source 列包含空值或非数值。');
        end

        NCI_source = values;
        [~,file_name,file_ext] = fileparts(data_file);
        source_name = [file_name,file_ext,':NCI_source'];
        return
    end
end

time_hour = w(:)/3600;
NCI_source_base = 0.30;
NCI_source_min = 0.2;
NCI_source_max = 0.3;
NCI_wave_amp1 = 0.12;
NCI_wave_amp2 = 0.07;
NCI_wave_amp3 = 0.04;
NCI_wave_period1 = 6;
NCI_wave_period2 = 2.4;
NCI_wave_period3 = 1.1;

NCI_source_raw = NCI_source_base ...
    + NCI_wave_amp1*sin(2*pi*time_hour/NCI_wave_period1) ...
    + NCI_wave_amp2*sin(2*pi*time_hour/NCI_wave_period2 + pi/5) ...
    + NCI_wave_amp3*sin(2*pi*time_hour/NCI_wave_period3 + pi/3);
NCI_source = NCI_source_min ...
    + (NCI_source_raw-min(NCI_source_raw))/(max(NCI_source_raw)-min(NCI_source_raw)) ...
    *(NCI_source_max-NCI_source_min);
end
