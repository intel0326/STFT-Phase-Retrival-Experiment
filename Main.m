
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 13 May. 2019.
%
%   �������e
%        ��������U���𒊏o�D
%        ���̐U������ʑ����Q�̎�@(GLA, GLA+ADMM, GLA+ADMM+prop, ��ʉ�ADMM)�Ő��肵�C��r����
%
%        ��������U���𒊏o�D
%        �U���Ƀm�C�Y�������Ă���ʑ����Q�̎�@(GLA, GLA+ADMM, GLA+ADMM+prop, ��ʉ�ADMM)�Ő��肵�C��r����
%
%   ���s���@
%       �R�}���h�E�B���h�E���ɁuMain�v�Ŏ��s
%       
%
%   �����l
%       Initialize.m�ɂď����l��ݒ肵�C./Variable/Initialize.mat�ɂĕۑ�
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% �p�X��ʂ�
addpath ./Tool
% �N���X�̌Ăяo��
ins_tool = tool();

% �����l����т����C���[�N�X�y�[�X�ɕۑ�
run('Initialize.m');
load('./Variable/Initialize.mat');

%��ڂ̎������J�n
fprintf('**********Experiment clear amp**********\n');

% admm�̃p�����[�^�� ( �� = 0.1, 0.2, 10, 100)
%for rho = 0.1:0.1:1.0
%for rho = [0.01, 0.1, 0.2, 10, 20, 100]
for rho = [0.1, 10]
    
    % 10��̃C�e���[�V�������ňʑ��ƃX�y�N�g����ۑ�
    [line, column]=size(spectrum);
    result_spe = zeros(line, column, 10);
    result_angle = zeros(line, column, 10);
    
    Mean_spe = zeros(line, column, 2);
    Var_spe = zeros(line, column, 2);
    
    Mean_angle = zeros(line, column, 2);
    Var_angle = zeros(line, column, 2);
    
    %�J�E���g
    count=1;
    
    % 10��̃C�e���[�V�����ŏ����ʑ��ɍ����ł邩����
    for q = 1:1:10
        fprintf('\n');
        fprintf('q = %d, rho = %d \n', q, rho);
        fprintf('\n');
        
        % �����_���Ȉʑ����擾
        phase_temp = rand(amp_FFTsize, frames);
        
        run('Experiment_clear_amp.m');
        save(sprintf('./Variable/result_rho_%.2f.mat', rho));
    end
    
    % ���ςƕ��U���o��
    Mean_spe(:, :, count) = mean(result_spe, 3);
    Var_spe(:, :, count) = var(result_spe, 0, 3);
    
    Mean_angle(:, :, count) = mean(result_angle, 3);
    Var_angle(:, :, count) = var(result_angle, 0, 3);
    
    count = count + 1;
    
    sell_angle = sell_angle + 1;
    sell_spe = sell_spe + 1;
    
end
    
%��ڂ̎������J�n
%fprintf('**********Experiment noise amp**********\n');
%run('Experiment_noise_amp.m');

% �p�X������
rmpath ./Tool


