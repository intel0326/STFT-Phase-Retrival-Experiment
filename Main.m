
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Corded by R.Nakatsu (dragonstar30210326@gmail.com) on 13 May. 2019.
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
addpath ./OptimizeMethod
% �N���X�̌Ăяo��
ins_tool = tool();

% �����l����т����C���[�N�X�y�[�X�ɕۑ�
run('Initialize.m');
load('./Variable/Initialize.mat');

%��ڂ̎������J�n
fprintf('**********Experiment clear amp**********\n');

% admm�̃p�����[�^�� ( �� = 0.1, 0.2, 10, 100)
%for rho = 0.1:0.1:1.0
%for rho = [0.001, 0.01, 0.1, 10, 100]
for rho = [0.001, 0.01]
    
    fprintf('\n');
    fprintf('rho = %d \n', rho);
    fprintf('\n');
        
    run('Experiment_clear_amp.m');
    save(sprintf('./Variable/result_rho_%.2f.mat', rho));
    
end
    
%�m�C�Y����̐U���X�y�N�g���O�����ɑ΂���������J�n
%fprintf('**********Experiment noise amp**********\n');
%run('Experiment_noise_amp.m');

%���𕷂��Ƃ�
%sound(y, freq);
%sound(true_sound, freq);
%sound(GLA_sound, freq);
%sound(ADMM_sound, freq);
%sound(Prop_sound, freq);
%sound(gADMM_sound, freq);

%�O���t
%{
figure(200)
plot(gADMM_A_error(100:end))
legend('GLA','ADMM','Prop','PropGeneral');
xlim([0 900]);
%}

%�O���t�̃t�H���g�T�C�Y��13���炢


