clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 22 Apr. 2019.
%
%   ��������U���𒊏o�D
%   �U���Ƀm�C�Y�������Ă���ʑ����Q�̎�@(GLA, GLA+ADMM)�Ő��肵�C��r����
%
%   ���s���@
%       �R�}���h�E�B���h�E���ɁuMain�v�Ŏ��s
%
%   �����l
%       Initialize.m�ɂď����l��ݒ肵�C./Variable/Initialize.mat�ɂĕۑ�
%       Initialize.m�͖{�v���O�����ɂĎ��s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% �����l����т����C���[�N�X�y�[�X�ɕۑ�
run('Initialize.m');
load('./Variable/Initialize.mat');

% �p�X��ʂ�
addpath ./Tool
% �N���X�̌Ăяo��
ins_tool = tool();

% �O����
% 1.�����̓ǂݍ���
% 2.�^�̕��f�X�y�N�g���O�������擾
[music, spectrum] = ins_tool.AudioReadMethod(filename, total_sec, freq, fftsize, shiftsize, window);

% ���]�̐U���ƈʑ����擾
amp_corr = abs(spectrum);
phase_corr = angle(spectrum);

% �U���Ƀm�C�Y��������
noise_amp_corr = imnoise(amp_corr, 'speckle');


%%%%%%%%%%%%%%%%%%%%
% GLA
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start GLA \n');
% �U������ʑ��𐄒肷��A���S���Y��
spectrum_est_GLA = ins_tool.GLA(noise_amp_corr, fftsize, shiftsize, window, iteration);
% �ʑ����擾
phase_est_GLA = angle(spectrum_est_GLA);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start GLA + ADMM \n');
% �U������ʑ��𐄒肷��A���S���Y��
spectrum_est_ADMM = ins_tool.GLA_ADMM(noise_amp_corr, rho, fftsize, shiftsize, window, iteration);
% �ʑ����擾
phase_est_ADMM = angle(spectrum_est_ADMM);


%%%%%%%%%%%%%%%%%%%%
% �]��
%%%%%%%%%%%%%%%%%%%%

% ��敽�ό덷�ɂ��]���̈�
fprintf('Result :  Root Mean Square Error \n');

% ���]�̈ʑ��ƁCGLA�Ɋ�Â����肵���ʑ��ԂŐU��1�̕��f��������C��敽�ό덷
%      ���]�̈ʑ��ƐU��1�ɂ�蕡�f��������
spectrum_amp1_corr = ones( size(amp_corr) ) .* exp( 1i * phase_corr );
%      GLA�Ɋ�Â����肵���ʑ��ƐU��1�ɂ�蕡�f��������
spectrum_amp1_GLA = ones( size(amp_corr) ) .* exp( 1i * phase_est_GLA );
%      ��敽�ό덷
err_GLA = immse(spectrum_amp1_corr, spectrum_amp1_GLA);


% ���]�̈ʑ��ƁCADMM�Ɋ�Â����肵���ʑ��ԂŐU��1�̕��f��������C��敽�ό덷
%      ADMM�Ɋ�Â����肵���ʑ��ƐU��1�ɂ�蕡�f��������
spectrum_amp1_ADMM = ones( size(amp_corr) ) .* exp( 1i * phase_est_ADMM );
%      ��敽�ό덷
err_ADMM = immse(spectrum_amp1_corr, spectrum_amp1_ADMM);


% 2�敽�ό덷�̌��ʂ���
fprintf('    GLA : %d,  ADMM : %d \n', err_GLA, err_ADMM);


%%%%%%%%%%%%%%%%%%%%
% �����̏o��
%%%%%%%%%%%%%%%%%%%%

% �o�͂��邱�Ƃ���
fprintf('Output :  Sound Source \n');
% ISTFT�C���Ԏ���
signal_corr = ISTFT(spectrum, shiftsize, window)';
signal_GLA = ISTFT(spectrum_est_GLA, shiftsize, window)';
signal_ADMM = ISTFT(spectrum_est_ADMM, shiftsize, window)';


% �p�X������
rmpath ./Tool


