%clear;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 13 May. 2019.
%
%   ��������U���𒊏o�D
%   ���̐U������ʑ����Q�̎�@(GLA, GLA+ADMM, GLA+ADMM+prop, ��ʉ�ADMM)�Ő��肵�C��r����
%
%   ���s���@
%       �R�}���h�E�B���h�E���ɁuMain�v�Ŏ��s
%
%   �����l
%       Initialize.m�ɂď����l��ݒ肵�C./Variable/Initialize.mat�ɂĕۑ�
%       Initialize.m�͖{�v���O�����ɂĎ��s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
% GLA
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start GLA \n');
% �U������ʑ��𐄒肷��A���S���Y��
spectrum_est_GLA = ins_tool.GLA(amp_corr, fftsize, shiftsize, win, windual, iteration, phase_temp, Ls, signal_len);
% �ʑ����擾
phase_est_GLA = angle(spectrum_est_GLA);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start GLA + ADMM \n');
% �U������ʑ��𐄒肷��A���S���Y��
spectrum_est_ADMM = ins_tool.GLA_ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len);
% �ʑ����擾
phase_est_ADMM = angle(spectrum_est_ADMM);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + prop
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start GLA + ADMM + prop \n');
% �U������ʑ��𐄒肷��A���S���Y��
spectrum_est_prop = ins_tool.Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len);
% �ʑ����擾
phase_est_prop = angle(spectrum_est_prop);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + �o�b�`�����ɂ��prop
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start GLA + ADMM + �o�b�`�����ɂ��Prop \n');
% �U������ʑ��𐄒肷��A���S���Y��
spectrum_est_Prop_batch = ins_tool.Prop_batch(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len);
% �ʑ����擾
phase_est_Prop_batch = angle(spectrum_est_Prop_batch);


%%%%%%%%%%%%%%%%%%%%
% ��ʉ�ADMM
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start ��ʉ�ADMM \n');
% �U������ʑ��𐄒肷��A���S���Y��
[spectrum_est_General, min_alpha_general] = ins_tool.General(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len);
% �ʑ����擾
phase_est_General = angle(spectrum_est_General);


%%%%%%%%%%%%%%%%%%%%
% Douglas-Rachford Splitting Algorithm
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start Douglas-Rachford Splitting Algorithm \n');
% �U������ʑ��𐄒肷��A���S���Y��
[spectrum_est_Douglas, min_alpha_Douglas] = ins_tool.DouglasRachfordSplitting(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len, gamma);
% �ʑ����擾
phase_est_Douglas= angle(spectrum_est_Douglas);


%%%%%%%%%%%%%%%%%%%%
% SDMM
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start SDMM \n');
[spectrum_est_SDMM, min_alpha_SDMM] = ins_tool.SDMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len);
% �ʑ����擾
phase_est_SDMM= angle(spectrum_est_SDMM);


%%%%%%%%%%%%%%%%%%%%
% PPXA
%%%%%%%%%%%%%%%%%%%%

% �X�^�[�g�̈�
fprintf('Start PPXA \n');
[spectrum_est_PPXA, min_alpha_PPXA] = ins_tool.PPXA(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len, gamma);
% �ʑ����擾
phase_est_PPXA= angle(spectrum_est_PPXA);


%%%%%%%%%%%%%%%%%%%%
% �]��
%%%%%%%%%%%%%%%%%%%%

fprintf('\n');

% ���]�̈ʑ��ƁC�e��@�Ɋ�Â����肵���ʑ��ԂŐU��1�̕��f��������C��敽�ό덷
%      ���]�̈ʑ��ƐU��1�ɂ�蕡�f��������i�����X�y�N�g���j
spectrum_amp1_corr = ones( size(amp_corr) ) .* exp( 1i * phase_corr );

% �ʑ����ƃX�y�N�g�������Z�o����֐����Ăяo��
[err_GLA, fro_GLA] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_GLA);
[err_ADMM, fro_ADMM] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_ADMM);
[err_prop, fro_prop] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_prop);
[err_Prop_batch, fro_Prop_batch] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_Prop_batch);
[err_General, fro_General ] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_General);
[err_Douglas, fro_Douglas] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_Douglas);
[err_SDMM, fro_SDMM] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_SDMM);
[err_PPXA, fro_PPXA] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_PPXA);


% ��敽�ό덷�ɂ���ĕ]��
fprintf('Result :  Mean Square Error \n');
% ���ʂ���
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, Prop_batch : %d, General : %d, Douglas : %d, SDMM : %d, PPXA : %d \n', err_GLA, err_ADMM, err_prop, err_Prop_batch, err_General, err_Douglas, err_SDMM, err_PPXA);
% ���z�I�ȐU���Ɛ��肵���ʑ��̕��f���ŕ��ςɂ���ĕ]��
fprintf('Result :  frobenius norm \n');
% ���ʂ���
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, Prop_batch : %d, General : %d, Douglas : %d, SDMM : %d, PPXA : %d \n', fro_GLA, fro_ADMM, fro_prop, fro_Prop_batch, fro_General, fro_Douglas, fro_SDMM, fro_PPXA);


%%%%%%%%%%%%%%%%%%%%
% �����̏o��
%%%%%%%%%%%%%%%%%%%%

% �o�͂��邱�Ƃ���
fprintf('Output :  Sound Source \n');

% �t�H���_�쐬
[status, msg, msgID] = mkdir(sprintf('%s/signal_rho_%.2f', outputDir, rho));

ins_tool.OutputMethod(spectrum, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'correct');
ins_tool.OutputMethod(spectrum_est_GLA, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'GLA');
ins_tool.OutputMethod(spectrum_est_ADMM, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'ADMM');
ins_tool.OutputMethod(spectrum_est_prop, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'prop');
ins_tool.OutputMethod(spectrum_est_Prop_batch, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'Prop_batch');
ins_tool.OutputMethod(spectrum_est_General, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'General');
ins_tool.OutputMethod(spectrum_est_Douglas, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'Douglas');
ins_tool.OutputMethod(spectrum_est_PPXA, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'PPXA');


%%%%%%%%%%%%%%%%%%%%
% �G�N�Z���V�[�g�ւ̏o��
%%%%%%%%%%%%%%%%%%%%

A = {rho, err_GLA, err_ADMM, err_prop, err_Prop_batch, err_General, err_Douglas, err_SDMM, err_PPXA, min_alpha_general, min_alpha_Douglas, min_alpha_SDMM, min_alpha_PPXA};
xlRange = sprintf('B%d', sell_angle);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

A = {rho, fro_GLA, fro_ADMM, fro_prop, fro_Prop_batch, fro_General, fro_Douglas, fro_SDMM, fro_PPXA, min_alpha_general, min_alpha_Douglas, min_alpha_SDMM, min_alpha_PPXA};
xlRange = sprintf('B%d', sell_spe);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

sell_angle = sell_angle + 1;
sell_spe = sell_spe + 1;


