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
% �U������ʑ��𐄒肷��A���S���Y��
[spectrum_est_SDMM, min_alpha_SDMM] = ins_tool.SDMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len);
% �ʑ����擾
phase_est_SDMM= angle(spectrum_est_SDMM);


%%%%%%%%%%%%%%%%%%%%
% �]��
%%%%%%%%%%%%%%%%%%%%

% ��敽�ό덷�ɂ��]���̈�
fprintf('\n');
fprintf('Result :  Root Mean Square Error \n');

% ���]�̈ʑ��ƁC�e��@�Ɋ�Â����肵���ʑ��ԂŐU��1�̕��f��������C��敽�ό덷
%      ���]�̈ʑ��ƐU��1�ɂ�蕡�f��������i�����X�y�N�g���j
spectrum_amp1_corr = ones( size(amp_corr) ) .* exp( 1i * phase_corr );

%      �e��@�Ɋ�Â����肵���ʑ��ƐU��1�ɂ�蕡�f��������
spectrum_amp1_GLA = ones( size(amp_corr) ) .* exp( 1i * phase_est_GLA );
spectrum_amp1_ADMM = ones( size(amp_corr) ) .* exp( 1i * phase_est_ADMM );
spectrum_amp1_prop = ones( size(amp_corr) ) .* exp( 1i * phase_est_prop );
spectrum_amp1_General = ones( size(amp_corr) ) .* exp( 1i * phase_est_General );
spectrum_amp1_Douglas = ones( size(amp_corr) ) .* exp( 1i * phase_est_Douglas );
spectrum_amp1_SDMM = ones( size(amp_corr) ) .* exp( 1i * phase_est_SDMM );

%      ��敽�ό덷
err_GLA = immse(spectrum_amp1_corr, spectrum_amp1_GLA);
err_ADMM = immse(spectrum_amp1_corr, spectrum_amp1_ADMM);
err_prop = immse(spectrum_amp1_corr, spectrum_amp1_prop);
err_General = immse(spectrum_amp1_corr, spectrum_amp1_General);
err_Douglas = immse(spectrum_amp1_corr, spectrum_amp1_Douglas);
err_SDMM= immse(spectrum_amp1_corr, spectrum_amp1_SDMM);

%      ��敽�ό덷�̌��ʂ���
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, General : %d, Douglas : %d, SDMM : %d\n', err_GLA, err_ADMM, err_prop, err_General, err_Douglas, err_SDMM);


% ���z�I�ȐU���Ɛ��肵���ʑ��̕��f���ŕ��ςɂ���ĕ]��
fprintf('Result :  frobenius norm \n');

fro_GLA = norm(spectrum - spectrum_est_GLA,'fro');
fro_ADMM = norm(spectrum - spectrum_est_ADMM,'fro');
fro_prop = norm(spectrum - spectrum_est_prop,'fro');
fro_General = norm(spectrum - spectrum_est_General,'fro');
fro_Douglas = norm(spectrum - spectrum_est_Douglas,'fro');
fro_SDMM = norm(spectrum - spectrum_est_SDMM,'fro');

% ���ʂ���
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, General : %d, Douglas : %d, SDMM : %d \n', fro_GLA, fro_ADMM, fro_prop, fro_General, fro_Douglas, fro_SDMM);


%%%%%%%%%%%%%%%%%%%%
% �����̏o��
%%%%%%%%%%%%%%%%%%%%

% �o�͂��邱�Ƃ���
fprintf('Output :  Sound Source \n');
% Normalize
Normalize = @(x) x/max(abs(x));
% ISTFT�C���Ԏ���
signal_corr = ISTFT(spectrum, windual, shiftsize, fftsize, Ls);
signal_GLA = ISTFT(spectrum_est_GLA, windual, shiftsize, fftsize, Ls);
signal_ADMM = ISTFT(spectrum_est_ADMM, windual, shiftsize, fftsize, Ls);
signal_prop = ISTFT(spectrum_est_prop, windual, shiftsize, fftsize, Ls);
signal_General = ISTFT(spectrum_est_General, windual, shiftsize, fftsize, Ls);
signal_Douglas = ISTFT(spectrum_est_Douglas, windual, shiftsize, fftsize, Ls);
signal_SDMM = ISTFT(spectrum_est_Douglas, windual, shiftsize, fftsize, Ls);

% �t�H���_�쐬
[status, msg, msgID] = mkdir(sprintf('%s/signal_rho_%.2f', outputDir, rho));
% �����̏o��
audiowrite(sprintf('%s/signal_rho_%.2f/signal_correct_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_corr), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_GLA_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_GLA), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_ADMM_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_ADMM), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_prop_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_prop), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_general_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_General), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_Douglas_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_Douglas), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_SDMM_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_SDMM), freq);


%%%%%%%%%%%%%%%%%%%%
% �G�N�Z���V�[�g�ւ̏o��
%%%%%%%%%%%%%%%%%%%%

%edit( sprintf('%s/result.xlsx', outputDir) );

A = {rho, err_GLA, err_ADMM, err_prop, err_General, err_Douglas, err_SDMM, min_alpha_general, min_alpha_Douglas, min_alpha_SDMM};
xlRange = sprintf('B%d', sell_angle);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

A = {rho, fro_GLA, fro_ADMM, fro_prop, fro_General, fro_Douglas, fro_SDMM, min_alpha_general, min_alpha_Douglas, min_alpha_SDMM};
xlRange = sprintf('B%d', sell_spe);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

sell_angle = sell_angle + 1;
sell_spe = sell_spe + 1;


