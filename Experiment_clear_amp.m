%clear;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 13 May. 2019.
%
%   音源から振幅を抽出．
%   その振幅から位相を２つの手法(GLA, GLA+ADMM, GLA+ADMM+prop, 一般化ADMM)で推定し，比較する
%
%   実行方法
%       コマンドウィンドウ内に「Main」で実行
%
%   初期値
%       Initialize.mにて初期値を設定し，./Variable/Initialize.matにて保存
%       Initialize.mは本プログラムにて実行
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
% GLA
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA \n');
% 振幅から位相を推定するアルゴリズム
%[GLA_spe, GLA_sound, GLA_amp_err, GLA_A_err] = ins_tool.GLA(amp_corr, fftsize, shiftsize, win, windual, iteration, phase_temp, signal_len, STFT_type, A_weight);
[GLA_spe, GLA_sound, GLA_amp_err, GLA_A_err] = GLA(amp_corr, fftsize, shiftsize, win, windual, iteration, phase_temp, signal_len, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% ADMM （矢田部法）
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM \n');
% 振幅から位相を推定するアルゴリズム
%[ADMM_spe, ADMM_sound, ADMM_amp_err, ADMM_A_err] = ins_tool.ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, STFT_type, A_weight);
[ADMM_spe, ADMM_sound, ADMM_amp_err, ADMM_A_err] = ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + prop
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + prop \n');
% 振幅から位相を推定するアルゴリズム
%[prop_spe, prop_sound, prop_amp_err, prop_A_err] = ins_tool.Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, STFT_type, A_weight);
[prop_spe, prop_sound, prop_amp_err, prop_A_err] = Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + バッチ処理によるprop
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + バッチ処理によるProp \n');
% 振幅から位相を推定するアルゴリズム
%[prop_batch_spe, prop_batch_sound, prop_batch_amp_err, prop_batch_A_err] = ins_tool.Prop_batch(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, STFT_type, A_weight);
[prop_batch_spe, prop_batch_sound, prop_batch_amp_err, prop_batch_A_err] = Prop_batch(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + prop + 振幅weight
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + prop + 振幅weight \n');
% 振幅から位相を推定するアルゴリズム
%[prop_weight_spe, prop_weight_sound, prop_weight_amp_err, prop_weight_A_err] = ins_tool.PropWeight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, Delta, STFT_type, A_weight);
[prop_weight_spe, prop_weight_sound, prop_weight_amp_err, prop_weight_A_err] = PropWeight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, Delta, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + バッチ処理によるprop + 振幅weight
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + バッチ処理によるProp + 振幅weight \n');
% 振幅から位相を推定するアルゴリズム
%[prop_batch_weight_spe, prop_batch_weight_sound, prop_batch_weight_amp_err, prop_batch_weight_A_err] = ins_tool.Prop_batch_weight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, Delta, STFT_type, A_weight);
[prop_batch_weight_spe, prop_batch_weight_sound, prop_batch_weight_amp_err, prop_batch_weight_A_err] = Prop_batch_weight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, signal_len, Delta, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% 一般化ADMM
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start 一般化ADMM \n');
% 振幅から位相を推定するアルゴリズム
%[Gprop_spe, Gprop_sound, Gprop_amp_err, Gprop_A_err, min_alpha_Gprop] = ins_tool.General(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, signal_len, STFT_type, A_weight);
[Gprop_spe, Gprop_sound, Gprop_amp_err, Gprop_A_err, min_alpha_Gprop] = General(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, signal_len, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% Douglas-Rachford Splitting Algorithm
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start Douglas-Rachford Splitting Algorithm \n');
% 振幅から位相を推定するアルゴリズム
%[Douglas_spe, Douglas_sound, Douglas_amp_err, Douglas_A_err, min_alpha_Douglas] = ins_tool.DouglasRachfordSplitting(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, signal_len, gamma, STFT_type, A_weight);
[Douglas_spe, Douglas_sound, Douglas_amp_err, Douglas_A_err, min_alpha_Douglas] = DouglasRachfordSplitting(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, signal_len, gamma, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% SDMM
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start SDMM \n');
%[SDMM_spe, SDMM_sound, SDMM_amp_err, SDMM_A_err, min_alpha_SDMM] = ins_tool.SDMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, signal_len, STFT_type, A_weight);
[SDMM_spe, SDMM_sound, SDMM_amp_err, SDMM_A_err, min_alpha_SDMM] = SDMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, signal_len, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% PPXA
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start PPXA \n');
%[PPXA_spe, PPXA_sound, PPXA_amp_err, PPXA_A_err, min_alpha_PPXA] = ins_tool.PPXA(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, signal_len, gamma, STFT_type, A_weight);
[PPXA_spe, PPXA_sound, PPXA_amp_err, PPXA_A_err, min_alpha_PPXA] = PPXA(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, signal_len, gamma, STFT_type, A_weight);


%%%%%%%%%%%%%%%%%%%%
% 評価
%%%%%%%%%%%%%%%%%%%%

fprintf('\n');

%{

% 所望の位相と，各手法に基づき推定した位相間で振幅1の複素数を仮定，二乗平均誤差
%      所望の位相と振幅1により複素数を仮定（正解スペクトル）
spectrum_amp1_corr = ones( size(amp_corr) ) .* exp( 1i * phase_corr );

% 位相差とスペクトル差を算出する関数を呼び出し

[err_GLA, fro_GLA] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, GLA_spe);
[err_ADMM, fro_ADMM] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_ADMM);
[err_prop, fro_prop] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_prop);
[err_Prop_batch, fro_Prop_batch] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_Prop_batch);
[err_General, fro_General ] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_General);
[err_Douglas, fro_Douglas] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_Douglas);
[err_SDMM, fro_SDMM] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_SDMM);
[err_PPXA, fro_PPXA] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_PPXA);
[err_prop_weight, fro_prop_weight] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_prop_weight);
[err_Prop_batch_weight, fro_Prop_batch_weight] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_Prop_batch_weight);

%}

%{

% 二乗平均誤差によって評価
fprintf('Result :  Mean Square Error \n');
% 結果を印字
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, Prop_batch : %d, General : %d, Douglas : %d, SDMM : %d, PPXA : %d, prop_weight : %d, Prop_batch_weight : %d\n', err_GLA, err_ADMM, err_prop, err_Prop_batch, err_General, err_Douglas, err_SDMM, err_PPXA, err_prop_weight, err_Prop_batch_weight);
% 理想的な振幅と推定した位相の複素数で平均によって評価
fprintf('Result :  frobenius norm \n');
% 結果を印字
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, Prop_batch : %d, General : %d, Douglas : %d, SDMM : %d, PPXA : %d, prop_weight : %d, Prop_batch_weight : %d\n', fro_GLA, fro_ADMM, fro_prop, fro_Prop_batch, fro_General, fro_Douglas, fro_SDMM, fro_PPXA, fro_prop_weight, fro_Prop_batch_weight);

%}

%{
% 振幅の誤差によって評価
fprintf('Result : Amp Error \n');
% 結果を印字
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, Prop_batch : %d, Prop_weight : %d, Prop_batch_weight : %d, General : %d, Douglas : %d, SDMM : %d, PPXA : %d\n', GLA_amp_err(end), ADMM_amp_err(end), prop_amp_err(end), prop_batch_amp_err(end), prop_weight_amp_err(end), prop_batch_weight_amp_err(end), Gprop_amp_err(end), Douglas_amp_err(end), SDMM_amp_err(end), PPXA_amp_err(end));
% A特性をかけた振幅によって評価
fprintf('Result : A filter Amp Error \n');
% 結果を印字
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, Prop_batch : %d, Prop_weight : %d, Prop_batch_weight : %d, General : %d, Douglas : %d, SDMM : %d, PPXA : %d\n', GLA_A_err(end), ADMM_A_err(end), prop_A_err(end), prop_batch_A_err(end), prop_weight_A_err(end), prop_batch_weight_A_err(end), Gprop_A_err(end), Douglas_A_err(end), SDMM_A_err(end), PPXA_A_err(end));
%}

% 振幅の誤差によって評価
fprintf('Result : Error \n');

ins_tool.OutputPrintingMethod('GLA', GLA_amp_err(end), GLA_A_err(end));
ins_tool.OutputPrintingMethod('ADMM', ADMM_amp_err(end), ADMM_A_err(end));
ins_tool.OutputPrintingMethod('prop', prop_amp_err(end), prop_A_err(end));
ins_tool.OutputPrintingMethod('prop_batch', prop_batch_amp_err(end), prop_batch_A_err(end));
ins_tool.OutputPrintingMethod('prop_weight', prop_weight_amp_err(end), prop_weight_A_err(end));
ins_tool.OutputPrintingMethod('prop_batch_weight', prop_batch_weight_amp_err(end), prop_batch_weight_A_err(end));
ins_tool.OutputPrintingMethod('Genaral', Gprop_amp_err(end), Gprop_A_err(end));
ins_tool.OutputPrintingMethod('Douglas', Douglas_amp_err(end), Douglas_A_err(end));
ins_tool.OutputPrintingMethod('SDMM', SDMM_amp_err(end), SDMM_A_err(end));
ins_tool.OutputPrintingMethod('PPXA', PPXA_amp_err(end), PPXA_A_err(end));

fprintf('\n');

%%%%%%%%%%%%%%%%%%%%
% 音源の出力
%%%%%%%%%%%%%%%%%%%%

% 出力することを印字
fprintf('Output :  Sound Source \n');

% フォルダ作成
[status, msg, msgID] = mkdir(sprintf('%s/signal_rho_%.2f', outputDir, rho));

% 音源をwave形式で出力
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'correct', rho), music, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'GLA', rho), GLA_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'ADMM', rho), ADMM_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'prop', rho), prop_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'prop_batch', rho), prop_batch_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'prop_weight', rho), prop_weight_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'prop_batch_weight', rho), prop_batch_weight_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'general', rho), Gprop_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'Douglas', rho), Douglas_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'SDMM', rho), SDMM_sound, freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, 'PPXA', rho), PPXA_sound, freq);

%{
ins_tool.OutputMethod(spectrum, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'correct');
ins_tool.OutputMethod(GLA_spe, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'GLA');
ins_tool.OutputMethod(spectrum_est_ADMM, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'ADMM');
ins_tool.OutputMethod(spectrum_est_prop, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'prop');
ins_tool.OutputMethod(spectrum_est_Prop_batch, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'Prop_batch');
ins_tool.OutputMethod(spectrum_est_General, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'General');
ins_tool.OutputMethod(spectrum_est_Douglas, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'Douglas');
ins_tool.OutputMethod(spectrum_est_PPXA, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'PPXA');
ins_tool.OutputMethod(spectrum_est_prop_weight, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'prop_weight');
ins_tool.OutputMethod(spectrum_est_Prop_batch_weight, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, 'Prop_batch_weight');
%}

%%%%%%%%%%%%%%%%%%%%
% エクセルシートへの出力
%%%%%%%%%%%%%%%%%%%%

A = {rho, GLA_amp_err(end), ADMM_amp_err(end), prop_amp_err(end), prop_batch_amp_err(end), prop_weight_amp_err(end), prop_batch_weight_amp_err(end), Gprop_amp_err(end), Douglas_amp_err(end), SDMM_amp_err(end), PPXA_amp_err(end), min_alpha_Gprop, min_alpha_Douglas, min_alpha_SDMM, min_alpha_PPXA};
xlRange = sprintf('B%d', sell_angle);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

A = {rho, GLA_A_err(end), ADMM_A_err(end), prop_A_err(end), prop_batch_A_err(end), prop_weight_A_err(end), prop_batch_weight_A_err(end), Gprop_A_err(end), Douglas_A_err(end), SDMM_A_err(end), PPXA_A_err(end), min_alpha_Gprop, min_alpha_Douglas, min_alpha_SDMM, min_alpha_PPXA};
xlRange = sprintf('B%d', sell_spe);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

sell_angle = sell_angle + 1;
sell_spe = sell_spe + 1;


