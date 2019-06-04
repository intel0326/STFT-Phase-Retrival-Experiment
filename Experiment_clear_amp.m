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
spectrum_est_GLA = ins_tool.GLA(amp_corr, fftsize, shiftsize, win, windual, iteration, phase_temp, Ls, signal_len);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_ADMM = ins_tool.GLA_ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + prop
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + prop \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_prop = ins_tool.Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + バッチ処理によるprop
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + バッチ処理によるProp \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_Prop_batch = ins_tool.Prop_batch(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + prop + 振幅weight
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + prop + 振幅weight \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_prop_weight = ins_tool.PropWeight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len, Delta);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + バッチ処理によるprop + 振幅weight
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + バッチ処理によるProp + 振幅weight \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_Prop_batch_weight = ins_tool.Prop_batch_weight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len, Delta);


%%%%%%%%%%%%%%%%%%%%
% 一般化ADMM
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start 一般化ADMM \n');
% 振幅から位相を推定するアルゴリズム
[spectrum_est_General, min_alpha_general] = ins_tool.General(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len);


%%%%%%%%%%%%%%%%%%%%
% Douglas-Rachford Splitting Algorithm
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start Douglas-Rachford Splitting Algorithm \n');
% 振幅から位相を推定するアルゴリズム
[spectrum_est_Douglas, min_alpha_Douglas] = ins_tool.DouglasRachfordSplitting(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len, gamma);


%%%%%%%%%%%%%%%%%%%%
% SDMM
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start SDMM \n');
[spectrum_est_SDMM, min_alpha_SDMM] = ins_tool.SDMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len);


%%%%%%%%%%%%%%%%%%%%
% PPXA
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start PPXA \n');
[spectrum_est_PPXA, min_alpha_PPXA] = ins_tool.PPXA(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len, gamma);


%%%%%%%%%%%%%%%%%%%%
% 評価
%%%%%%%%%%%%%%%%%%%%

fprintf('\n');

% 所望の位相と，各手法に基づき推定した位相間で振幅1の複素数を仮定，二乗平均誤差
%      所望の位相と振幅1により複素数を仮定（正解スペクトル）
spectrum_amp1_corr = ones( size(amp_corr) ) .* exp( 1i * phase_corr );

% 位相差とスペクトル差を算出する関数を呼び出し
[err_GLA, fro_GLA] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_GLA);
[err_ADMM, fro_ADMM] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_ADMM);
[err_prop, fro_prop] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_prop);
[err_Prop_batch, fro_Prop_batch] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_Prop_batch);
[err_General, fro_General ] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_General);
[err_Douglas, fro_Douglas] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_Douglas);
[err_SDMM, fro_SDMM] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_SDMM);
[err_PPXA, fro_PPXA] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_PPXA);
[err_prop_weight, fro_prop_weight] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_prop_weight);
[err_Prop_batch_weight, fro_Prop_batch_weight] = ins_tool.evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est_Prop_batch_weight);

% 二乗平均誤差によって評価
fprintf('Result :  Mean Square Error \n');
% 結果を印字
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, Prop_batch : %d, General : %d, Douglas : %d, SDMM : %d, PPXA : %d, prop_weight : %d, Prop_batch_weight : %d\n', err_GLA, err_ADMM, err_prop, err_Prop_batch, err_General, err_Douglas, err_SDMM, err_PPXA, err_prop_weight, err_Prop_batch_weight);
% 理想的な振幅と推定した位相の複素数で平均によって評価
fprintf('Result :  frobenius norm \n');
% 結果を印字
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, Prop_batch : %d, General : %d, Douglas : %d, SDMM : %d, PPXA : %d, prop_weight : %d, Prop_batch_weight : %d\n', fro_GLA, fro_ADMM, fro_prop, fro_Prop_batch, fro_General, fro_Douglas, fro_SDMM, fro_PPXA, fro_prop_weight, fro_Prop_batch_weight);


%%%%%%%%%%%%%%%%%%%%
% 音源の出力
%%%%%%%%%%%%%%%%%%%%

% 出力することを印字
fprintf('Output :  Sound Source \n');

% フォルダ作成
[status, msg, msgID] = mkdir(sprintf('%s/signal_rho_%.2f', outputDir, rho));

ins_tool.OutputMethod(spectrum, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'correct');
ins_tool.OutputMethod(spectrum_est_GLA, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'GLA');
ins_tool.OutputMethod(spectrum_est_ADMM, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'ADMM');
ins_tool.OutputMethod(spectrum_est_prop, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'prop');
ins_tool.OutputMethod(spectrum_est_Prop_batch, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'Prop_batch');
ins_tool.OutputMethod(spectrum_est_General, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'General');
ins_tool.OutputMethod(spectrum_est_Douglas, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'Douglas');
ins_tool.OutputMethod(spectrum_est_PPXA, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'PPXA');
ins_tool.OutputMethod(spectrum_est_prop_weight, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'prop_weight');
ins_tool.OutputMethod(spectrum_est_Prop_batch_weight, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, 'Prop_batch_weight');


%%%%%%%%%%%%%%%%%%%%
% エクセルシートへの出力
%%%%%%%%%%%%%%%%%%%%

A = {rho, err_GLA, err_ADMM, err_prop, err_Prop_batch, err_General, err_Douglas, err_SDMM, err_PPXA, err_prop_weight, err_Prop_batch_weight, min_alpha_general, min_alpha_Douglas, min_alpha_SDMM, min_alpha_PPXA};
xlRange = sprintf('B%d', sell_angle);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

A = {rho, fro_GLA, fro_ADMM, fro_prop, fro_Prop_batch, fro_General, fro_Douglas, fro_SDMM, fro_PPXA, fro_prop_weight, fro_Prop_batch_weight, min_alpha_general, min_alpha_Douglas, min_alpha_SDMM, min_alpha_PPXA};
xlRange = sprintf('B%d', sell_spe);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

sell_angle = sell_angle + 1;
sell_spe = sell_spe + 1;


