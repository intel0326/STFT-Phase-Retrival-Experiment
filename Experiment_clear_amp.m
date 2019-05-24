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
% 位相を取得
phase_est_GLA = angle(spectrum_est_GLA);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_ADMM = ins_tool.GLA_ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len);
% 位相を取得
phase_est_ADMM = angle(spectrum_est_ADMM);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM + prop
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM + prop \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_prop = ins_tool.Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, frames, Ls, signal_len);
% 位相を取得
phase_est_prop = angle(spectrum_est_prop);


%%%%%%%%%%%%%%%%%%%%
% 一般化ADMM
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start 一般化ADMM \n');
% 振幅から位相を推定するアルゴリズム
[spectrum_est_General, min_alpha] = ins_tool.General(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, amp_FFTsize, spectrum, frames, Ls, signal_len);
% 位相を取得
phase_est_General = angle(spectrum_est_General);


%%%%%%%%%%%%%%%%%%%%
% 評価
%%%%%%%%%%%%%%%%%%%%

% 二乗平均誤差による評価の印字
fprintf('Result :  Root Mean Square Error \n');

% 所望の位相と，GLAに基づき推定した位相間で振幅1の複素数を仮定，二乗平均誤差
%      所望の位相と振幅1により複素数を仮定
spectrum_amp1_corr = ones( size(amp_corr) ) .* exp( 1i * phase_corr );
%      GLAに基づき推定した位相と振幅1により複素数を仮定
spectrum_amp1_GLA = ones( size(amp_corr) ) .* exp( 1i * phase_est_GLA );
%      二乗平均誤差
err_GLA = immse(spectrum_amp1_corr, spectrum_amp1_GLA);

% 所望の位相と，ADMMに基づき推定した位相間で振幅1の複素数を仮定，二乗平均誤差
%      ADMMに基づき推定した位相と振幅1により複素数を仮定
spectrum_amp1_ADMM = ones( size(amp_corr) ) .* exp( 1i * phase_est_ADMM );
%      二乗平均誤差
err_ADMM = immse(spectrum_amp1_corr, spectrum_amp1_ADMM);

% 所望の位相と，Propに基づき推定した位相間で振幅1の複素数を仮定，二乗平均誤差
%      Propに基づき推定した位相と振幅1により複素数を仮定
spectrum_amp1_prop = ones( size(amp_corr) ) .* exp( 1i * phase_est_prop );
%      二乗平均誤差
err_prop = immse(spectrum_amp1_corr, spectrum_amp1_prop);

% 所望の位相と，generalに基づき推定した位相間で振幅1の複素数を仮定，二乗平均誤差
%      generalに基づき推定した位相と振幅1により複素数を仮定
spectrum_amp1_General = ones( size(amp_corr) ) .* exp( 1i * phase_est_General );
%      二乗平均誤差
err_General = immse(spectrum_amp1_corr, spectrum_amp1_General);


% 2乗平均誤差の結果を印字
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, General : %d \n', err_GLA, err_ADMM, err_prop, err_General);

% 理想的な振幅と推定した位相の複素数で平均によって評価
fprintf('Result :  frobenius norm \n');

fro_GLA = norm(spectrum - spectrum_est_GLA,'fro');
fro_ADMM = norm(spectrum - spectrum_est_ADMM,'fro');
fro_prop = norm(spectrum - spectrum_est_prop,'fro');
fro_General = norm(spectrum - spectrum_est_General,'fro');


% 結果を印字
fprintf('    GLA : %d,  ADMM : %d,  Prop : %d, General : %d \n', fro_GLA, fro_ADMM, fro_prop, fro_General);


%%%%%%%%%%%%%%%%%%%%
% 音源の出力
%%%%%%%%%%%%%%%%%%%%

% 出力することを印字
fprintf('Output :  Sound Source \n');
% Normalize
Normalize = @(x) x/max(abs(x));
% ISTFT，時間軸に
signal_corr = ISTFT(spectrum, windual, shiftsize, fftsize, Ls);
signal_GLA = ISTFT(spectrum_est_GLA, windual, shiftsize, fftsize, Ls);
signal_ADMM = ISTFT(spectrum_est_ADMM, windual, shiftsize, fftsize, Ls);
signal_prop = ISTFT(spectrum_est_prop, windual, shiftsize, fftsize, Ls);
signal_General = ISTFT(spectrum_est_General, windual, shiftsize, fftsize, Ls);
% フォルダ作成
[status, msg, msgID] = mkdir(sprintf('%s/signal_rho_%.2f', outputDir, rho));
% 音源の出力
audiowrite(sprintf('%s/signal_rho_%.2f/signal_correct_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_corr), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_GLA_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_GLA), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_ADMM_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_ADMM), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_prop_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_prop), freq);
audiowrite(sprintf('%s/signal_rho_%.2f/signal_general_rho=%.2f.wav', outputDir, rho, rho), Normalize(signal_General), freq);


%%%%%%%%%%%%%%%%%%%%
% エクセルシートへの出力
%%%%%%%%%%%%%%%%%%%%

%edit( sprintf('%s/result.xlsx', outputDir) );

A = {rho, err_GLA, err_ADMM, err_prop, err_General, min_alpha};
xlRange = sprintf('B%d', sell_angle);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

A = {rho, fro_GLA, fro_ADMM, fro_prop, fro_General, min_alpha};
xlRange = sprintf('B%d', sell_spe);
xlswrite(sprintf('%s/result.xlsx', outputDir), A, 1, xlRange);

sell_angle = sell_angle + 1;
sell_spe = sell_spe + 1;


