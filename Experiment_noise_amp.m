clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 22 Apr. 2019.
%
%   音源から振幅を抽出．
%   振幅にノイズを加えてから位相を２つの手法(GLA, GLA+ADMM)で推定し，比較する
%
%   実行方法
%       コマンドウィンドウ内に「Main」で実行
%
%   初期値
%       Initialize.mにて初期値を設定し，./Variable/Initialize.matにて保存
%       Initialize.mは本プログラムにて実行
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 初期値をよびだし，ワークスペースに保存
run('Initialize.m');
load('./Variable/Initialize.mat');

% パスを通す
addpath ./Tool
% クラスの呼び出し
ins_tool = tool();

% 前処理
% 1.音源の読み込み
% 2.真の複素スペクトログラムを取得
[music, spectrum] = ins_tool.AudioReadMethod(filename, total_sec, freq, fftsize, shiftsize, window);

% 所望の振幅と位相を取得
amp_corr = abs(spectrum);
phase_corr = angle(spectrum);

% 振幅にノイズを加える
noise_amp_corr = imnoise(amp_corr, 'speckle');


%%%%%%%%%%%%%%%%%%%%
% GLA
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_GLA = ins_tool.GLA(noise_amp_corr, fftsize, shiftsize, window, iteration);
% 位相を取得
phase_est_GLA = angle(spectrum_est_GLA);


%%%%%%%%%%%%%%%%%%%%
% GLA + ADMM
%%%%%%%%%%%%%%%%%%%%

% スタートの印字
fprintf('Start GLA + ADMM \n');
% 振幅から位相を推定するアルゴリズム
spectrum_est_ADMM = ins_tool.GLA_ADMM(noise_amp_corr, rho, fftsize, shiftsize, window, iteration);
% 位相を取得
phase_est_ADMM = angle(spectrum_est_ADMM);


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


% 2乗平均誤差の結果を印字
fprintf('    GLA : %d,  ADMM : %d \n', err_GLA, err_ADMM);


%%%%%%%%%%%%%%%%%%%%
% 音源の出力
%%%%%%%%%%%%%%%%%%%%

% 出力することを印字
fprintf('Output :  Sound Source \n');
% ISTFT，時間軸に
signal_corr = ISTFT(spectrum, shiftsize, window)';
signal_GLA = ISTFT(spectrum_est_GLA, shiftsize, window)';
signal_ADMM = ISTFT(spectrum_est_ADMM, shiftsize, window)';


% パスを消す
rmpath ./Tool


