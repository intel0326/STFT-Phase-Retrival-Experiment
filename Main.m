
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 13 May. 2019.
%
%   実験内容
%        音源から振幅を抽出．
%        その振幅から位相を２つの手法(GLA, GLA+ADMM, GLA+ADMM+prop, 一般化ADMM)で推定し，比較する
%
%        音源から振幅を抽出．
%        振幅にノイズを加えてから位相を２つの手法(GLA, GLA+ADMM, GLA+ADMM+prop, 一般化ADMM)で推定し，比較する
%
%   実行方法
%       コマンドウィンドウ内に「Main」で実行
%       
%
%   初期値
%       Initialize.mにて初期値を設定し，./Variable/Initialize.matにて保存
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% パスを通す
addpath ./Tool
% クラスの呼び出し
ins_tool = tool();

% 初期値をよびだし，ワークスペースに保存
run('Initialize.m');
load('./Variable/Initialize.mat');

%一つ目の実験を開始
fprintf('**********Experiment clear amp**********\n');

% admmのパラメータρ ( ρ = 0.1, 0.2, 10, 100)
%for rho = 0.1:0.1:1.0
%for rho = [0.01, 0.1, 0.2, 10, 20, 100]
for rho = [0.1, 10]
    
    % 10回のイテレーション内で位相とスペクトルを保存
    [line, column]=size(spectrum);
    result_spe = zeros(line, column, 10);
    result_angle = zeros(line, column, 10);
    
    Mean_spe = zeros(line, column, 2);
    Var_spe = zeros(line, column, 2);
    
    Mean_angle = zeros(line, column, 2);
    Var_angle = zeros(line, column, 2);
    
    %カウント
    count=1;
    
    % 10回のイテレーションで初期位相に差がでるか検証
    for q = 1:1:10
        fprintf('\n');
        fprintf('q = %d, rho = %d \n', q, rho);
        fprintf('\n');
        
        % ランダムな位相を取得
        phase_temp = rand(amp_FFTsize, frames);
        
        run('Experiment_clear_amp.m');
        save(sprintf('./Variable/result_rho_%.2f.mat', rho));
    end
    
    % 平均と分散を出力
    Mean_spe(:, :, count) = mean(result_spe, 3);
    Var_spe(:, :, count) = var(result_spe, 0, 3);
    
    Mean_angle(:, :, count) = mean(result_angle, 3);
    Var_angle(:, :, count) = var(result_angle, 0, 3);
    
    count = count + 1;
    
    sell_angle = sell_angle + 1;
    sell_spe = sell_spe + 1;
    
end
    
%二つ目の実験を開始
%fprintf('**********Experiment noise amp**********\n');
%run('Experiment_noise_amp.m');

% パスを消す
rmpath ./Tool


