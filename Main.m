
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Corded by R.Nakatsu (dragonstar30210326@gmail.com) on 13 May. 2019.
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
addpath ./OptimizeMethod
% クラスの呼び出し
ins_tool = tool();

% 初期値をよびだし，ワークスペースに保存
run('Initialize.m');
load('./Variable/Initialize.mat');

%一つ目の実験を開始
fprintf('**********Experiment clear amp**********\n');

% admmのパラメータρ ( ρ = 0.1, 0.2, 10, 100)
%for rho = 0.1:0.1:1.0
%for rho = [0.001, 0.01, 0.1, 10, 100]
for rho = [0.001, 0.01]
    
    fprintf('\n');
    fprintf('rho = %d \n', rho);
    fprintf('\n');
        
    run('Experiment_clear_amp.m');
    save(sprintf('./Variable/result_rho_%.2f.mat', rho));
    
end
    
%ノイズありの振幅スペクトログラムに対する実験を開始
%fprintf('**********Experiment noise amp**********\n');
%run('Experiment_noise_amp.m');

%音を聞くとき
%sound(y, freq);
%sound(true_sound, freq);
%sound(GLA_sound, freq);
%sound(ADMM_sound, freq);
%sound(Prop_sound, freq);
%sound(gADMM_sound, freq);

%グラフ
%{
figure(200)
plot(gADMM_A_error(100:end))
legend('GLA','ADMM','Prop','PropGeneral');
xlim([0 900]);
%}

%グラフのフォントサイズは13くらい


