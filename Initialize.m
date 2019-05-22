
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   本コードの説明
%       「Main_program.m」で使用する変数の初期値を設定
%       変数の初期値の設定に伴い，「initialize.mat」を作成
%
%   実行方法
%       コマンドウィンドウ内で「Main_program」を実行
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 標本化周波数
freq = 16000;
% 音を何秒間使うか
total_sec = 10;
% 短時間フーリエ変換のフレームフレーム幅
fftsize = 1024;
% 短時間フーリエ変換のフレームシフト量
shiftsize = 128;
% ADMMのイテレーション回数を指定
iteration = 1000;
% 窓の種類
window = 'hann';
% 対象音
filename = './Sound_source/mixture.wav';
% 出力先
outputDir = './Output';
% エクセルシートの初期値
sell_angle = 3;
sell_spe = 16;


%%%%%%%%%%%%%%%%%%%%
% 前処理
%%%%%%%%%%%%%%%%%%%%

% 1.音源の読み込み
% 2.真の複素スペクトログラムを取得
[music, spectrum] = ins_tool.AudioReadMethod(filename, total_sec, freq, fftsize, shiftsize, window);

% 所望の振幅と位相を取得
amp_corr = abs(spectrum);
phase_corr = angle(spectrum);

% 所望の振幅amp_corrから各種パラメータを取得
%       amp_FFTsize = STFT後の折り返しを考慮した，STFT/2の点数
%           例.8192点フーリエ / 2 = 4092点
%       frames = STFT/2の点数が何フレームあるか
% 位相の行列サイズを取得
[amp_FFTsize, frames] = size(amp_corr);
% ランダムな位相を取得
phase_temp = rand(amp_FFTsize, frames);


save('./Variable/Initialize')

