
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
% 6/17 40-50秒を使う
start_sec = 40;
end_sec = 50;
% 短時間フーリエ変換のフレームフレーム幅
fftsize = 1024;
% 短時間フーリエ変換のフレームシフト量
shiftsize = 256;
% ADMMのイテレーション回数を指定
iteration = 1000;
% 窓の種類
win = hann(fftsize,'periodic'); % ハニング窓
% 対象音
%filename = './Sound_source/mixture.wav';
filename = './Sound_source/01 Roundabout.mp3';
% 出力先
outputDir = './Output';
% STFTのタイプ(矢田部さんの通常のやつ=1, 理想的なSTFT=2)
STFT_type = 1;
% A特性のフィルタを作成
[~,A_weight] = filterA(ones(fftsize,1), freq);
% Douglas-Rachford Splitting Algorithm のγの値
gamma = 1.0;
% 振幅の発散を防ぐ重みweight の初期化
Delta = 0.01;
% エクセルシートの初期値
sell_angle = 3;
sell_spe = 18;


%%%%%%%%%%%%%%%%%%%%
% 前処理
%%%%%%%%%%%%%%%%%%%%

% 1.ISTFTに利用する逆の窓を合成
% 2.真の複素スペクトログラムを取得
[windual, spectrum, signal_len, music] = ins_tool.AudioReadMethod(filename, start_sec, end_sec, freq, fftsize, shiftsize, win, STFT_type);

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
phase_temp = zeros(amp_FFTsize, frames);


save('./Variable/Initialize')

