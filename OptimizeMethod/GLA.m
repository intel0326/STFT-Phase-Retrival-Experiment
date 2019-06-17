function [spectrum, sound, amp_error, A_error] = GLA(amp_corr, fftsize, shiftsize, win, windual, iteration, phase_temp, signal_len, STFT_type, A_weight)
   %
   % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 18 Apr. 2019.
   %
   % [inputs]
   %   amp_corr: correct amplitude ( FrequencyBin * Frames)
   %   fftsize: FFT length (=8192)
   %   shiftsize: frame shift (default: fftSize/2)
   %   window: window function used in STFT (fftSize x 1) or choose used
   %   iteration: iteration number
   %   phase_temp: Random phase used as initial value
   %
   % [outputs]
   %   spectrum: frequency-domain (fftSize/2+1 x frames)
   %
   
   %%%%%%%%%%%%%%%%%%%%
   % GLA
   %%%%%%%%%%%%%%%%%%%%
   
   % 初期値
   %       spectrum = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
   %       amp_error = 0 : 振幅間のフロベニウスノルム
   %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
   spectrum = amp_corr .* exp(1i * phase_temp);
   amp_error = zeros(iteration,1);
   A_error = zeros(iteration,1);
   
   
   % イテレーションの最後を値域空間への射影で終わらしたいので，
   % イテレーション1回目はループの外側で
   sound = ISTFT(spectrum, windual, shiftsize, fftsize, signal_len, STFT_type);
   spectrum = STFT(sound, win, shiftsize, fftsize, STFT_type);
   
   % 最初の誤差を出力
   amp_error(1) = norm(abs(spectrum)-amp_corr,'fro');
   A_error(1) = norm(diag(A_weight)*(abs(spectrum)-amp_corr),'fro');
   
   
   for i = 2:iteration
       
       % イテレーション回数の印字
       %fprintf('    Iteration : %d \n', i);
       
       % 振幅集合への射影，ampは所望に変更，位相だけ保管
       % アダマール積に注意
       spectrum = amp_corr .* exp( 1i * angle(spectrum) );
       
       % 位相を更新
       sound = ISTFT(spectrum, windual, shiftsize, fftsize, signal_len, STFT_type);
       spectrum = STFT(sound, win, shiftsize, fftsize, STFT_type);
       
       % 誤差を出力
       amp_error(i) = norm(abs(spectrum)-amp_corr,'fro');
       A_error(i) = norm(diag(A_weight)*(abs(spectrum)-amp_corr),'fro');
       
   end
   
end