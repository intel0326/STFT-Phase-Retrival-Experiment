function [spectrum, sound, amp_error, A_error] = ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, signal_len, STFT_type, A_weight)
    %
    % Corded by R.Nakatsu (dragonstar30210326@gmail.com) on 16 Apr. 2019.
    %
    % [inputs]
    %   amp_corr: correct amplitude ( FrequencyBin * Frames)
    %   rho: ADMM Parameter ( ρ = 0.1, 0.2, 10, 100)
    %   fftsize: FFT length (=8192)
    %   shiftsize: frame shift (default: fftSize/2)
    %   window: window function used in STFT (fftSize x 1) or choose used
    %   iteration: iteration number
    %   phase_temp: Random phase used as initial value
    %   freq: fft length/2 (8192/2)
    %
    % [outputs]
    %   spectrum: frequency-domain (fftSize/2+1 x frames)
    %
    
    %%%%%%%%%%%%%%%%%%%%
    % ADMM (升山法)
    %%%%%%%%%%%%%%%%%%%%
    
    % 初期値
    %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
    %       z =  x: ADMMで解く最適化問題に落とし込むため，新たに変数zを用意
    %       u = 0 : ラグランジュ未定乗数
    %       amp_error = 0 : 振幅間のフロベニウスノルム
    %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
    x = amp_corr .* exp(1i * phase_temp);
    z = x;
    u = zeros(freq, frames);
    amp_error = zeros(iteration,1);
    A_error = zeros(iteration,1);
    
    for i = 1:iteration
        % イテレーション回数の印字
        %fprintf('    Iteration : %d \n', i);
        
        % スペクトログラムの更新，ampは所望に変更，位相はだけ保管
        % アダマール積に注意
        x = amp_corr .* exp( 1i * angle(z - u) );
        
        % 新たな変数であるzを更新
        %   下準備として STFT( ISTFT() )をおこなう
        temp_sound = ISTFT(x+u, windual, shiftsize, fftsize, signal_len, STFT_type);
        temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
        %   zの更新
        z = ( rho*temp_spectrum + x + u ) / (1 + rho);
        
        % ラグランジュ未定乗数uの更新
        u = u + x - z;
        
        % イテレーション毎に評価するために一度STFTの値域空間に射影
        sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
        spectrum = STFT(sound, win, shiftsize, fftsize, STFT_type);
        
        % 誤差を出力
        amp_error(i) = norm(abs(spectrum)-amp_corr,'fro');
        A_error(i) = norm(diag(A_weight)*(abs(spectrum)-amp_corr),'fro');
        
    end
    
end