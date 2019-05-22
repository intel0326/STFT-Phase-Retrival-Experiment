classdef tool
    methods (Static)
        function [windual, Spe, Ls] = AudioReadMethod(filename, total_sec, freq, fftsize, shiftsize, win)
            %サンプリングレートのみ取得
            %いらない音源の部分は「~」で消しておく
            [~,fs] = audioread(filename);
            %読み込む長さを計算
            samples = [1, total_sec*fs];
            %音源取得
            [music, fs] = audioread(filename, samples);
            %ダウンサンプリング
            music = resample(music, freq, fs);
            %逆の窓を合成
            windual = winDual(win, shiftsize);
            % !! Ls must be even number due to our STFT/iSTFT implementation !!
            Ls = ceil((length(music)+2*(fftsize-shiftsize)-fftsize)/shiftsize)*shiftsize+fftsize;
            % zero パディング：信号の両端を0詰め
            music = [ zeros(fftsize-shiftsize,1); music; zeros(Ls-length(music)-2*(fftsize-shiftsize),1); zeros(fftsize-shiftsize,1) ];
            % STFT
            Spe = STFT(music, win, shiftsize, fftsize, Ls);
        end
        
        function spectrum = GLA(amp_corr, fftsize, shiftsize, win, windual, iteration, phase_temp, Ls)
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
            %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
            x = amp_corr .* exp(1i * phase_temp);
            
            for i = 1:iteration
                
                % イテレーション回数の印字
                %fprintf('    Iteration : %d \n', i);
                
                % 位相を更新
                S1 = ISTFT(x, windual, shiftsize, fftsize, Ls);
                S2 = STFT(S1, win, shiftsize, fftsize, Ls);
                
                % スペクトログラムの更新，ampは所望に変更，位相だけ保管
                % アダマール積に注意
                x = amp_corr .* exp( 1i * angle(S2) );
                
            end
            
            spectrum = x;
            
        end
        
        function spectrum = GLA_ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, Ls)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 16 Apr. 2019.
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
            % GLA + ADMM
            %%%%%%%%%%%%%%%%%%%%
            
            % 初期値
            %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
            %       z =  x: ADMMで解く最適化問題に落とし込むため，新たに変数zを用意
            %       u = 0 : ラグランジュ未定乗数
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            
            for i = 1:iteration
                % イテレーション回数の印字
                %fprintf('    Iteration : %d \n', i);
                
                % スペクトログラムの更新，ampは所望に変更，位相はだけ保管
                % アダマール積に注意
                % x = amp_corr .* ( (z - u) / abs(z - u) ); 間違ってる気がする
                x = amp_corr .* exp( 1i * angle(z - u) );
                
                % 新たな変数であるzを更新
                %   下準備として STFT( ISTFT() )をおこなう
                S1 = ISTFT(x+u, windual, shiftsize, fftsize, Ls);
                S2 = STFT(S1, win, shiftsize, fftsize, Ls);
                %   zの更新
                z = ( rho*S2 + x + u ) / (1 + rho);
                
                % ラグランジュ未定乗数uの更新
                u = u + x - z;
                
            end
            
            spectrum = x;
            
        end
        
        function spectrum = Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, Ls)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 11 May. 2019.
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
            % GLA + ADMM + prop
            %%%%%%%%%%%%%%%%%%%%
            
            % 初期値
            %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
            %       z =  x: ADMMで解く最適化問題に落とし込むため，新たに変数zを用意
            %       u = 0 : ラグランジュ未定乗数
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            
            for i = 1:iteration
                
                % イテレーション回数の印字
                %fprintf('    Iteration : %d \n', i);
                
                % xの更新
                %   STFT( ISTFT() )をおこなう
                S1 = ISTFT(z-u, windual, shiftsize, fftsize, Ls);
                x = STFT(S1, win, shiftsize, fftsize, Ls);
                
                % 新たな変数であるzを更新
                %   アダマール積に注意
                z = ( ( rho*amp_corr + abs( x + u ) ) / (1 + rho) ) .* exp( 1i * angle(x + u) );
                
                % ラグランジュ未定乗数uの更新
                u = u + x - z;
                
            end
            
            spectrum = x;
            
        end
        
        function [spectrum, min_alpha] = General(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, Ls)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 13 May. 2019.
            %
            % [inputs]
            %   amp_corr: correct amplitude ( FrequencyBin * Frames)
            %   rho: ADMM Parameter ( ρ = 0.1, 0.2, 10, 100)
            %   fftsize: FFT length (=8192)
            %   shiftsize: frame shift (default: fftSize/2)
            %   window: window function used in STFT (fftSize x 1) or choose used
            %   iteration: iteration number
            %   phase_temp: Random phase used as initial value
            %   spectrum_corr: correct spectrum
            %   freq: fft length/2 (8192/2)
            %
            % [outputs]
            %   spectrum: frequency-domain (fftSize/2+1 x frames)
            %
            
            
            % 100.1 + 100.1i の複素数を用意し，スペクトル距離最小となるスペクトルを得る
            temp_comp = complex(double(100.1), double(100.1));
            min_x = repmat(temp_comp, freq, frames);
            
            
            % alphaの更新
            %       0.1ずつインクリメントしながらalphaの値を更新
            %for alpha = 0.05:0.05:0.95
            for alpha = 0.05:0.05:0.95
                
                % alphaの更新回数の印字
                %fprintf('    alpha : %d \n', alpha);
            
                % 初期値
                %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
                %       z =  x: ADMMで解く最適化問題に落とし込むため，新たに変数zを用意
                %       u = 0 : ラグランジュ未定乗数
                %       min_x = 0 : alpha更新に伴ってスペクトル距離最小となるスペクトルを得る
                x = amp_corr .* exp(1i * phase_temp);
                z = x;
                u = zeros(freq, frames);
                
                
                %%%%%%%%%%%%%%%%%%%%
                % 一般化ADMM
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % xを更新(アダマール積に注意)
                    x = ( ( rho*amp_corr + ( 1 - alpha )*abs(z - u) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(z - u) );
                    
                    % 新たな変数であるzを更新
                    %   下準備として STFT( ISTFT() )をおこなう
                    S1 = ISTFT(x+u, windual, shiftsize, fftsize, Ls);
                    S2 = STFT(S1, win, shiftsize, fftsize, Ls);
                    %   zの更新
                    z = ( alpha*(x + u) + rho*S2 ) / ( alpha + rho );
                    
                    % ラグランジュ未定乗数uの更新
                    u = u + x - z;
                    
                end
                
                % 正解スペクトルとの差をみる
                temp_err_x = norm(spectrum_corr - x, 'fro');
                
                % 距離の印字
                %fprintf('    distance : %d \n', temp_err_x);
                
                % 差が小さければ更新
                if norm(spectrum_corr - min_x, 'fro') > temp_err_x
                    min_x = x;
                    min_alpha = alpha;
                    min_err_x = temp_err_x;
                end
                    
            end
            
            spectrum = min_x;
            
            % 印字
            fprintf('  min distance : %d ,   min alpha : %d \n', min_err_x, min_alpha);
            
        end
        
    end
end