function [spectrum, sound, amp_error, A_error] = PropWeight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, signal_len, Delta, STFT_type, A_weight)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 4 June. 2019.
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
            % GLA + ADMM + prop + 振幅weight
            %%%%%%%%%%%%%%%%%%%%
            
            % 初期値
            %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
            %       z =  x: ADMMで解く最適化問題に落とし込むため，新たに変数zを用意
            %       u = 0 : ラグランジュ未定乗数
            %       w : 振幅の発散を防ぐ重み(「.^」は要素単位のべき乗)
            %       amp_error = 0 : 振幅間のフロベニウスノルム
            %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            w = ( amp_corr + Delta ).^-1;
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            for i = 1:iteration
                
                % xの更新
                %   STFT( ISTFT() )をおこなう
                sound = ISTFT(z-u, windual, shiftsize, fftsize, signal_len, STFT_type);
                x = STFT(sound, win, shiftsize, fftsize, STFT_type);
                
                % 新たな変数であるzを更新
                %   アダマール積「.*」と要素毎の除算「./」に注意
                temp1_A = w.*(rho*amp_corr) + abs( x + u );
                temp2_A = 1 + rho*w;
                A2 = temp1_A ./ temp2_A;
                z = A2 .* exp( 1i * angle(x + u) );
                
                % ラグランジュ未定乗数uの更新
                u = u + x - z;
                
                % 誤差を出力
                amp_error(i) = norm(abs(x)-amp_corr,'fro');
                A_error(i) = norm(diag(A_weight)*(abs(x)-amp_corr),'fro');
                
            end
            
            spectrum = x;
            
end