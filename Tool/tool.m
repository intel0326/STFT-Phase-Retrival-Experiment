classdef tool
    methods (Static)
        function [windual, Spe, signal_len, music] = AudioReadMethod(filename, start_sec, end_sec, freq, fftsize, shiftsize, win, STFT_type)
            %サンプリングレートのみ取得
            %いらない音源の部分は「~」で消しておく
            [~,fs] = audioread(filename);
            %読み込む長さを計算
            samples = [1, end_sec*fs];
            %音源取得
            [music, fs] = audioread(filename, samples);
            %ダウンサンプリング
            music = resample(music, freq, fs);
            %ステレオをモノラル化
            music = music(:, 1); %mean(true_sound, 2);
            music = music(start_sec*freq+1:end);
            %musicの長さを取得，後にSTFTで使用
            signal_len=length(music);
            %逆の窓を合成
            windual = winDual(win, shiftsize);
            % STFT(0パディングもSTFT.m内で行う)
            Spe = STFT(music, win, shiftsize, fftsize, STFT_type);
        end
        
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
        
        function [spectrum, sound, amp_error, A_error] = ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, signal_len, STFT_type, A_weight)
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
            % ADMM (矢田部法)
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
        
        function [spectrum, sound, amp_error, A_error] = Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, signal_len, STFT_type, A_weight)
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
            %       amp_error = 0 : 振幅間のフロベニウスノルム
            %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            for i = 1:iteration
                
                % xの更新
                %   STFT( ISTFT() )をおこなう
                sound = ISTFT(z-u, windual, shiftsize, fftsize, signal_len, STFT_type);
                x = STFT(sound, win, shiftsize, fftsize, STFT_type);
                
                % 新たな変数であるzを更新
                %   アダマール積に注意
                z = ( ( rho*amp_corr + abs(x + u) ) / (1 + rho) ) .* exp( 1i * angle(x + u) );
                
                % ラグランジュ未定乗数uの更新
                u = u + x - z;
                
                % 誤差を出力
                amp_error(i) = norm(abs(x)-amp_corr,'fro');
                A_error(i) = norm(diag(A_weight)*(abs(x)-amp_corr),'fro');
                
            end
            
            spectrum = x;
            
        end

        function [spectrum, sound, amp_error, A_error] = Prop_batch(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, signal_len, STFT_type, A_weight)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 27 May. 2019.
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
            % GLA + ADMM + prop + バッチ処理
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
            
            % 今回の音源では，freq*frames = 322677点．
            % したがってバッチサイズを約数のb=999にする．
            % さらに，バッチサイズが999であるため，イテレーションは 322677/999=323 となる
            % 6/11 イテレーションを323になるように調整しないか
            b = 999;
            batch_iteration = ( freq*frames ) / b;
            
            for i = 1:iteration
            
                batch = randperm(freq*frames);
                
                for k = 1 : batch_iteration
                    
                    % xの更新
                    %   STFT( ISTFT() )をおこなう
                    sound = ISTFT(z-u, windual, shiftsize, fftsize, signal_len, STFT_type);
                    x = STFT(sound, win, shiftsize, fftsize, STFT_type);
                    
                    % zを更新
                    %   zに一度x+uを代入し，バッチ処理するインデックスのみを後ほど取り出す
                    z = x + u;
                    %   ある行列からバッチ処理する要素のインデックスをランダムに抽出したbatchベクトルから，
                    %   ループ毎にバッチサイズb分インデックスを取得
                    Extract = batch( b*(k-1)+1 : b*k );
                    %   アダマール積に注意しながらzを更新
                    z( Extract ) = ( ( rho*amp_corr( Extract ) + abs( z(Extract) ) ) / (1 + rho) ) .* exp( 1i * angle( z(Extract) ) );
                    
                    % ラグランジュ未定乗数uの更新
                    u = u + x - z;
                    
                end
                
                % 誤差を出力
                amp_error(i) = norm(abs(x)-amp_corr,'fro');
                A_error(i) = norm(diag(A_weight)*(abs(x)-amp_corr),'fro');
                
            end
            
            spectrum = x;
            
        end
        
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

        function [spectrum, sound, amp_error, A_error] = Prop_batch_weight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, signal_len, Delta, STFT_type, A_weight)
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
            % GLA + ADMM + prop + バッチ処理 + 振幅weight
            %%%%%%%%%%%%%%%%%%%%
            
            % 初期値
            %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
            %       z =  x: ADMMで解く最適化問題に落とし込むため，新たに変数zを用意
            %       u = 0 : ラグランジュ未定乗数
            %       w : 振幅の発散を防ぐ重み
            %       amp_error = 0 : 振幅間のフロベニウスノルム
            %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            w = ( amp_corr + Delta ).^-1;
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            % 今回の音源では，freq*frames = 322677点．
            % したがってバッチサイズを約数のb=999にする．
            % さらに，バッチサイズが999であるため，イテレーションは 322677/999=323 となる
            % 6/11 イテレーションを323になるように調整しないか
            %freq*frames
            b = 999;
            batch_iteration = ( freq*frames ) / b;
            
            for i = 1:iteration
            
                batch = randperm(freq*frames);
                
                for k = 1 : batch_iteration
                    
                    % xの更新
                    %   STFT( ISTFT() )をおこなう
                    sound = ISTFT(z-u, windual, shiftsize, fftsize, signal_len, STFT_type);
                    x = STFT(sound, win, shiftsize, fftsize, STFT_type);
                    
                    % zを更新
                    %   zに一度x+uを代入し，バッチ処理するインデックスのみを後ほど取り出す
                    z = x + u;
                    %   ある行列からバッチ処理する要素のインデックスをランダムに抽出したbatchベクトルから，
                    %   ループ毎にバッチサイズb分インデックスを取得
                    Extract = batch( b*(k-1)+1 : b*k );
                    %   アダマール積「.*」と要素毎の除算「./」に注意
                    temp1_A = w( Extract ).*( rho*amp_corr( Extract ) ) + abs( z(Extract) );
                    temp2_A = 1 + rho*w( Extract );
                    A2 = temp1_A ./ temp2_A;
                    z( Extract ) = A2 .* exp( 1i * angle( z(Extract) ) );
                    
                    % ラグランジュ未定乗数uの更新
                    u = u + x - z;
                    
                end
                
                % 誤差を出力
                amp_error(i) = norm(abs(x)-amp_corr,'fro');
                A_error(i) = norm(diag(A_weight)*(abs(x)-amp_corr),'fro');
                
            end
            
            spectrum = x;
            
        end
        
        
        function [spectrum, sound, min_amp_error, min_A_error, min_alpha] = General(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, signal_len, STFT_type, A_weight)
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
            
            
            % 1.e1000(無限) + 1.e1000i(無限) の複素数を用意し，スペクトル距離最小となるスペクトルを得る
            temp_comp = complex(double(1.e1000), double(1.e1000));
            min_x = repmat(temp_comp, freq, frames);
            min_alpha = 0.0;
            min_amp_error = zeros(iteration,1);
            min_A_error = zeros(iteration,1);
            min_sound = zeros(signal_len,1);
            
            % alphaの更新
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]
            
                % alphaの更新回数の印字
                %fprintf('    alpha : %d \n', alpha);
            
                % 初期値
                %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
                %       z =  x: ADMMで解く最適化問題に落とし込むため，新たに変数zを用意
                %       u = 0 : ラグランジュ未定乗数
                %       min_x = 0 : alpha更新に伴ってスペクトル距離最小となるスペクトルを得る
                %       amp_error = 0 : 振幅間のフロベニウスノルム
                %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
                x = amp_corr .* exp(1i * phase_temp);
                z = x;
                u = zeros(freq, frames);
                amp_error = zeros(iteration,1);
                A_error = zeros(iteration,1);
                
                %%%%%%%%%%%%%%%%%%%%
                % 一般化ADMM
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % xを更新(アダマール積に注意)
                    x = ( ( rho*amp_corr + ( 1 - alpha )*abs(z - u) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(z - u) );
                    
                    % 新たな変数であるzを更新
                    %   下準備として STFT( ISTFT() )をおこなう
                    S1 = ISTFT(x+u, windual, shiftsize, fftsize, signal_len, STFT_type);
                    S2 = STFT(S1, win, shiftsize, fftsize, STFT_type);
                    %   zの更新
                    z = ( alpha*(x + u) + rho*S2 ) / ( alpha + rho );
                    
                    % ラグランジュ未定乗数uの更新
                    u = u + x - z;
                    
                    % イテレーション毎に評価するために一度STFTの値域空間に射影
                    temp_sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                    temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                    
                    % 誤差を出力
                    amp_error(i) = norm(abs(temp_spectrum)-amp_corr,'fro');
                    A_error(i) = norm(diag(A_weight)*(abs(temp_spectrum)-amp_corr),'fro');
                    
                end
                
                % 正解スペクトルとの差をみる
                %temp_err_x = norm(spectrum_corr - x, 'fro');
                temp_err_x = A_error(end);
                
                % 距離の印字
                %fprintf('    distance : %d \n', temp_err_x);
                
                % 差が小さければ更新
                %if norm(diag(A_weight)*(abs(min_x)-amp_corr),'fro') > temp_err_x
                if min_A_error(end) > temp_err_x
                    min_x = temp_spectrum;
                    min_sound = temp_sound;
                    min_alpha = alpha;
                    % 誤差を更新
                    min_amp_error = amp_error;
                    min_A_error = A_error;
                end
                    
            end
            
            spectrum = min_x;
            sound = min_sound;
            
            % 印字
            fprintf('  min distance A filter : %d ,   min alpha : %d \n', min_A_error, min_alpha);
            
        end
        
        function [spectrum, sound, min_amp_error, min_A_error, min_alpha] = DouglasRachfordSplitting(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, signal_len, gamma, STFT_type, A_weight)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 24 May. 2019.
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
                
            % 1.e1000(無限) + 1.e1000i(無限) の複素数を用意し，スペクトル距離最小となるスペクトルを得る
            temp_comp = complex(double(1.e1000), double(1.e1000));
            min_x = repmat(temp_comp, freq, frames);
            min_alpha = 0.0;
            min_amp_error = zeros(iteration,1);
            min_A_error = zeros(iteration,1);
            min_sound = zeros(signal_len,1);
            
            % alphaの更新
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]

                % 初期値
                %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
                %       s =  x: Douglas-Rachford Splitting Algorithmで解く最適化問題に落とし込むため，変数sを用意
                %       amp_error = 0 : 振幅間のフロベニウスノルム
                %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
                x = amp_corr .* exp(1i * phase_temp);
                s = x;
                amp_error = zeros(iteration,1);
                A_error = zeros(iteration,1);
                
                %%%%%%%%%%%%%%%%%%%%
                % Douglas-Rachford Splitting Algorithm
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % xを更新(アダマール積に注意)
                    x = ( ( rho*amp_corr + ( 1 - alpha )*abs(s) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(s) );
                    
                    % sを更新
                    %   下準備として STFT( ISTFT() )をおこなう
                    S1 = ISTFT(2*x-s, windual, shiftsize, fftsize, signal_len, STFT_type);
                    S2 = STFT(S1, win, shiftsize, fftsize, STFT_type);
                    %   sの更新
                    s = s + gamma * ( ( alpha*(2*x-s) + rho*S2 ) / ( alpha + rho ) - x );
                    
                    % イテレーション毎に評価するために一度STFTの値域空間に射影
                    temp_sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                    temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                    
                    % 誤差を出力
                    amp_error(i) = norm(abs(temp_spectrum)-amp_corr,'fro');
                    A_error(i) = norm(diag(A_weight)*(abs(temp_spectrum)-amp_corr),'fro');
                                   
                end
                
                % 正解スペクトルとの差をみる
                %temp_err_x = norm(spectrum_corr - x, 'fro');
                temp_err_x = A_error(end);
                
                % 距離の印字
                %fprintf('    distance : %d \n', temp_err_x);
                
                % 差が小さければ更新
                if min_A_error(end) > temp_err_x
                    min_x = temp_spectrum;
                    min_sound = temp_sound;
                    min_alpha = alpha;
                    % 誤差を更新
                    min_amp_error = amp_error;
                    min_A_error = A_error;
                end
                    
            end
            
            spectrum = min_x;
            sound = min_sound;
            
            % 印字
            fprintf('  min distance A filter : %d ,   min alpha : %d \n', min_A_error, min_alpha);
            
        end
        
        function [spectrum, sound, min_amp_error, min_A_error, min_alpha] = SDMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, signal_len, STFT_type, A_weight)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 24 May. 2019.
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
                
            % 1.e1000(無限) + 1.e1000i(無限) の複素数を用意し，スペクトル距離最小となるスペクトルを得る
            temp_comp = complex(double(1.e1000), double(1.e1000));
            min_x = repmat(temp_comp, freq, frames);
            min_alpha = 0.0;
            min_amp_error = zeros(iteration,1);
            min_A_error = zeros(iteration,1);
            min_sound = zeros(signal_len,1);
            
            % alphaの更新
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]

                % 初期値
                %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
                %       y1, y2 =  x: SDMMで解く最適化問題に落とし込むため，変数y1,y2を用意
                %       amp_error = 0 : 振幅間のフロベニウスノルム
                %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
                x = amp_corr .* exp(1i * phase_temp);
                y1 = x;
                y2 = x;
                z1 = zeros(freq, frames);
                z2 = zeros(freq, frames);
                amp_error = zeros(iteration,1);
                A_error = zeros(iteration,1);
                
                %%%%%%%%%%%%%%%%%%%%
                % SDMM
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % xを更新
                    x = ( ( y1 - z1 ) + ( y2 - z2 ) ) / 2;
                    
                    %y1を更新
                    %   理想的な振幅Aの集合へ写像することにより振幅を修正
                    y1 = ( ( rho*amp_corr + ( 1 - alpha )*abs(x+z1) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(x+z1) );
                    
                    %y2を更新
                    %   下準備として STFT( ISTFT() )をおこなう
                    S1 = ISTFT(x+z2, windual, shiftsize, fftsize, signal_len, STFT_type);
                    S2 = STFT(S1, win, shiftsize, fftsize, STFT_type);
                    %   スペクトログラムの値域集合への写像により位相を修正
                    y2 = ( alpha*(x+z2) + rho*S2 ) / ( alpha + rho );
                    
                    % z1,z2の更新
                    z1 = z1 + x - y1;
                    z2 = z2 + x - y2;
                    
                    % イテレーション毎に評価するために一度STFTの値域空間に射影
                    temp_sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                    temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                    
                    % 誤差を出力
                    amp_error(i) = norm(abs(temp_spectrum)-amp_corr,'fro');
                    A_error(i) = norm(diag(A_weight)*(abs(temp_spectrum)-amp_corr),'fro');
                                   
                end
                
                % 正解スペクトルとの差をみる
                %temp_err_x = norm(spectrum_corr - x, 'fro');
                temp_err_x = A_error(end);
                
                % 距離の印字
                %fprintf('    distance : %d \n', temp_err_x);
                
                % 差が小さければ更新
                %if norm(diag(A_weight)*(abs(min_x)-amp_corr),'fro') > temp_err_x
                if min_A_error(end) > temp_err_x
                    min_x = temp_spectrum;
                    min_sound = temp_sound;
                    min_alpha = alpha;
                    % 誤差を更新
                    min_amp_error = amp_error;
                    min_A_error = A_error;
                end
                    
            end
            
            spectrum = min_x;
            sound = min_sound;
            
            % 印字
            fprintf('  min distance A filter : %d ,   min alpha : %d \n', min_A_error, min_alpha);
            
        end
        
        function [spectrum, sound, min_amp_error, min_A_error, min_alpha] = PPXA(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, signal_len, gamma, STFT_type, A_weight)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 25 May. 2019.
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
                
            % 1.e1000(無限) + 1.e1000i(無限) の複素数を用意し，スペクトル距離最小となるスペクトルを得る
            temp_comp = complex(double(1.e1000), double(1.e1000));
            min_x = repmat(temp_comp, freq, frames);
            min_alpha = 0.0;
            min_amp_error = zeros(iteration,1);
            min_A_error = zeros(iteration,1);
            min_sound = zeros(signal_len,1);
            
            % 重みωを決定
            w1 = 1/2;
            w2 = 1 - w1;
            
            % alphaの更新
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]

                % 初期値
                %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
                %       y1, y2 =  x: SDMMで解く最適化問題に落とし込むため，変数y1,y2を用意
                %       amp_error = 0 : 振幅間のフロベニウスノルム
                %       A_error = 0 : A特性付きの振幅間のフロベニウスノルム
                x = amp_corr .* exp(1i * phase_temp);
                y1 = x;
                y2 = x;
                amp_error = zeros(iteration,1);
                A_error = zeros(iteration,1);
                
                %%%%%%%%%%%%%%%%%%%%
                % APXX
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % p1を更新
                    %   理想的な振幅Aの集合へ写像することにより振幅を修正
                    p1 = ( ( rho*amp_corr + ( w1 - alpha*w1 )*abs(y1) ) / ( w1 - alpha*w1 + rho ) ) .* exp( 1i * angle(y1) );

                    % p2を更新
                    %   下準備として STFT( ISTFT() )をおこなう
                    S1 = ISTFT(y2, windual, shiftsize, fftsize, signal_len, STFT_type);
                    S2 = STFT(S1, win, shiftsize, fftsize, STFT_type);
                    %   スペクトログラムの値域集合への写像により位相を修正
                    p2 = ( alpha*w2*(y2) + rho*S2 ) / ( alpha*w2 + rho );
                    
                    % pの更新
                    p = w1*p1 + w2*p2;
                    
                    % y1,y2の更新
                    y1 = y1 + gamma*(2*p - x - p1);
                    y2 = y2 + gamma*(2*p - x - p2);
                    
                    % xを更新
                    x = x + gamma*(p - x);
                    
                    % イテレーション毎に評価するために一度STFTの値域空間に射影
                    temp_sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                    temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                    
                    % 誤差を出力
                    amp_error(i) = norm(abs(temp_spectrum)-amp_corr,'fro');
                    A_error(i) = norm(diag(A_weight)*(abs(temp_spectrum)-amp_corr),'fro');
                   
                end
                
                % 正解スペクトルとの差をみる
                %temp_err_x = norm(spectrum_corr - x, 'fro');
                temp_err_x = A_error(end);
                
                % 距離の印字
                %fprintf('    distance : %d \n', temp_err_x);
                
                % 差が小さければ更新
                %if norm(diag(A_weight)*(abs(min_x)-amp_corr),'fro') > temp_err_x
                if min_A_error(end) > temp_err_x
                    min_x = temp_spectrum;
                    min_sound = temp_sound;
                    min_alpha = alpha;
                    % 誤差を更新
                    min_amp_error = amp_error;
                    min_A_error = A_error;
                end
                    
            end
            
            spectrum = min_x;
            sound = min_sound;
            
            % 印字
            fprintf('  min distance A filter : %d ,   min alpha : %d \n', min_A_error, min_alpha);
            
        end
        
        function [err, fro] = evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 28 May. 2019.
            %
            
            % 推定した位相を算出
            phase_est= angle(spectrum_est);
            % 振幅1かつ推定した位相のスペクトルを算出
            spectrum_amp1_est = ones( size(amp_corr) ) .* exp( 1i * phase_est );
            % 位相差を算出し，評価
            err = immse(spectrum_amp1_corr, spectrum_amp1_est);
            % スペクトル差を算出し，評価
            fro = norm(spectrum - spectrum_est,'fro');
            
        end

        function OutputMethod(spectrum, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, name)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 27 May. 2019.
            %
            
            % 振幅を正規化する関数
            Normalize = @(x) x/max(abs(x));
            % スペクトルを時間信号に変換
            signal = ISTFT(spectrum, windual, shiftsize, fftsize, signal_len, STFT_type);
            % 音源をwave形式で出力
            audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, name, rho), Normalize(signal), freq);
        
        end
        
        function OutputPrintingMethod(name, amp_error, A_error, alpha)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 17 June. 2019.
            %
            
            if nargin < 4
                fprintf('    %s\n', name);
                fprintf('        amp Error : %d, Afilter Error : %d\n', amp_error, A_error);
            else
                fprintf('    %s', name);
                fprintf('        amp Error : %d, Afilter Error : %d, alpha : %d\n', amp_error, A_error, alpha);
            end
            

        end
        
    end
end




