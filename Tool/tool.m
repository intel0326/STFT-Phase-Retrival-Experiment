classdef tool
    methods (Static)
        function [windual, Spe, Ls, signal_len] = AudioReadMethod(filename, total_sec, freq, fftsize, shiftsize, win)
            %サンプリングレートのみ取得
            %いらない音源の部分は「~」で消しておく
            [~,fs] = audioread(filename);
            %読み込む長さを計算
            samples = [1, total_sec*fs];
            %音源取得
            [music, fs] = audioread(filename, samples);
            %ダウンサンプリング
            music = resample(music, freq, fs);
            %ステレオをモノラル化
            music=mean(music, 2);
            %musicの長さを取得，後にSTFTで使用
            signal_len=length(music);
            %逆の窓を合成
            windual = winDual(win, shiftsize);
            % !! Ls must be even number due to our STFT/iSTFT implementation !!
            Ls = ceil((signal_len+2*(fftsize-shiftsize)-fftsize)/shiftsize)*shiftsize+fftsize;
            % STFT(0パディングもSTFT.m内で行う)
            Spe = STFT(music, win, shiftsize, fftsize, Ls, signal_len);
        end
        
        function spectrum = GLA(amp_corr, fftsize, shiftsize, win, windual, iteration, phase_temp, Ls, signal_len)
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
                S1 = ISTFT(x, windual, shiftsize, fftsize, Ls, signal_len);
                S2 = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                
                % スペクトログラムの更新，ampは所望に変更，位相だけ保管
                % アダマール積に注意
                x = amp_corr .* exp( 1i * angle(S2) );
                
            end
            
            spectrum = x;
            
        end
        
        function spectrum = GLA_ADMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, Ls, signal_len)
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
                S1 = ISTFT(x+u, windual, shiftsize, fftsize, Ls, signal_len);
                S2 = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                %   zの更新
                z = ( rho*S2 + x + u ) / (1 + rho);
                
                % ラグランジュ未定乗数uの更新
                u = u + x - z;
                
            end
            
            spectrum = x;
            
        end
        
        function spectrum = Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, Ls, signal_len)
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
                
                % xの更新
                %   STFT( ISTFT() )をおこなう
                S1 = ISTFT(z-u, windual, shiftsize, fftsize, Ls, signal_len);
                x = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                
                % 新たな変数であるzを更新
                %   アダマール積に注意
                z = ( ( rho*amp_corr + abs( x + u ) ) / (1 + rho) ) .* exp( 1i * angle(x + u) );
                
                % ラグランジュ未定乗数uの更新
                u = u + x - z;
                
            end
            
            spectrum = x;
            
        end

        function spectrum = Prop_batch(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, Ls, signal_len)
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
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            
            % 今回の音源では，freq*frames = 322164点．
            % したがってバッチサイズを約数のb=26847にする．
            % さらに，バッチサイズが26847であるため，イテレーションは 322164/26847=12 となる
            %freq*frames
            b = 26847;
            batch_iteration = ( freq*frames ) / b;
            
            for i = 1:iteration
            
                batch = randperm(freq*frames);
                
                for k = 1 : batch_iteration
                    
                    % xの更新
                    %   STFT( ISTFT() )をおこなう
                    S1 = ISTFT(z-u, windual, shiftsize, fftsize, Ls, signal_len);
                    x = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                    
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
                
            end
            
            spectrum = x;
            
        end
        
        function spectrum = PropWeight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, Ls, signal_len, Delta)
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
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            w = ( amp_corr + Delta ).^-1;
            
            for i = 1:iteration
                
                % xの更新
                %   STFT( ISTFT() )をおこなう
                S1 = ISTFT(z-u, windual, shiftsize, fftsize, Ls, signal_len);
                x = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                
                % 新たな変数であるzを更新
                %   アダマール積「.*」と要素毎の除算「./」に注意
                temp1_A = w.*(rho*amp_corr) + abs( x + u );
                temp2_A = 1 + rho*w;
                A2 = temp1_A ./ temp2_A;
                z = A2 .* exp( 1i * angle(x + u) );
                
                % ラグランジュ未定乗数uの更新
                u = u + x - z;
                
            end
            
            spectrum = x;
            
        end

        function spectrum = Prop_batch_weight(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, Ls, signal_len, Delta)
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
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            w = ( amp_corr + Delta ).^-1;
            
            % 今回の音源では，freq*frames = 322164点．
            % したがってバッチサイズを約数のb=26847にする．
            % さらに，バッチサイズが26847であるため，イテレーションは 322164/26847=12 となる
            %freq*frames
            b = 26847;
            batch_iteration = ( freq*frames ) / b;
            
            for i = 1:iteration
            
                batch = randperm(freq*frames);
                
                for k = 1 : batch_iteration
                    
                    % xの更新
                    %   STFT( ISTFT() )をおこなう
                    S1 = ISTFT(z-u, windual, shiftsize, fftsize, Ls, signal_len);
                    x = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                    
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
                
            end
            
            spectrum = x;
            
        end
        
        
        function [spectrum, min_alpha] = General(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, Ls, signal_len)
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
            min_err_x = norm(spectrum_corr - spectrum_corr, 'fro');
            min_alpha = 0.0;
            
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
                    S1 = ISTFT(x+u, windual, shiftsize, fftsize, Ls, signal_len);
                    S2 = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
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
        
        function [spectrum, min_alpha] = DouglasRachfordSplitting(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, Ls, signal_len, gamma)
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
            min_err_x = norm(spectrum_corr - spectrum_corr, 'fro');
            min_alpha = 0.0;
            
            % alphaの更新
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]

                % 初期値
                %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
                %       s =  x: Douglas-Rachford Splitting Algorithmで解く最適化問題に落とし込むため，変数sを用意
                x = amp_corr .* exp(1i * phase_temp);
                s = x;
                
                %%%%%%%%%%%%%%%%%%%%
                % Douglas-Rachford Splitting Algorithm
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % xを更新(アダマール積に注意)
                    x = ( ( rho*amp_corr + ( 1 - alpha )*abs(s) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(s) );
                    
                    % sを更新
                    %   下準備として STFT( ISTFT() )をおこなう
                    S1 = ISTFT(2*x-s, windual, shiftsize, fftsize, Ls, signal_len);
                    S2 = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                    %   sの更新
                    s = s + gamma * ( ( alpha*(2*x-s) + rho*S2 ) / ( alpha + rho ) - x );
                                   
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
        
        function [spectrum, min_alpha] = SDMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, Ls, signal_len)
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
            min_err_x = norm(spectrum_corr - spectrum_corr, 'fro');
            min_alpha = 0.0;
            
            % alphaの更新
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]

                % 初期値
                %       x = amp_corr .* exp(1i * phase_temp) : 所望の振幅とランダムな位相によるスペクトル
                %       y1, y2 =  x: SDMMで解く最適化問題に落とし込むため，変数y1,y2を用意
                x = amp_corr .* exp(1i * phase_temp);
                y1 = x;
                y2 = x;
                z1 = zeros(freq, frames);
                z2 = zeros(freq, frames);
                
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
                    S1 = ISTFT(x+z2, windual, shiftsize, fftsize, Ls, signal_len);
                    S2 = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                    %   スペクトログラムの値域集合への写像により位相を修正
                    y2 = ( alpha*(x+z2) + rho*S2 ) / ( alpha + rho );
                    
                    % z1,z2の更新
                    z1 = z1 + x - y1;
                    z2 = z2 + x - y2;
                                   
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
        
        function [spectrum, min_alpha] = PPXA(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, Ls, signal_len, gamma)
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
            min_err_x = norm(spectrum_corr - spectrum_corr, 'fro');
            min_alpha = 0.0;
            
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
                x = amp_corr .* exp(1i * phase_temp);
                y1 = x;
                y2 = x;
                %p1 = zeros(freq, frames);
                %p2 = zeros(freq, frames);
                
                %%%%%%%%%%%%%%%%%%%%
                % APXX
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % p1を更新
                    %   理想的な振幅Aの集合へ写像することにより振幅を修正
                    p1 = ( ( rho*amp_corr + ( w1 - alpha*w1 )*abs(y1) ) / ( w1 - alpha*w1 + rho ) ) .* exp( 1i * angle(y1) );

                    % p2を更新
                    %   下準備として STFT( ISTFT() )をおこなう
                    S1 = ISTFT(y2, windual, shiftsize, fftsize, Ls, signal_len);
                    S2 = STFT(S1, win, shiftsize, fftsize, Ls, signal_len);
                    %   スペクトログラムの値域集合への写像により位相を修正
                    p2 = ( alpha*w2*(y2) + rho*S2 ) / ( alpha*w2 + rho );
                    
                    % pの更新
                    p = w1*p1 + w2*p2;
                    
                    % y1,y2の更新
                    y1 = y1 + gamma*(2*p - x - p1);
                    y2 = y2 + gamma*(2*p - x - p2);
                    
                    % xを更新
                    x = x + gamma*(p - x);
                   
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

        function OutputMethod(spectrum, windual, shiftsize, fftsize, Ls, freq, rho, outputDir, signal_len, name)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 27 May. 2019.
            %
            
            % 振幅を正規化する関数
            Normalize = @(x) x/max(abs(x));
            % スペクトルを時間信号に変換
            signal = ISTFT(spectrum, windual, shiftsize, fftsize, Ls, signal_len);
            % 音源をwave形式で出力
            audiowrite(sprintf('%s/signal_rho_%.2f/signal_%s_rho=%.2f.wav', outputDir, rho, name, rho), Normalize(signal), freq);
        
        end
        
    end
end