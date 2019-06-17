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




