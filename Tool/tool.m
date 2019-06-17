classdef tool
    methods (Static)
        function [windual, Spe, signal_len, music] = AudioReadMethod(filename, start_sec, end_sec, freq, fftsize, shiftsize, win, STFT_type)
            %�T���v�����O���[�g�̂ݎ擾
            %����Ȃ������̕����́u~�v�ŏ����Ă���
            [~,fs] = audioread(filename);
            %�ǂݍ��ޒ������v�Z
            samples = [1, end_sec*fs];
            %�����擾
            [music, fs] = audioread(filename, samples);
            %�_�E���T���v�����O
            music = resample(music, freq, fs);
            %�X�e���I�����m������
            music = music(:, 1); %mean(true_sound, 2);
            music = music(start_sec*freq+1:end);
            %music�̒������擾�C���STFT�Ŏg�p
            signal_len=length(music);
            %�t�̑�������
            windual = winDual(win, shiftsize);
            % STFT(0�p�f�B���O��STFT.m���ōs��)
            Spe = STFT(music, win, shiftsize, fftsize, STFT_type);
        end
        
        function [err, fro] = evaluation(spectrum, spectrum_amp1_corr, amp_corr, spectrum_est)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 28 May. 2019.
            %
            
            % ���肵���ʑ����Z�o
            phase_est= angle(spectrum_est);
            % �U��1�����肵���ʑ��̃X�y�N�g�����Z�o
            spectrum_amp1_est = ones( size(amp_corr) ) .* exp( 1i * phase_est );
            % �ʑ������Z�o���C�]��
            err = immse(spectrum_amp1_corr, spectrum_amp1_est);
            % �X�y�N�g�������Z�o���C�]��
            fro = norm(spectrum - spectrum_est,'fro');
            
        end

        function OutputMethod(spectrum, windual, shiftsize, fftsize, STFT_type, freq, rho, outputDir, signal_len, name)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 27 May. 2019.
            %
            
            % �U���𐳋K������֐�
            Normalize = @(x) x/max(abs(x));
            % �X�y�N�g�������ԐM���ɕϊ�
            signal = ISTFT(spectrum, windual, shiftsize, fftsize, signal_len, STFT_type);
            % ������wave�`���ŏo��
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




