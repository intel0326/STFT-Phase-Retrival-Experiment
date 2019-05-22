classdef tool
    methods (Static)
        function [music, Spe] = AudioReadMethod(filename, total_sec, freq, fftsize, shiftsize, window)
            %�T���v�����O���[�g�̂ݎ擾
            %����Ȃ������̕����́u~�v�ŏ����Ă���
            [~,fs] = audioread(filename);
            %�ǂݍ��ޒ������v�Z
            samples = [1, total_sec*fs];
            %�����擾
            [music, fs] = audioread(filename, samples);
            %�_�E���T���v�����O
            music = resample(music, freq, fs);
            %FFT�Ɠǂݍ��ޒ����𒲐�
            %floor( (10�b*16000 - 8192 + 2048) / 2048 )
            temp_len = floor( (total_sec*freq - fftsize + shiftsize) / shiftsize );
            length = temp_len*shiftsize + fftsize;
            %������FFT�p�ɃJ�b�g
            music = music(1:length);
            %�����M���ɃV���[�g�^�C���t�[���G�ϊ���������B�i�^�̕��f�X�y�N�g���O����)
            [Spe,~] = STFT(music', fftsize, shiftsize, window);
        end
        
        function spectrum = GLA(amp_corr, fftsize, shiftsize, window, iteration, phase_temp)
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
            
            % �����l
            %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
            x = amp_corr .* exp(1i * phase_temp);
            
            for i = 1:iteration
                
                % �C�e���[�V�����񐔂̈�
                %fprintf('    Iteration : %d \n', i);
                
                % �ʑ����X�V
                S1 = ISTFT(x, shiftsize, window);
                [S2,~] = STFT(S1, fftsize, shiftsize, window );
                
                % �X�y�N�g���O�����̍X�V�Camp�͏��]�ɕύX�C�ʑ������ۊ�
                % �A�_�}�[���ςɒ���
                x = amp_corr .* exp( 1i * angle(S2) );
                
            end
            
            spectrum = x;
            
        end
        
        function spectrum = GLA_ADMM(amp_corr, rho, fftsize, shiftsize, window, iteration, phase_temp, freq,  frames)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 16 Apr. 2019.
            %
            % [inputs]
            %   amp_corr: correct amplitude ( FrequencyBin * Frames)
            %   rho: ADMM Parameter ( �� = 0.1, 0.2, 10, 100)
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
            
            % �����l
            %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
            %       z =  x: ADMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�V���ɕϐ�z��p��
            %       u = 0 : ���O�����W������搔
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            
            for i = 1:iteration
                % �C�e���[�V�����񐔂̈�
                %fprintf('    Iteration : %d \n', i);
                
                % �X�y�N�g���O�����̍X�V�Camp�͏��]�ɕύX�C�ʑ��͂����ۊ�
                % �A�_�}�[���ςɒ���
                % x = amp_corr .* ( (z - u) / abs(z - u) ); �Ԉ���Ă�C������
                x = amp_corr .* exp( 1i * angle(z - u) );
                
                % �V���ȕϐ��ł���z���X�V
                %   �������Ƃ��� STFT( ISTFT() )�������Ȃ�
                S1 = ISTFT(x + u, shiftsize, window);
                [S2,~] = STFT(S1, fftsize, shiftsize, window );
                %   z�̍X�V
                z = ( rho*S2 + x + u ) / (1 + rho);
                
                % ���O�����W������搔u�̍X�V
                u = u + x - z;
                
            end
            
            spectrum = x;
            
        end
        
        function spectrum = Prop(amp_corr, rho, fftsize, shiftsize, window, iteration, phase_temp, freq,  frames)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 11 May. 2019.
            %
            % [inputs]
            %   amp_corr: correct amplitude ( FrequencyBin * Frames)
            %   rho: ADMM Parameter ( �� = 0.1, 0.2, 10, 100)
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
            
            % �����l
            %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
            %       z =  x: ADMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�V���ɕϐ�z��p��
            %       u = 0 : ���O�����W������搔
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            
            for i = 1:iteration
                
                % �C�e���[�V�����񐔂̈�
                %fprintf('    Iteration : %d \n', i);
                
                % x�̍X�V
                %   STFT( ISTFT() )�������Ȃ�
                S1 = ISTFT(z - u, shiftsize, window);
                [x,~] = STFT(S1, fftsize, shiftsize, window );
                
                % �V���ȕϐ��ł���z���X�V
                %   �A�_�}�[���ςɒ���
                z = ( ( rho*amp_corr + abs( x + u ) ) / (1 + rho) ) .* exp( 1i * angle(x + u) );
                
                % ���O�����W������搔u�̍X�V
                u = u + x - z;
                
            end
            
            spectrum = x;
            
        end
        
        function [spectrum, min_alpha] = General(amp_corr, rho, fftsize, shiftsize, window, iteration, phase_temp, freq, spectrum_corr,  frames)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 13 May. 2019.
            %
            % [inputs]
            %   amp_corr: correct amplitude ( FrequencyBin * Frames)
            %   rho: ADMM Parameter ( �� = 0.1, 0.2, 10, 100)
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
            
            
            % 100.1 + 100.1i �̕��f����p�ӂ��C�X�y�N�g�������ŏ��ƂȂ�X�y�N�g���𓾂�
            temp_comp = complex(double(100.1), double(100.1));
            min_x = repmat(temp_comp, freq, frames);
            
            
            % alpha�̍X�V
            %       0.1���C���N�������g���Ȃ���alpha�̒l���X�V
            %for alpha = 0.05:0.05:0.95
            for alpha = 0.05:0.05:0.95
                
                % alpha�̍X�V�񐔂̈�
                %fprintf('    alpha : %d \n', alpha);
            
                % �����l
                %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
                %       z =  x: ADMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�V���ɕϐ�z��p��
                %       u = 0 : ���O�����W������搔
                %       min_x = 0 : alpha�X�V�ɔ����ăX�y�N�g�������ŏ��ƂȂ�X�y�N�g���𓾂�
                x = amp_corr .* exp(1i * phase_temp);
                z = x;
                u = zeros(freq, frames);
                
                
                %%%%%%%%%%%%%%%%%%%%
                % ��ʉ�ADMM
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % x���X�V(�A�_�}�[���ςɒ���)
                    x = ( ( rho*amp_corr + ( 1 - alpha )*abs(z - u) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(z - u) );
                    
                    % �V���ȕϐ��ł���z���X�V
                    %   �������Ƃ��� STFT( ISTFT() )�������Ȃ�
                    S1 = ISTFT(x + u, shiftsize, window);
                    [S2,~] = STFT(S1, fftsize, shiftsize, window );
                    %   z�̍X�V
                    z = ( alpha*(x + u) + rho*S2 ) / ( alpha + rho );
                    
                    % ���O�����W������搔u�̍X�V
                    u = u + x - z;
                    
                end
                
                % �����X�y�N�g���Ƃ̍����݂�
                temp_err_x = norm(spectrum_corr - x, 'fro');
                
                % �����̈�
                %fprintf('    distance : %d \n', temp_err_x);
                
                % ������������΍X�V
                if norm(spectrum_corr - min_x, 'fro') > temp_err_x
                    min_x = x;
                    min_alpha = alpha;
                    min_err_x = temp_err_x;
                end
                    
            end
            
            spectrum = min_x;
            
            % ��
            fprintf('  min distance : %d ,   min alpha : %d \n', min_err_x, min_alpha);
            
        end
        
    end
end