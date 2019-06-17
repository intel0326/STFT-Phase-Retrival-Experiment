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
            
            % �����l
            %       spectrum = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
            %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
            %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
            spectrum = amp_corr .* exp(1i * phase_temp);
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            
            % �C�e���[�V�����̍Ō��l���Ԃւ̎ˉe�ŏI��炵�����̂ŁC
            % �C�e���[�V����1��ڂ̓��[�v�̊O����
            sound = ISTFT(spectrum, windual, shiftsize, fftsize, signal_len, STFT_type);
            spectrum = STFT(sound, win, shiftsize, fftsize, STFT_type);
            
            % �ŏ��̌덷���o��
            amp_error(1) = norm(abs(spectrum)-amp_corr,'fro');
            A_error(1) = norm(diag(A_weight)*(abs(spectrum)-amp_corr),'fro');
            
            
            for i = 2:iteration
                
                % �C�e���[�V�����񐔂̈�
                %fprintf('    Iteration : %d \n', i);
                
                % �U���W���ւ̎ˉe�Camp�͏��]�ɕύX�C�ʑ������ۊ�
                % �A�_�}�[���ςɒ���
                spectrum = amp_corr .* exp( 1i * angle(spectrum) );
                
                % �ʑ����X�V
                sound = ISTFT(spectrum, windual, shiftsize, fftsize, signal_len, STFT_type);
                spectrum = STFT(sound, win, shiftsize, fftsize, STFT_type);
                
                % �덷���o��
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
            % ADMM (��c���@)
            %%%%%%%%%%%%%%%%%%%%
            
            % �����l
            %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
            %       z =  x: ADMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�V���ɕϐ�z��p��
            %       u = 0 : ���O�����W������搔
            %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
            %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            for i = 1:iteration
                % �C�e���[�V�����񐔂̈�
                %fprintf('    Iteration : %d \n', i);
                
                % �X�y�N�g���O�����̍X�V�Camp�͏��]�ɕύX�C�ʑ��͂����ۊ�
                % �A�_�}�[���ςɒ���
                x = amp_corr .* exp( 1i * angle(z - u) );
                
                % �V���ȕϐ��ł���z���X�V
                %   �������Ƃ��� STFT( ISTFT() )�������Ȃ�
                temp_sound = ISTFT(x+u, windual, shiftsize, fftsize, signal_len, STFT_type);
                temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                %   z�̍X�V
                z = ( rho*temp_spectrum + x + u ) / (1 + rho);
                
                % ���O�����W������搔u�̍X�V
                u = u + x - z;
                
                % �C�e���[�V�������ɕ]�����邽�߂Ɉ�xSTFT�̒l���ԂɎˉe
                sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                spectrum = STFT(sound, win, shiftsize, fftsize, STFT_type);
                
                % �덷���o��
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
            %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
            %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            for i = 1:iteration
                
                % x�̍X�V
                %   STFT( ISTFT() )�������Ȃ�
                sound = ISTFT(z-u, windual, shiftsize, fftsize, signal_len, STFT_type);
                x = STFT(sound, win, shiftsize, fftsize, STFT_type);
                
                % �V���ȕϐ��ł���z���X�V
                %   �A�_�}�[���ςɒ���
                z = ( ( rho*amp_corr + abs(x + u) ) / (1 + rho) ) .* exp( 1i * angle(x + u) );
                
                % ���O�����W������搔u�̍X�V
                u = u + x - z;
                
                % �덷���o��
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
            % GLA + ADMM + prop + �o�b�`����
            %%%%%%%%%%%%%%%%%%%%
            
            % �����l
            %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
            %       z =  x: ADMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�V���ɕϐ�z��p��
            %       u = 0 : ���O�����W������搔
            %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
            %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            % ����̉����ł́Cfreq*frames = 322677�_�D
            % ���������ăo�b�`�T�C�Y��񐔂�b=999�ɂ���D
            % ����ɁC�o�b�`�T�C�Y��999�ł��邽�߁C�C�e���[�V������ 322677/999=323 �ƂȂ�
            % 6/11 �C�e���[�V������323�ɂȂ�悤�ɒ������Ȃ���
            b = 999;
            batch_iteration = ( freq*frames ) / b;
            
            for i = 1:iteration
            
                batch = randperm(freq*frames);
                
                for k = 1 : batch_iteration
                    
                    % x�̍X�V
                    %   STFT( ISTFT() )�������Ȃ�
                    sound = ISTFT(z-u, windual, shiftsize, fftsize, signal_len, STFT_type);
                    x = STFT(sound, win, shiftsize, fftsize, STFT_type);
                    
                    % z���X�V
                    %   z�Ɉ�xx+u�������C�o�b�`��������C���f�b�N�X�݂̂���قǎ��o��
                    z = x + u;
                    %   ����s�񂩂�o�b�`��������v�f�̃C���f�b�N�X�������_���ɒ��o����batch�x�N�g������C
                    %   ���[�v���Ƀo�b�`�T�C�Yb���C���f�b�N�X���擾
                    Extract = batch( b*(k-1)+1 : b*k );
                    %   �A�_�}�[���ςɒ��ӂ��Ȃ���z���X�V
                    z( Extract ) = ( ( rho*amp_corr( Extract ) + abs( z(Extract) ) ) / (1 + rho) ) .* exp( 1i * angle( z(Extract) ) );
                    
                    % ���O�����W������搔u�̍X�V
                    u = u + x - z;
                    
                end
                
                % �덷���o��
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
            % GLA + ADMM + prop + �U��weight
            %%%%%%%%%%%%%%%%%%%%
            
            % �����l
            %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
            %       z =  x: ADMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�V���ɕϐ�z��p��
            %       u = 0 : ���O�����W������搔
            %       w : �U���̔��U��h���d��(�u.^�v�͗v�f�P�ʂׂ̂���)
            %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
            %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            w = ( amp_corr + Delta ).^-1;
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            for i = 1:iteration
                
                % x�̍X�V
                %   STFT( ISTFT() )�������Ȃ�
                sound = ISTFT(z-u, windual, shiftsize, fftsize, signal_len, STFT_type);
                x = STFT(sound, win, shiftsize, fftsize, STFT_type);
                
                % �V���ȕϐ��ł���z���X�V
                %   �A�_�}�[���ρu.*�v�Ɨv�f���̏��Z�u./�v�ɒ���
                temp1_A = w.*(rho*amp_corr) + abs( x + u );
                temp2_A = 1 + rho*w;
                A2 = temp1_A ./ temp2_A;
                z = A2 .* exp( 1i * angle(x + u) );
                
                % ���O�����W������搔u�̍X�V
                u = u + x - z;
                
                % �덷���o��
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
            % GLA + ADMM + prop + �o�b�`���� + �U��weight
            %%%%%%%%%%%%%%%%%%%%
            
            % �����l
            %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
            %       z =  x: ADMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�V���ɕϐ�z��p��
            %       u = 0 : ���O�����W������搔
            %       w : �U���̔��U��h���d��
            %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
            %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
            x = amp_corr .* exp(1i * phase_temp);
            z = x;
            u = zeros(freq, frames);
            w = ( amp_corr + Delta ).^-1;
            amp_error = zeros(iteration,1);
            A_error = zeros(iteration,1);
            
            % ����̉����ł́Cfreq*frames = 322677�_�D
            % ���������ăo�b�`�T�C�Y��񐔂�b=999�ɂ���D
            % ����ɁC�o�b�`�T�C�Y��999�ł��邽�߁C�C�e���[�V������ 322677/999=323 �ƂȂ�
            % 6/11 �C�e���[�V������323�ɂȂ�悤�ɒ������Ȃ���
            %freq*frames
            b = 999;
            batch_iteration = ( freq*frames ) / b;
            
            for i = 1:iteration
            
                batch = randperm(freq*frames);
                
                for k = 1 : batch_iteration
                    
                    % x�̍X�V
                    %   STFT( ISTFT() )�������Ȃ�
                    sound = ISTFT(z-u, windual, shiftsize, fftsize, signal_len, STFT_type);
                    x = STFT(sound, win, shiftsize, fftsize, STFT_type);
                    
                    % z���X�V
                    %   z�Ɉ�xx+u�������C�o�b�`��������C���f�b�N�X�݂̂���قǎ��o��
                    z = x + u;
                    %   ����s�񂩂�o�b�`��������v�f�̃C���f�b�N�X�������_���ɒ��o����batch�x�N�g������C
                    %   ���[�v���Ƀo�b�`�T�C�Yb���C���f�b�N�X���擾
                    Extract = batch( b*(k-1)+1 : b*k );
                    %   �A�_�}�[���ρu.*�v�Ɨv�f���̏��Z�u./�v�ɒ���
                    temp1_A = w( Extract ).*( rho*amp_corr( Extract ) ) + abs( z(Extract) );
                    temp2_A = 1 + rho*w( Extract );
                    A2 = temp1_A ./ temp2_A;
                    z( Extract ) = A2 .* exp( 1i * angle( z(Extract) ) );
                    
                    % ���O�����W������搔u�̍X�V
                    u = u + x - z;
                    
                end
                
                % �덷���o��
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
            
            
            % 1.e1000(����) + 1.e1000i(����) �̕��f����p�ӂ��C�X�y�N�g�������ŏ��ƂȂ�X�y�N�g���𓾂�
            temp_comp = complex(double(1.e1000), double(1.e1000));
            min_x = repmat(temp_comp, freq, frames);
            min_alpha = 0.0;
            min_amp_error = zeros(iteration,1);
            min_A_error = zeros(iteration,1);
            min_sound = zeros(signal_len,1);
            
            % alpha�̍X�V
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]
            
                % alpha�̍X�V�񐔂̈�
                %fprintf('    alpha : %d \n', alpha);
            
                % �����l
                %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
                %       z =  x: ADMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�V���ɕϐ�z��p��
                %       u = 0 : ���O�����W������搔
                %       min_x = 0 : alpha�X�V�ɔ����ăX�y�N�g�������ŏ��ƂȂ�X�y�N�g���𓾂�
                %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
                %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
                x = amp_corr .* exp(1i * phase_temp);
                z = x;
                u = zeros(freq, frames);
                amp_error = zeros(iteration,1);
                A_error = zeros(iteration,1);
                
                %%%%%%%%%%%%%%%%%%%%
                % ��ʉ�ADMM
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % x���X�V(�A�_�}�[���ςɒ���)
                    x = ( ( rho*amp_corr + ( 1 - alpha )*abs(z - u) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(z - u) );
                    
                    % �V���ȕϐ��ł���z���X�V
                    %   �������Ƃ��� STFT( ISTFT() )�������Ȃ�
                    S1 = ISTFT(x+u, windual, shiftsize, fftsize, signal_len, STFT_type);
                    S2 = STFT(S1, win, shiftsize, fftsize, STFT_type);
                    %   z�̍X�V
                    z = ( alpha*(x + u) + rho*S2 ) / ( alpha + rho );
                    
                    % ���O�����W������搔u�̍X�V
                    u = u + x - z;
                    
                    % �C�e���[�V�������ɕ]�����邽�߂Ɉ�xSTFT�̒l���ԂɎˉe
                    temp_sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                    temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                    
                    % �덷���o��
                    amp_error(i) = norm(abs(temp_spectrum)-amp_corr,'fro');
                    A_error(i) = norm(diag(A_weight)*(abs(temp_spectrum)-amp_corr),'fro');
                    
                end
                
                % �����X�y�N�g���Ƃ̍����݂�
                %temp_err_x = norm(spectrum_corr - x, 'fro');
                temp_err_x = A_error(end);
                
                % �����̈�
                %fprintf('    distance : %d \n', temp_err_x);
                
                % ������������΍X�V
                %if norm(diag(A_weight)*(abs(min_x)-amp_corr),'fro') > temp_err_x
                if min_A_error(end) > temp_err_x
                    min_x = temp_spectrum;
                    min_sound = temp_sound;
                    min_alpha = alpha;
                    % �덷���X�V
                    min_amp_error = amp_error;
                    min_A_error = A_error;
                end
                    
            end
            
            spectrum = min_x;
            sound = min_sound;
            
            % ��
            fprintf('  min distance A filter : %d ,   min alpha : %d \n', min_A_error, min_alpha);
            
        end
        
        function [spectrum, sound, min_amp_error, min_A_error, min_alpha] = DouglasRachfordSplitting(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, signal_len, gamma, STFT_type, A_weight)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 24 May. 2019.
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
                
            % 1.e1000(����) + 1.e1000i(����) �̕��f����p�ӂ��C�X�y�N�g�������ŏ��ƂȂ�X�y�N�g���𓾂�
            temp_comp = complex(double(1.e1000), double(1.e1000));
            min_x = repmat(temp_comp, freq, frames);
            min_alpha = 0.0;
            min_amp_error = zeros(iteration,1);
            min_A_error = zeros(iteration,1);
            min_sound = zeros(signal_len,1);
            
            % alpha�̍X�V
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]

                % �����l
                %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
                %       s =  x: Douglas-Rachford Splitting Algorithm�ŉ����œK�����ɗ��Ƃ����ނ��߁C�ϐ�s��p��
                %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
                %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
                x = amp_corr .* exp(1i * phase_temp);
                s = x;
                amp_error = zeros(iteration,1);
                A_error = zeros(iteration,1);
                
                %%%%%%%%%%%%%%%%%%%%
                % Douglas-Rachford Splitting Algorithm
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % x���X�V(�A�_�}�[���ςɒ���)
                    x = ( ( rho*amp_corr + ( 1 - alpha )*abs(s) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(s) );
                    
                    % s���X�V
                    %   �������Ƃ��� STFT( ISTFT() )�������Ȃ�
                    S1 = ISTFT(2*x-s, windual, shiftsize, fftsize, signal_len, STFT_type);
                    S2 = STFT(S1, win, shiftsize, fftsize, STFT_type);
                    %   s�̍X�V
                    s = s + gamma * ( ( alpha*(2*x-s) + rho*S2 ) / ( alpha + rho ) - x );
                    
                    % �C�e���[�V�������ɕ]�����邽�߂Ɉ�xSTFT�̒l���ԂɎˉe
                    temp_sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                    temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                    
                    % �덷���o��
                    amp_error(i) = norm(abs(temp_spectrum)-amp_corr,'fro');
                    A_error(i) = norm(diag(A_weight)*(abs(temp_spectrum)-amp_corr),'fro');
                                   
                end
                
                % �����X�y�N�g���Ƃ̍����݂�
                %temp_err_x = norm(spectrum_corr - x, 'fro');
                temp_err_x = A_error(end);
                
                % �����̈�
                %fprintf('    distance : %d \n', temp_err_x);
                
                % ������������΍X�V
                if min_A_error(end) > temp_err_x
                    min_x = temp_spectrum;
                    min_sound = temp_sound;
                    min_alpha = alpha;
                    % �덷���X�V
                    min_amp_error = amp_error;
                    min_A_error = A_error;
                end
                    
            end
            
            spectrum = min_x;
            sound = min_sound;
            
            % ��
            fprintf('  min distance A filter : %d ,   min alpha : %d \n', min_A_error, min_alpha);
            
        end
        
        function [spectrum, sound, min_amp_error, min_A_error, min_alpha] = SDMM(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, signal_len, STFT_type, A_weight)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 24 May. 2019.
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
                
            % 1.e1000(����) + 1.e1000i(����) �̕��f����p�ӂ��C�X�y�N�g�������ŏ��ƂȂ�X�y�N�g���𓾂�
            temp_comp = complex(double(1.e1000), double(1.e1000));
            min_x = repmat(temp_comp, freq, frames);
            min_alpha = 0.0;
            min_amp_error = zeros(iteration,1);
            min_A_error = zeros(iteration,1);
            min_sound = zeros(signal_len,1);
            
            % alpha�̍X�V
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]

                % �����l
                %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
                %       y1, y2 =  x: SDMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�ϐ�y1,y2��p��
                %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
                %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
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
                     
                    % x���X�V
                    x = ( ( y1 - z1 ) + ( y2 - z2 ) ) / 2;
                    
                    %y1���X�V
                    %   ���z�I�ȐU��A�̏W���֎ʑ����邱�Ƃɂ��U�����C��
                    y1 = ( ( rho*amp_corr + ( 1 - alpha )*abs(x+z1) ) / ( 1 - alpha + rho ) ) .* exp( 1i * angle(x+z1) );
                    
                    %y2���X�V
                    %   �������Ƃ��� STFT( ISTFT() )�������Ȃ�
                    S1 = ISTFT(x+z2, windual, shiftsize, fftsize, signal_len, STFT_type);
                    S2 = STFT(S1, win, shiftsize, fftsize, STFT_type);
                    %   �X�y�N�g���O�����̒l��W���ւ̎ʑ��ɂ��ʑ����C��
                    y2 = ( alpha*(x+z2) + rho*S2 ) / ( alpha + rho );
                    
                    % z1,z2�̍X�V
                    z1 = z1 + x - y1;
                    z2 = z2 + x - y2;
                    
                    % �C�e���[�V�������ɕ]�����邽�߂Ɉ�xSTFT�̒l���ԂɎˉe
                    temp_sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                    temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                    
                    % �덷���o��
                    amp_error(i) = norm(abs(temp_spectrum)-amp_corr,'fro');
                    A_error(i) = norm(diag(A_weight)*(abs(temp_spectrum)-amp_corr),'fro');
                                   
                end
                
                % �����X�y�N�g���Ƃ̍����݂�
                %temp_err_x = norm(spectrum_corr - x, 'fro');
                temp_err_x = A_error(end);
                
                % �����̈�
                %fprintf('    distance : %d \n', temp_err_x);
                
                % ������������΍X�V
                %if norm(diag(A_weight)*(abs(min_x)-amp_corr),'fro') > temp_err_x
                if min_A_error(end) > temp_err_x
                    min_x = temp_spectrum;
                    min_sound = temp_sound;
                    min_alpha = alpha;
                    % �덷���X�V
                    min_amp_error = amp_error;
                    min_A_error = A_error;
                end
                    
            end
            
            spectrum = min_x;
            sound = min_sound;
            
            % ��
            fprintf('  min distance A filter : %d ,   min alpha : %d \n', min_A_error, min_alpha);
            
        end
        
        function [spectrum, sound, min_amp_error, min_A_error, min_alpha] = PPXA(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, spectrum_corr, frames, signal_len, gamma, STFT_type, A_weight)
            %
            % Corded by R.Nakatsu (is0269rx@ed.ritsumei.ac.jp) on 25 May. 2019.
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
                
            % 1.e1000(����) + 1.e1000i(����) �̕��f����p�ӂ��C�X�y�N�g�������ŏ��ƂȂ�X�y�N�g���𓾂�
            temp_comp = complex(double(1.e1000), double(1.e1000));
            min_x = repmat(temp_comp, freq, frames);
            min_alpha = 0.0;
            min_amp_error = zeros(iteration,1);
            min_A_error = zeros(iteration,1);
            min_sound = zeros(signal_len,1);
            
            % �d�݃ւ�����
            w1 = 1/2;
            w2 = 1 - w1;
            
            % alpha�̍X�V
            %for alpha = 0.05:0.05:0.95
            %for alpha = [0.1, 0.15, 0.2, 0.4, 0.8]
            for alpha = [0.05, 0.07, 0.1, 0.12, 0.15, 0.2]

                % �����l
                %       x = amp_corr .* exp(1i * phase_temp) : ���]�̐U���ƃ����_���Ȉʑ��ɂ��X�y�N�g��
                %       y1, y2 =  x: SDMM�ŉ����œK�����ɗ��Ƃ����ނ��߁C�ϐ�y1,y2��p��
                %       amp_error = 0 : �U���Ԃ̃t���x�j�E�X�m����
                %       A_error = 0 : A�����t���̐U���Ԃ̃t���x�j�E�X�m����
                x = amp_corr .* exp(1i * phase_temp);
                y1 = x;
                y2 = x;
                amp_error = zeros(iteration,1);
                A_error = zeros(iteration,1);
                
                %%%%%%%%%%%%%%%%%%%%
                % APXX
                %%%%%%%%%%%%%%%%%%%%
                
                for i = 1:iteration
                     
                    % p1���X�V
                    %   ���z�I�ȐU��A�̏W���֎ʑ����邱�Ƃɂ��U�����C��
                    p1 = ( ( rho*amp_corr + ( w1 - alpha*w1 )*abs(y1) ) / ( w1 - alpha*w1 + rho ) ) .* exp( 1i * angle(y1) );

                    % p2���X�V
                    %   �������Ƃ��� STFT( ISTFT() )�������Ȃ�
                    S1 = ISTFT(y2, windual, shiftsize, fftsize, signal_len, STFT_type);
                    S2 = STFT(S1, win, shiftsize, fftsize, STFT_type);
                    %   �X�y�N�g���O�����̒l��W���ւ̎ʑ��ɂ��ʑ����C��
                    p2 = ( alpha*w2*(y2) + rho*S2 ) / ( alpha*w2 + rho );
                    
                    % p�̍X�V
                    p = w1*p1 + w2*p2;
                    
                    % y1,y2�̍X�V
                    y1 = y1 + gamma*(2*p - x - p1);
                    y2 = y2 + gamma*(2*p - x - p2);
                    
                    % x���X�V
                    x = x + gamma*(p - x);
                    
                    % �C�e���[�V�������ɕ]�����邽�߂Ɉ�xSTFT�̒l���ԂɎˉe
                    temp_sound = ISTFT(x, windual, shiftsize, fftsize, signal_len, STFT_type);
                    temp_spectrum = STFT(temp_sound, win, shiftsize, fftsize, STFT_type);
                    
                    % �덷���o��
                    amp_error(i) = norm(abs(temp_spectrum)-amp_corr,'fro');
                    A_error(i) = norm(diag(A_weight)*(abs(temp_spectrum)-amp_corr),'fro');
                   
                end
                
                % �����X�y�N�g���Ƃ̍����݂�
                %temp_err_x = norm(spectrum_corr - x, 'fro');
                temp_err_x = A_error(end);
                
                % �����̈�
                %fprintf('    distance : %d \n', temp_err_x);
                
                % ������������΍X�V
                %if norm(diag(A_weight)*(abs(min_x)-amp_corr),'fro') > temp_err_x
                if min_A_error(end) > temp_err_x
                    min_x = temp_spectrum;
                    min_sound = temp_sound;
                    min_alpha = alpha;
                    % �덷���X�V
                    min_amp_error = amp_error;
                    min_A_error = A_error;
                end
                    
            end
            
            spectrum = min_x;
            sound = min_sound;
            
            % ��
            fprintf('  min distance A filter : %d ,   min alpha : %d \n', min_A_error, min_alpha);
            
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




