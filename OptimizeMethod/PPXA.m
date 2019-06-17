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