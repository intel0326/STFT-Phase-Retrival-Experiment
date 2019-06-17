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