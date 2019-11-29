function [spectrum, sound, amp_error, A_error] = Prop(amp_corr, rho, fftsize, shiftsize, win, windual, iteration, phase_temp, freq, frames, signal_len, STFT_type, A_weight)
    %
    % Corded by R.Nakatsu (dragonstar30210326@gmail.com) on 11 May. 2019.
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