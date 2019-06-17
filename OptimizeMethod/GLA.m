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