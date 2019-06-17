
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   �{�R�[�h�̐���
%       �uMain_program.m�v�Ŏg�p����ϐ��̏����l��ݒ�
%       �ϐ��̏����l�̐ݒ�ɔ����C�uinitialize.mat�v���쐬
%
%   ���s���@
%       �R�}���h�E�B���h�E���ŁuMain_program�v�����s
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% �W�{�����g��
freq = 16000;
% �������b�Ԏg����
% 6/17 40-50�b���g��
start_sec = 40;
end_sec = 50;
% �Z���ԃt�[���G�ϊ��̃t���[���t���[����
fftsize = 1024;
% �Z���ԃt�[���G�ϊ��̃t���[���V�t�g��
shiftsize = 256;
% ADMM�̃C�e���[�V�����񐔂��w��
iteration = 1000;
% ���̎��
win = hann(fftsize,'periodic'); % �n�j���O��
% �Ώۉ�
%filename = './Sound_source/mixture.wav';
filename = './Sound_source/01 Roundabout.mp3';
% �o�͐�
outputDir = './Output';
% STFT�̃^�C�v(��c������̒ʏ�̂��=1, ���z�I��STFT=2)
STFT_type = 1;
% A�����̃t�B���^���쐬
[~,A_weight] = filterA(ones(fftsize,1), freq);
% Douglas-Rachford Splitting Algorithm �̃��̒l
gamma = 1.0;
% �U���̔��U��h���d��weight �̏�����
Delta = 0.01;
% �G�N�Z���V�[�g�̏����l
sell_angle = 3;
sell_spe = 18;


%%%%%%%%%%%%%%%%%%%%
% �O����
%%%%%%%%%%%%%%%%%%%%

% 1.ISTFT�ɗ��p����t�̑�������
% 2.�^�̕��f�X�y�N�g���O�������擾
[windual, spectrum, signal_len, music] = ins_tool.AudioReadMethod(filename, start_sec, end_sec, freq, fftsize, shiftsize, win, STFT_type);

% ���]�̐U���ƈʑ����擾
amp_corr = abs(spectrum);
phase_corr = angle(spectrum);

% ���]�̐U��amp_corr����e��p�����[�^���擾
%       amp_FFTsize = STFT��̐܂�Ԃ����l�������CSTFT/2�̓_��
%           ��.8192�_�t�[���G / 2 = 4092�_
%       frames = STFT/2�̓_�������t���[�����邩
% �ʑ��̍s��T�C�Y���擾
[amp_FFTsize, frames] = size(amp_corr);
% �����_���Ȉʑ����擾
phase_temp = zeros(amp_FFTsize, frames);


save('./Variable/Initialize')

