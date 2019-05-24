
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
total_sec = 10;
% �Z���ԃt�[���G�ϊ��̃t���[���t���[����
fftsize = 1024;
% �Z���ԃt�[���G�ϊ��̃t���[���V�t�g��
shiftsize = 256;
% ADMM�̃C�e���[�V�����񐔂��w��
iteration = 1000;
% ���̎��
win = hann(fftsize,'periodic'); % �n�j���O��
% �Ώۉ�
filename = './Sound_source/mixture.wav';
% �o�͐�
outputDir = './Output';
% Douglas-Rachford Splitting Algorithm �̃��̒l
gamma = 1.9;
% �G�N�Z���V�[�g�̏����l
sell_angle = 3;
sell_spe = 16;


%%%%%%%%%%%%%%%%%%%%
% �O����
%%%%%%%%%%%%%%%%%%%%

% 1.ISTFT�ɗ��p����t�̑�������
% 2.�^�̕��f�X�y�N�g���O�������擾
[windual, spectrum, Ls, signal_len] = ins_tool.AudioReadMethod(filename, total_sec, freq, fftsize, shiftsize, win);

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

