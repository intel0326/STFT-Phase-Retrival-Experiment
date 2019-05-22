
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
shiftsize = 128;
% ADMM�̃C�e���[�V�����񐔂��w��
iteration = 1000;
% ���̎��
window = 'hann';
% �Ώۉ�
filename = './Sound_source/mixture.wav';
% �o�͐�
outputDir = './Output';
% �G�N�Z���V�[�g�̏����l
sell_angle = 3;
sell_spe = 16;


%%%%%%%%%%%%%%%%%%%%
% �O����
%%%%%%%%%%%%%%%%%%%%

% 1.�����̓ǂݍ���
% 2.�^�̕��f�X�y�N�g���O�������擾
[music, spectrum] = ins_tool.AudioReadMethod(filename, total_sec, freq, fftsize, shiftsize, window);

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
phase_temp = rand(amp_FFTsize, frames);


save('./Variable/Initialize')

