function C = STFT(sig,win,skip,winLen,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                                  STFT
%
%   升山氏のSTFTプログラムを拝借，一部改編しております
%   https://ieeexplore.ieee.org/document/8552369/algorithms#algorithms
%
%
%%% -- Input --------------------------------------------------------------
% sig   : signal (samples x 1)
% win   : analysis window (winLen x 1)
% skip  : skipping samples (1 x 1)
% winLen: window length (1 x 1)
% flag: 1: standard STFT, 2: special STFT
% Ls    : signal length (1 x 1)
% signal_len : Original signal length (1 x 1)
%
% !! Attention !!
% length(sig) = Ls = win + skip x N
% where N in natural numbers 
%
%
%%% -- Output -------------------------------------------------------------
% C     : spectrograms (freq x time)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal_len = length(sig);
% !! Ls must be even number due to our STFT/iSTFT implementation !!
Ls = ceil((signal_len+2*(winLen-skip)-winLen)/skip)*skip+winLen;
% zero パディング：信号の両端を0詰め
sig = [zeros(winLen-skip,1); sig; zeros(Ls-signal_len-2*(winLen-skip),1); zeros(winLen-skip,1)];

idx = (1:winLen)' + (0:skip:Ls-winLen);
size(idx);
C = fft(sig(idx).*win);
hWL = floor(winLen/2);
if flag == 1
    C = C(1:hWL+1,:); % 通常のSTFT
elseif flag == 2
    C = C(1:hWL+1,:).*exp(-2i*pi*(mod((0:hWL)'*(0:size(C,2)-1)*skip,winLen)/winLen)); % 理論的に好ましいSTFT
end
end