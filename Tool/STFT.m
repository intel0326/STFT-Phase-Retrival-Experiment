function C = STFT(sig,win,skip,winLen,Ls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                                  STFT
%
%%% -- Input --------------------------------------------------------------
% sig   : signal (samples x 1)
% win   : analysis window (winLen x 1)
% skip  : skipping samples (1 x 1)
% winLen: window length (1 x 1)
% Ls    : signal length (1 x 1)
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

idx = (1:winLen)' + (0:skip:Ls-winLen);
size(idx);
C = fft(sig(idx).*win);
hWL = floor(winLen/2);
%C = C(1:hWL+1,:); % í èÌÇÃSTFT
C = C(1:hWL+1,:).*exp(-2i*pi*(mod((0:hWL)'*(0:size(C,2)-1)*skip,winLen)/winLen)); % óùò_ìIÇ…çDÇ‹ÇµÇ¢STFT
end