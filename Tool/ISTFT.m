function sigr = ISTFT(C,windual,skip,winLen,signal_len,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             inverse STFT
%
%%% -- Input --------------------------------------------------------------
% C     : spectrograms (freq x time)
% windual   : synthesis window (winLen x 1)
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
% sigr  : signal (samples x 1)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !! Ls must be even number due to our STFT/iSTFT implementation !!
Ls = ceil((signal_len+2*(winLen-skip)-winLen)/skip)*skip+winLen;
if flag == 2
    hWL = floor(winLen/2);
    C = C.*exp(+2i*pi*(mod((0:hWL)'*(0:size(C,2)-1)*skip,winLen)/winLen)); % 通常のiSTFTならばここをコメントアウト
end
sigr = ifft([C;zeros(size(C)-[2,0])],'symmetric').*windual;
idx = (1:winLen)' + (0:skip:Ls-winLen);
idx2 = repmat(1:size(C,2),winLen,1);
sigr = full(sum(sparse(idx(:),idx2(:),sigr(:)),2));
sigr = sigr(winLen-skip+1:winLen-skip+signal_len); % 信号の両端に0を埋めた部分を取り除く
end