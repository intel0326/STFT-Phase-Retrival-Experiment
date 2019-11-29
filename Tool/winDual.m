function windual = winDual(win,skip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                  Calculating the canonical dual window
%
%   升山氏のSTFTプログラムを拝借，一部改編しております
%   https://ieeexplore.ieee.org/document/8552369/algorithms#algorithms
%
%
%%% -- Input --------------------------------------------------------------
% win    : analysis window (winLen x 1)
% skip   : skipping samples (1 x 1)
%
%%% -- Output -------------------------------------------------------------
% windual: synthesis window (winLen x 1)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

windual = [win; zeros(skip*ceil(length(win)/skip)-length(win),1)];
windual = reshape(windual,skip,[]);
windual = windual./sum(abs(windual).^2,2);
windual = reshape(windual(1:length(win)),[],1);
end