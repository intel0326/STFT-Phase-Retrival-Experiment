function waveform = ISTFT( S, shiftSize, window, length )
%
% Corded by D. Kitamura (d-kitamura@nii.ac.jp) on 9 Aug. 2016.
%
% This function calculates inverse short-time Fourier transform (ISTFT) of
% input complex-valued spectrogram. Both single and multi-channel signals
% are supported.
%
% see also
% http://d-kitamura.sakura.ne.jp/index.html
%
% [syntax]
%   waveform = ISTFT(S)
%   waveform = ISTFT(S, shiftSize)
%   waveform = ISTFT(S, shiftSize, window)
%   waveform = ISTFT(S, shiftSize, window, length)
%
% [inputs]
%           S: STFT of input signal (fftSize/2+1 x frames x nch)
%   shiftSize: frame shift (default: fftSize/2)
%      window: window function used in STFT (fftSize x 1) or choose used
%              function from below.
%              'hamming'    : Hamming window (default)
%              'hann'       : von Hann window
%              'rectangular': rectangular window
%              'blackman'   : Blackman window
%              'sine'       : sine window
%      length: length of original signal (before STFT) (default: the same with output signal)
%
% [outputs]
%   waveform: time-domain waveform of the input spectrogram (length x nch)
%

% Check errors and set default values
if (nargin < 1)
    error('Too few input arguments.\n');
end

[freq, frames, nch] = size(S);

if (nch > freq)
    error('The input spectrogram might be wrong. The size of it must be (freq x frame x ch).\n');
end
if (isreal(S) == 1)
    error('The input spectrogram might be wrong. It does not complex-valued matrix.\n');
end
if (mod(freq,2) == 0)
    error('The number of rows of the first argument must be an odd number because it is (frame length /2)+1.\n');
end

fftSize = (freq-1) * 2;

if (nargin < 2)
    shiftSize = fftSize/2;
elseif (mod(fftSize,shiftSize) ~= 0)
    error('The frame length must be dividable by the second argument.\n');
end
if (nargin<3)
    window = hamming_local(fftSize); % use default window
else
    if (isnumeric(window))
        if (size(window, 1) ~= fftSize)
            error('The length of the third argument must be the same as FFT size used in STFT.\n');
        end
    else
        switch window
            case 'hamming'
                window = hamming_local(fftSize);
            case 'hann'
                window = hann_local(fftSize);
            case 'rectangular'
                window = rectangular_local(fftSize);
            case 'blackman'
                window = blackman_local(fftSize);
            case 'sine'
                window = sine_local(fftSize);
            otherwise
                error('Unsupported window is required. Type "help STFT" and check options.\n')
        end
    end
end
invWindow = optSynWnd_local( window, shiftSize );

% Inverse STFT
spectrum = zeros(fftSize,1);
tmpSignal = zeros((frames-1)*shiftSize+fftSize,nch);
for ch = 1:nch
    for i = 1:3
        spectrum(1:fftSize/2+1,1) = S(:,i,ch);
        spectrum(1,1) = spectrum(1,1)/2;
        spectrum(fftSize/2+1,1) = spectrum(fftSize/2+1,1)/2;
        sp = (i-1)*shiftSize;
        tmpSignal(sp+1:sp+fftSize,ch) = tmpSignal(sp+1:sp+fftSize,ch) + real((ifft(spectrum,fftSize)).*invWindow(:,:,i))*2;
    end
    for i = 4:frames-3
        spectrum(1:fftSize/2+1,1) = S(:,i,ch);
        spectrum(1,1) = spectrum(1,1)/2;
        spectrum(fftSize/2+1,1) = spectrum(fftSize/2+1,1)/2;
        sp = (i-1)*shiftSize;
        tmpSignal(sp+1:sp+fftSize,ch) = tmpSignal(sp+1:sp+fftSize,ch) + real((ifft(spectrum,fftSize)).*invWindow(:,:,4))*2;
    end
    ii = 5;
    for i = frames-2:frames
        spectrum(1:fftSize/2+1,1) = S(:,i,ch);
        spectrum(1,1) = spectrum(1,1)/2;
        spectrum(fftSize/2+1,1) = spectrum(fftSize/2+1,1)/2;
        sp = (i-1)*shiftSize;
        tmpSignal(sp+1:sp+fftSize,ch) = tmpSignal(sp+1:sp+fftSize,ch) + real((ifft(spectrum,fftSize)).*invWindow(:,:,ii))*2;
        ii = ii + 1;
    end
    
end
waveform = tmpSignal;%(fftSize-shiftSize+1:(frames-1)*shiftSize+fftSize,:); % discard padded zeros in STFT

% Discarding padded zeros in the end of the signal
if (nargin==4)
    waveform = waveform(1:length,:);
end
end

%% Local functions
function window = hamming_local(fftSize)
t = linspace(0,1,fftSize+1).'; % periodic (produce L+1 window and return L window)
window = 0.54*ones(fftSize,1) - 0.46*cos(2.0*pi*t(1:fftSize));
end

function window = hann_local(fftSize)
t = linspace(0,1,fftSize+1).'; % periodic (produce L+1 window and return L window)
window = max(0.5*ones(fftSize,1) - 0.5*cos(2.0*pi*t(1:fftSize)),eps);
end

function window = rectangular_local(fftSize)
window = ones(fftSize,1);
end

function window = blackman_local(fftSize)
t = linspace(0,1,fftSize+1).'; % periodic (produce L+1 window and return L window)
window = max(0.42*ones(fftSize,1) - 0.5*cos(2.0*pi*t(1:fftSize)) + 0.08*cos(4.0*pi*t(1:fftSize)),eps);
end

function window = sine_local(fftSize)
t = linspace(0,1,fftSize+1).'; % periodic (produce L+1 window and return L window)
window = max(sin(pi*t(1:fftSize)),eps);
end

function synthesizedWindow = optSynWnd_local(analysisWindow,shiftSize)
fftSize = size(analysisWindow,1);
synthesizedWindow = zeros(fftSize,1,7);

for i = 1:shiftSize
    amp = 0.0;
    for j = 1:fftSize/shiftSize
        amp = amp + analysisWindow(i+(j-1)*shiftSize,1)^2;
    end
    for j = 1:fftSize/shiftSize
        synthesizedWindow(i+(j-1)*shiftSize,1,4) = analysisWindow(i+(j-1)*shiftSize,1)/amp;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:shiftSize
    synthesizedWindow(i,1,1) = 1/analysisWindow(i,1);
end
for i = shiftSize+1:2*shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i-shiftSize,1)^2;
    synthesizedWindow(i,1,1) = analysisWindow(i,1)/amp;
end
for i = 2*shiftSize+1:3*shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i-shiftSize,1)^2+analysisWindow(i-2*shiftSize,1)^2;
    synthesizedWindow(i,1,1) = analysisWindow(i,1)/amp;
end  
for i = 3*shiftSize+1:4*shiftSize
    synthesizedWindow(i,1,1) = synthesizedWindow(i,1,4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i+shiftSize,1)^2;
    synthesizedWindow(i,1,2) = analysisWindow(i,1)/amp;
end
for i = shiftSize+1:2*shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i-shiftSize,1)^2+analysisWindow(i+shiftSize,1)^2;
    synthesizedWindow(i,1,2) = analysisWindow(i,1)/amp;
end 
for i = 2*shiftSize+1:4*shiftSize
    synthesizedWindow(i,1,2) = synthesizedWindow(i,1,4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i+shiftSize,1)^2+analysisWindow(i+2*shiftSize,1)^2;
    synthesizedWindow(i,1,3) = analysisWindow(i,1)/amp;
end
for i = shiftSize+1:4*shiftSize
    synthesizedWindow(i,1,3) = synthesizedWindow(i,1,4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:3*shiftSize
    synthesizedWindow(i,1,5) = synthesizedWindow(i,1,4);
end
for i = 3*shiftSize+1:4*shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i-shiftSize,1)^2+analysisWindow(i-2*shiftSize,1)^2;
    synthesizedWindow(i,1,5) = analysisWindow(i,1)/amp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:2*shiftSize
    synthesizedWindow(i,1,6) = synthesizedWindow(i,1,4);
end
for i = 2*shiftSize+1:3*shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i-shiftSize,1)^2+analysisWindow(i+shiftSize,1)^2;
    synthesizedWindow(i,1,6) = analysisWindow(i,1)/amp;
end
for i = 3*shiftSize+1:4*shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i-shiftSize,1)^2;
    synthesizedWindow(i,1,6) = analysisWindow(i,1)/amp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:shiftSize
    synthesizedWindow(i,1,7) = synthesizedWindow(i,1,4);
end
for i = shiftSize+1:2*shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i+shiftSize,1)^2+analysisWindow(i+2*shiftSize,1)^2;
    synthesizedWindow(i,1,7) = analysisWindow(i,1)/amp;
end
for i = 2*shiftSize+1:3*shiftSize
    amp = analysisWindow(i,1)^2+analysisWindow(i+shiftSize,1)^2;
    synthesizedWindow(i,1,7) = analysisWindow(i,1)/amp;
end
for i = 3*shiftSize+1:4*shiftSize
    synthesizedWindow(i,1,7) = 1/analysisWindow(i,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%