%This function simulates sending a message through a channel that
%introduces noise to it.
%channel_type: Could be 'bsc' for Binary Symmetric Channel, 'bec' for
%Binary Erasure Channel, and 'awgn' for Additive White Gaussian Noise
%Channel
%message: (Nx1) matrix, represents the message being sent through the
%channel
%param: in 'bsc', param specifies the probability of the bit getting
%flipped, in 'bec', it represents the probability of erasure, and in 'awgn'
%it represents the standard deviation of the distribution
function result = introduce_noise(channel_type, message, param)
  if (strcmp(channel_type, 'bsc'))
    result = introduce_bsc_noise(message, param);
  elseif (strcmp(channel_type, 'bec'))
    result = introduce_bec_noise(message, param);
  elseif (strcmp(channel_type, 'awgn'))
    result = introduce_awgn_noise(message, param);
  end
end

function result = introduce_bsc_noise(message, flip_probability)
  result = zeros(size(message, 1), 1);
  for i = 1:size(message, 1)
    p = rand;
    if (p <= flip_probability)
      result(i) = ~message(i);    
    else
      result(i) = message(i);
    end
  end
end

%In BEC, we represent and erasure by -1.
function result = introduce_bec_noise(message, flip_probability)
  result = zeros(size(message, 1), 1);
  for i = 1:size(message, 1)
    p = rand;
    if (p <= flip_probability)
      result(i) = -1;    
    else
      result(i) = message(i);
    end
  end
end

%By convention, BAWGNC takes two values as input, -1 and 1. We could make it 0 and 1 but we keep it as is for convention's sake. 
%Note: in this channel, the received results are samples from two normal
%distributions whose mean is the actual value of the bit, and whose
%standard deviation is sigma
function result = introduce_awgn_noise(message, sigma)
  % First, we transform all 0 to -1, for convention's sake
  message(message == 0) = -1;
  result = normrnd(message, sigma);
end