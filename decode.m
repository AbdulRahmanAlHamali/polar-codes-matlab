function result = decode(received, frozen_indices, frozen_bits, channel_type, param)
  length = size(received, 1);
  %We start by initializing storage elements to hold our intermediate
  %results
  likelihoods = cell(log2(length) + 1);
  xors = cell(log2(length));
  for i = 1:(log2(length)+1)
    likelihoods{i} = NaN(length/(2^(i-1)), 2^(i-1));
  end
  
  for i = 1:log2(length)
    xors{i} = NaN(length/(2^(i-1)), 2^(i-1));
  end
  
  for i = 1:size(frozen_indices, 1)
    frozen_indices(i) = reverse_index(frozen_indices(i), log2(length));
  end
  
  u = NaN(size(received, 1), 1);
  for i = 1:size(received, 1)
    [l, likelihoods] = calculate_likelihood(received, likelihoods, xors, 1, i, 1, log2(length) + 1, channel_type, param);
    frozen = find(frozen_indices == i);
    if (~isempty(frozen))
      u(i) = frozen_bits(frozen(1));
    else
      if (l > 1)
        u(i) = 0; 
      else
        u(i) = 1;
      end
    end
    format short g
    disp([i, l, u(i)]);
    xors = propagate_value(u(i), xors, 1, i, 1, log2(length) + 1); 
  end
  result = bitrevorder(u);
end

function [likelihood, storage] = calculate_likelihood(received, likelihoods, xors, layer, idx1, idx2, max_depth, channel_type, param)
  if (~isnan(likelihoods{layer}(idx1, idx2)))
    likelihood = likelihoods{layer}(idx1, idx2);
    storage = likelihoods;
    return;
  end
  
  if (layer == max_depth)
    l = get_channel_likelihood(received(reverse_index(idx2, log2(size(received, 1)))), channel_type, param);
    likelihood = l;
    likelihoods{layer}(idx1, idx2) = l;
    storage = likelihoods;
    return;
  end
  
  [l1, likelihoods] = calculate_likelihood(received, likelihoods, xors, layer+1, (floor((idx1 - 1)/2)) + 1, (2*(idx2 - 1)) + 1, max_depth, channel_type, param);
  [l2, likelihoods] = calculate_likelihood(received, likelihoods, xors, layer+1, (floor((idx1 - 1)/2)) + 1, (2*(idx2 - 1) + 1) + 1, max_depth, channel_type, param);
  if (mod(idx1 - 1, 2) == 0)
    if (strcmp(channel_type, 'bec'))
      if ((isinf(l1) && isinf(l2)) || (l1 == 0 && l2 == 0))
        likelihood = Inf;
      elseif ((isinf(l1) && l2 == 0) || (l1 == 0 && isinf(l2)))
        likelihood = 0;
      else
        likelihood = 1;
      end
    else
      likelihood = (l1*l2 + 1)/(l1 + l2);
    end
  else
    if (strcmp(channel_type, 'bec'))
      if (l2 ~= 1)
        likelihood = l2;
      elseif (l1 ~= 1)
        likelihood = xor(l1, xors{layer}(idx1 - 1, idx2)); 
      else
        likelihood = 1;
      end
    else
      if (xors{layer}(idx1 - 1, idx2) == 0)
        likelihood = l1*l2;    
      else
        likelihood = l2/l1;
      end
    end
  end
  likelihoods{layer}(idx1, idx2) = likelihood;
  storage = likelihoods;
end

function likelihood = get_channel_likelihood(value, channel_type, param)
  if (strcmp(channel_type, 'bsc'))
    likelihood = get_bsc_likelihood(value, param);
  elseif (strcmp(channel_type, 'bec'))
    likelihood = get_bec_likelihood(value);
  elseif (strcmp(channel_type, 'awgn'))
    likelihood = get_awgn_likelihood(value, param);
  end
end

function likelihood = get_bsc_likelihood(value, flip_probability)
  if (value == 0)
    likelihood = (1-flip_probability)/flip_probability; 
  else
    likelihood = flip_probability/(1-flip_probability);
  end
end

function likelihood = get_bec_likelihood(value)
  if (value == 0)
    likelihood = Inf;
  elseif (value == 1)
    likelihood = 0;
  else
    likelihood = 1;
  end
end

function likelihood = get_awgn_likelihood(value, sigma)
  likelihood = exp((-(value+1)^2/(2*sigma^2)))/exp((-(value-1)^2/(2*sigma^2)));
end

function reversed_index = reverse_index(idx, n)
  reversed_index = bin2dec(fliplr(dec2bin(idx - 1,n))) + 1;
end

function new_values = propagate_value(value, xors, layer, idx1, idx2, max_depth)
  xors{layer}(idx1, idx2) = value;
  
  if (layer == max_depth)
    new_values = xors;
    return;
  end
  
  if (mod(idx1 - 1, 2) == 0)
    new_values = xors;
    return;
  end
  
  xors = propagate_value(xor(value, xors{layer}(idx1-1, idx2)), xors, layer+1, (floor((idx1 - 1)/2)) + 1, 2*(idx2 - 1) + 1, max_depth);
  xors = propagate_value(value, xors, layer+1, (floor((idx1 - 1)/2)) + 1, (2*(idx2 - 1) + 1) + 1, max_depth);
  
  new_values = xors;
end