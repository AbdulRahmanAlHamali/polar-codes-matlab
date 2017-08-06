function z = calculate_z_parameters(block_length, sample_size, channel_type, param)

  z = zeros(block_length, 1);
  for i = 1:sample_size
    message = randi([0 1], block_length, 1);
    encoded_message = encode([], transpose(1:block_length), message);
    received_message = introduce_noise(channel_type, encoded_message, param);
    z_instance = calculate_z_instance(received_message, message, channel_type, param);
    z = z + z_instance;
  end
  
  z = z / sample_size;
end

function result = calculate_z_instance(received, frozen_bits, channel_type, param)
  length = size(received, 1);
  
  %We start by initializing storage elements to hold our intermediate
  %results. One of them will hold intermediate likelihoods of elements in
  %the middle to be a 0. The other holds the value of intermediate xor
  %operations that have already been calculated
  likelihoods = cell(log2(length) + 1);
  xors = cell(log2(length));
  for i = 1:(log2(length)+1)
    likelihoods{i} = NaN(length/(2^(i-1)), 2^(i-1));
  end
  
  for i = 1:log2(length)
    xors{i} = NaN(length/(2^(i-1)), 2^(i-1));
  end
  
  u = NaN(size(received, 1), 1);
  %we go here bit by bit, calculate the likelihood of the bit being a zero.
  %Then, we use its likelihood to calculate the Z parameter, and then we set it to its frozen value. 
  %Finally, we propagate the decision that we have made in order to help
  %coming levels calculate their values
  z_instance = NaN(length, 1);
  for i = 1:size(received, 1)
    [l, likelihoods] = calculate_likelihood(received, likelihoods, xors, 1, i, 1, log2(length) + 1, channel_type, param);
    
    u(i) = frozen_bits(reverse_index(i, log2(length)));
    
    if (u(i) == 1)
      z_instance(i) = sqrt(l);
    else
      z_instance(i) = sqrt(1/l);
    end
    
    xors = propagate_value(u(i), xors, 1, i, 1, log2(length) + 1); 
  end
  %When we are done, we have been doing everything in bit-reversed order,
  %so we reverse again to restore the original order
  result = bitrevorder(z_instance);
end

%calculates the likelihood of the bit at a specific level to be a zero, and
%updates the likelihood table accordingly
function [likelihood, storage] = calculate_likelihood(received, likelihoods, xors, layer, idx1, idx2, max_depth, channel_type, param)
  %if a value has already been calculated, return
  if (~isnan(likelihoods{layer}(idx1, idx2)))
    likelihood = likelihoods{layer}(idx1, idx2);
    storage = likelihoods;
    return;
  end
  
  %If we have reached the channel level, get the likelihood from the
  %channel and return
  if (layer == max_depth)
    l = get_channel_likelihood(received(reverse_index(idx2, log2(size(received, 1)))), channel_type, param);
    likelihood = l;
    likelihoods{layer}(idx1, idx2) = l;
    storage = likelihoods;
    return;
  end
  
  [l1, likelihoods] = calculate_likelihood(received, likelihoods, xors, layer+1, (floor((idx1 - 1)/2)) + 1, (2*(idx2 - 1)) + 1, max_depth, channel_type, param);
  [l2, likelihoods] = calculate_likelihood(received, likelihoods, xors, layer+1, (floor((idx1 - 1)/2)) + 1, (2*(idx2 - 1) + 1) + 1, max_depth, channel_type, param);
  %BEC requires special treatment because it haas zeroes and infinities
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

%propagates the xor values and updates the xor table accordingly
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