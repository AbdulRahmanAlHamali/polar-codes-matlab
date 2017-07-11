function result = introduce_noise(channel_type, message, param)
  if (strcmp(channel_type, 'bsc'))
    result = introduce_bsc_noise(message, param);
  elseif (strcmp(channel_type, 'bec'))
    result = introduce_bec_noise(message, param);
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