function crc = calculate_crc(message, polynomial, initial_crc)
  r = size(polynomial, 1) - 1;
  
  appended_message = [message;initial_crc];
  
  dividend = appended_message(1:r + 1, 1);
  for i = 1:size(message, 1)
    if (dividend(1) == 1)
      xor_result = xor(dividend, polynomial);
    else
      xor_result = dividend;
    end
    if (i ~= size(message, 1))
      dividend = [xor_result(2:r+1, 1); appended_message(r+1+i, 1)];
    else
      dividend = xor_result;
    end
  end
  
  crc = dividend(2:r+1);
end