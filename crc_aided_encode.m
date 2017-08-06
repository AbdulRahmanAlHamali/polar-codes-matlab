%This function applies crc-aided polar encoding on a given message. 
%message: ((K-r)x1) matrix, represents the information bits
%crc_polynomial: (rx1) matrix, represents the CRC polynomial
%frozen_indices: ((N-K)x1) matrix, specifies which bits are frozen
%frozen_bits: ((N-K)x1) matrix specifies the values of the frozen bits
function codeword = crc_aided_encode(message, crc_polynomial, frozen_indices, frozen_bits)
  %calculate the crc
  crc = calculate_crc(message, crc_polynomial, zeros(size(crc_polynomial, 1) - 1, 1));
  appended_message = [message;crc];
  
  %Prepare the generator matrix
  generator_matrix = get_classic_generator_matrix(size(appended_message, 1) + size(frozen_bits, 1));
  
  %Prepare the message to be encoded
  u(1: size(appended_message, 1) + size(frozen_bits, 1), 1) = NaN;
  u(frozen_indices) = frozen_bits;
  u(isnan(u)) = appended_message; 
  
  %Multiply the message with the generator matrix in GF(2)
  codeword = transpose(u) * generator_matrix;
  codeword = transpose(mod(codeword, 2)); 
end 

%The classic polar coding generator matrix (with no CRC bits), is acquired
%by taking the kronecker power of [1 0; 1 1], and then permuting its
%columns in a bit-reversal order
function G = get_classic_generator_matrix(sz)
  F = [1 0; 1 1];
  K = F;
  for i = 1 : log(sz)/log(2) - 1    
    K = kron(K, F);
  end
  G = bitrevorder(K);
end

