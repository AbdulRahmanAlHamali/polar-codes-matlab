function codeword = encode(message, frozen_indices, frozen_bits)
  %Prepare the generator matrix
  generator_matrix = get_classic_generator_matrix(size(message, 1) + size(frozen_bits, 1));
  
  %Prepare the message to be encoded
  u(1: size(message, 1) + size(frozen_bits, 1), 1) = -1;
  u(frozen_indices) = frozen_bits;
  u(u == -1) = message; 
  
  codeword = transpose(u) * generator_matrix;
  codeword = transpose(mod(codeword, 2)); 
end 

function G = get_classic_generator_matrix(sz)
  F = [1 0; 1 1];
  K = F;
  for i = 1 : log(sz)/log(2) - 1    
    K = kron(K, F);
  end
  G = bitrevorder(K);
end

