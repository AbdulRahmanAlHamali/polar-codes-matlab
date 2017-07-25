number_of_tests = 500;
N = 128;
frozen_indices = transpose(setdiff(1:128, [48 56 60 62:64 80 88 90:96 102:104 106:112 114:128]));
K = N - size(frozen_indices, 1);
frozen_bits = zeros(size(frozen_indices, 1), 1);
channel_type = 'awgn';
param = 0.987;
list_size = 8;

sc_count = 0;
lsc_count = 0;
for i = 1:number_of_tests
  message = randi([0 1], K, 1);
  encoded_message = encode(message, frozen_indices, frozen_bits);
  received_message = introduce_noise(channel_type, encoded_message, param);
  sc_result = decode(received_message, frozen_indices, frozen_bits, channel_type, param);
  l_sc_result = list_decode(received_message, frozen_indices, frozen_bits, channel_type, param, list_size);
  sc_result(frozen_indices) = [];
  l_sc_result(frozen_indices) = [];
  if (sc_result == message)
    sc_count = sc_count + 1;
  end
  if (l_sc_result == message)
    lsc_count = lsc_count + 1;
  end
end

disp('SC succeeded:')
disp(sc_count);
disp('LSC succeeded:')
disp(lsc_count);