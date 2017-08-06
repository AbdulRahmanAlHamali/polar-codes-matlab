number_of_tests = 500;
N = 128;
K = 40;
channel_type = 'awgn';
param = 0.987;
list_size = 16;
crc_polynomial = [1;1;1;0;1;0;1;0;1];
r = size(crc_polynomial, 1) - 1;
k = K - r;

z = calculate_z_parameters(N, 1000, channel_type, param);
[s, ind] = sort(z);
frozen_indices = ind(K+1:N);
%frozen_indices = transpose(setdiff(1:128, [48 56 60 62:64 80 88 90:96 102:104 106:112 114:128]));
%K = N - size(frozen_indices, 1);
frozen_bits = zeros(size(frozen_indices, 1), 1);

sc_count = 0;
lsc_count = 0;
crc_lsc_count = 0;
no_crc = 0;
for i = 1:number_of_tests
  message = randi([0 1], k, 1);
  crc_appended_message = [message;calculate_crc(message, crc_polynomial, zeros(size(crc_polynomial, 1) - 1, 1))];
  crc_encoded_message = crc_aided_encode(message, crc_polynomial, frozen_indices, frozen_bits);
  received_message = introduce_noise(channel_type, crc_encoded_message, param);
  sc_result = decode(received_message, frozen_indices, frozen_bits, channel_type, param);
  l_sc_result = list_decode(received_message, frozen_indices, frozen_bits, channel_type, param, list_size);
  crc_l_sc_result = crc_aided_list_decode(received_message, crc_polynomial, frozen_indices, frozen_bits, channel_type, param, list_size);
  sc_result(frozen_indices) = [];
  l_sc_result(frozen_indices) = [];
  crc_l_sc_result(frozen_indices) = [];

  if (sc_result == crc_appended_message)
    sc_count = sc_count + 1;
  end
  if (l_sc_result == crc_appended_message)
    lsc_count = lsc_count + 1;
    disp('LSC!!');
  end
  if (crc_l_sc_result == crc_appended_message)
    crc_lsc_count = crc_lsc_count + 1;
    disp('CRC_LSC!!');
  elseif (crc_l_sc_result(1) == -1)
    no_crc = no_crc + 1;
  end
end

disp('SC succeeded:')
disp(sc_count);
disp('LSC succeeded:')
disp(lsc_count);
disp('CRC_LSC succeeded:')
disp(crc_lsc_count);
disp('Instances where no one had a correct CRC');
disp(no_crc);