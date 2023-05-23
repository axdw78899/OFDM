clc;
clear all;

%setting
nfft = 64;
n_fft = 64;
n_cpe = 16;
snr = 20;

%image
im=imread('image0-1.png');
im_bin = dec2bin(im(:))';
im_bin = im_bin(:);

%binary stream to symbols
sym_rem = mod(2-mod(length(im_bin),2),2);
padding = repmat('0',sym_rem,1);
im_bin_padded = [im_bin;padding];

cons_data = reshape(im_bin_padded,2,length(im_bin_padded)/2)';
cons_sym_id = bin2dec(cons_data);

%QPSK symbol book
symbol_book=[0.7071-0.7071i; -0.7071-0.7071i;-0.7071+0.7071i;0.7071+0.7071i];
X = symbol_book(cons_sym_id+1);

%IFFT
fft_rem = mod(n_fft-mod(length(X),n_fft),n_fft);
X_padded = [X;zeros(fft_rem,1)];
X_blocks = reshape(X_padded,nfft,length(X_padded)/nfft);
x = ifft(X_blocks);

%CP
x_cpe = [x(end-n_cpe+1:end,:);x];
x_s   = x_cpe(:);

%AWGN
data_pwr = mean(abs(x_s.^2));
noise_pwr = data_pwr/10^(snr/10);
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;
x_s_noise = x_s + noise;
snr_meas = 10*log10(mean(abs(x_s.^2))/mean(abs(noise.^2)));

%apply fading channel
g = exp(-(0:7));
g = g/norm(g);
x_s_noise_fading = conv(x_s_noise,g,'same');

%fft 
x_p = reshape(x_s_noise_fading,nfft+n_cpe,length(x_s_noise_fading)/(nfft+n_cpe));
x_p_cpr = x_p(n_cpe+1:end,:);
X_hat_blocks = fft(x_p_cpr);

%estimate
G = X_hat_blocks(:,1)./X_blocks(:,1);
X_hat_blocks = X_hat_blocks./repmat(G,1,size(X_hat_blocks,2));

%remove padding
X_hat = X_hat_blocks(:);
X_hat = X_hat(1:end-fft_rem);

%recover
rec_syms = knnsearch([real(symbol_book) imag(symbol_book)],[real(X_hat) imag(X_hat)])-1;
rec_syms_cons = dec2bin(rec_syms);

rec_im_bin = reshape(rec_syms_cons',numel(rec_syms_cons),1);
rec_im_bin = rec_im_bin(1:end-sym_rem);

ber = sum(abs(rec_im_bin-im_bin))/length(im_bin);

rec_im = reshape(rec_im_bin,8,numel(rec_im_bin)/8);
rec_im = uint8(bin2dec(rec_im'));
rec_im = reshape(rec_im,size(im));

%plot
plot(X,'x','linewidth',2,'markersize',10);
xlim([-2 2]);
ylim([-2 2]);
xlabel('In phase')
ylabel('Qudrature')
imshow(rec_im);
