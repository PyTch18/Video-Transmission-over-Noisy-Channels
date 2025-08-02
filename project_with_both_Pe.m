% Load video
V = VideoReader('highway.avi');
a = read(V);
frames = get(V,'NumberOfFrames');

% Initialize frame storage
ArrayOfFrames(frames) = struct('cdata', []);
ArrayOfFrames(1).cdata = a(:,:,:,1);  % Store first frame

% Get frame dimensions
s = size(ArrayOfFrames(1).cdata);

% Define coding parameters
error_probs = [0.01, 0.1];  % Both probabilities to test

% Define all trellis structures
trellises = struct(...
    'rate_1_2', poly2trellis(3, [5 7]), ...    % Rate 1/2 code
    'rate_1_3', poly2trellis(3, [5 7 7]), ...  % Rate 1/3 code
    'rate_1_4', poly2trellis(3, [5 7 7 7]), ... % Rate 1/4 code
    'rate_2_3', poly2trellis([2 2], [3 1 3; 1 2 2]) ... % Rate 2/3 code with K=3
);

% Define all puncturing patterns
puncture_patterns = struct(...
    'rate_8_9', [1 1 1 1 0 1 1 1 1 0 0 0 1 0 0 0], ... % Keeps 9/16
    'rate_4_5', [1 1 1 1 1 1 1 1 1 0 0 0 1 0 0 0], ... % Keeps 10/16
    'rate_2_3', [1 1 1 1 1 1 1 1 1 0 1 0 1 0 1 0]  ... % Keeps 12/16
);

% Process for each error probability
for p_idx = 1:length(error_probs)
    p = error_probs(p_idx);
    
    % Initialize output videos with p in filename
    output_no_coding = VideoWriter(sprintf('no_coding_p%.2f.mp4', p), 'MPEG-4');
    output_conv_1_2 = VideoWriter(sprintf('conv_1_2_p%.2f.mp4', p), 'MPEG-4');
    output_conv_1_3 = VideoWriter(sprintf('conv_1_3_p%.2f.mp4', p), 'MPEG-4');
    output_conv_1_4 = VideoWriter(sprintf('conv_1_4_p%.2f.mp4', p), 'MPEG-4');
    output_conv_2_3 = VideoWriter(sprintf('conv_2_3_p%.2f.mp4', p), 'MPEG-4');
    output_punc_8_9 = VideoWriter(sprintf('punc_8_9_p%.2f.mp4', p), 'MPEG-4');
    output_punc_4_5 = VideoWriter(sprintf('punc_4_5_p%.2f.mp4', p), 'MPEG-4');
    output_punc_2_3 = VideoWriter(sprintf('punc_2_3_p%.2f.mp4', p), 'MPEG-4');

    % Open all video writers
    open(output_no_coding);
    open(output_conv_1_2);
    open(output_conv_1_3);
    open(output_conv_1_4);
    open(output_conv_2_3);
    open(output_punc_8_9);
    open(output_punc_4_5);
    open(output_punc_2_3);

    for i = 1:frames
        % Extract current frame
        ArrayOfFrames(i).cdata = a(:,:,:,i);
        
        % Convert frame to binary
        R = ArrayOfFrames(i).cdata(:,:,1);
        G = ArrayOfFrames(i).cdata(:,:,2);
        B = ArrayOfFrames(i).cdata(:,:,3); 
        Rdouble = double(R);
        Gdouble = double(G);
        Bdouble = double(B);
        Rbin = de2bi(Rdouble, 8);  
        Gbin = de2bi(Gdouble, 8);
        Bbin = de2bi(Bdouble, 8);
        R_bin_reshaped = reshape(Rbin, [144, 176, 8]);
        G_bin_reshaped = reshape(Gbin, [144, 176, 8]);
        B_bin_reshaped = reshape(Bbin, [144, 176, 8]);
        R_stream = reshape(R_bin_reshaped, [], 1);
        G_stream = reshape(G_bin_reshaped, [], 1);
        B_stream = reshape(B_bin_reshaped, [], 1);
        binary_stream = [R_stream; G_stream; B_stream];
        
        % Process with all methods
        [no_coding_frame, conv_1_2_frame, conv_1_3_frame, conv_1_4_frame, ...
         conv_2_3_frame, punc_8_9_frame, punc_4_5_frame, punc_2_3_frame] = ...
            process_all_schemes(binary_stream, s, p, trellises, puncture_patterns);
       
        % Write to all output videos (using your exact frame reconstruction)
        writeVideo(output_no_coding, no_coding_frame);
        writeVideo(output_conv_1_2, conv_1_2_frame);
        writeVideo(output_conv_1_3, conv_1_3_frame);
        writeVideo(output_conv_1_4, conv_1_4_frame);
        writeVideo(output_conv_2_3, conv_2_3_frame);
        writeVideo(output_punc_8_9, punc_8_9_frame);
        writeVideo(output_punc_4_5, punc_4_5_frame);
        writeVideo(output_punc_2_3, punc_2_3_frame);
    end
    
    % Close all video writers
    close(output_no_coding);
    close(output_conv_1_2);
    close(output_conv_1_3);
    close(output_conv_1_4);
    close(output_conv_2_3);
    close(output_punc_8_9);
    close(output_punc_4_5);
    close(output_punc_2_3);
end
% Main processing function with your exact reconstruction logic
function [no_coding_frame, conv_1_2_frame, conv_1_3_frame, conv_1_4_frame, ...
          conv_2_3_frame, punc_8_9_frame, punc_4_5_frame, punc_2_3_frame] = ...
         process_all_schemes(binary_stream, s, p, trellises, patterns)
    
    % 1. No channel coding (direct BSC)
    noisy_bits = bsc(binary_stream, p);
    no_coding_frame = reconstruct_frame(noisy_bits, s);
    
    % 2. Conventional coding (all rates)
    % Rate 1/2
    encoded_1_2 = convenc(binary_stream, trellises.rate_1_2);
    noisy_1_2 = bsc(encoded_1_2, p);
    decoded_1_2 = vitdec(noisy_1_2, trellises.rate_1_2, 34, 'trunc', 'hard');
    conv_1_2_frame = reconstruct_frame(decoded_1_2(1:length(binary_stream)), s);
    
    % Rate 1/3
    encoded_1_3 = convenc(binary_stream, trellises.rate_1_3);
    noisy_1_3 = bsc(encoded_1_3, p);
    decoded_1_3 = vitdec(noisy_1_3, trellises.rate_1_3, 34, 'trunc', 'hard');
    conv_1_3_frame = reconstruct_frame(decoded_1_3(1:length(binary_stream)), s);
    
    % Rate 1/4
    encoded_1_4 = convenc(binary_stream, trellises.rate_1_4);
    noisy_1_4 = bsc(encoded_1_4, p);
    decoded_1_4 = vitdec(noisy_1_4, trellises.rate_1_4, 34, 'trunc', 'hard');
    conv_1_4_frame = reconstruct_frame(decoded_1_4(1:length(binary_stream)), s);
    
    % Rate 2/3 (special handling for multiple inputs)
    % Rate 2/3 encoding (fixed)
    if mod(length(binary_stream), 2) ~= 0
        binary_stream = [binary_stream; 0];  % Pad if odd
    end
    encoded_2_3 = convenc(binary_stream, trellises.rate_2_3);
    noisy_2_3 = bsc(encoded_2_3, p);
    decoded_2_3 = vitdec(noisy_2_3, trellises.rate_2_3, 34, 'trunc', 'hard');
    % Remove padding if added
    conv_2_3_frame = reconstruct_frame(decoded_2_3(1:length(binary_stream)), s);
    
    % 3. Punctured coding (all patterns on rate 1/2)
    padded = [encoded_1_2; zeros(16-mod(length(encoded_1_2),16),1)];
    
    % Punctured 8/9
    punctured_8_9 = padded(logical(repmat(patterns.rate_8_9, 1, ceil(length(encoded_1_2)/16))));
    noisy_punc_8_9 = bsc(punctured_8_9, p);
    decoded_punc_8_9 = vitdec(noisy_punc_8_9, trellises.rate_1_2, 34, 'trunc', 'hard', patterns.rate_8_9);
    punc_8_9_frame = reconstruct_frame(decoded_punc_8_9(1:length(binary_stream)), s);
    
    % Punctured 4/5
    punctured_4_5 = padded(logical(repmat(patterns.rate_4_5, 1, ceil(length(encoded_1_2)/16))));
    noisy_punc_4_5 = bsc(punctured_4_5, p);
    decoded_punc_4_5 = vitdec(noisy_punc_4_5, trellises.rate_1_2, 34, 'trunc', 'hard', patterns.rate_4_5);
    punc_4_5_frame = reconstruct_frame(decoded_punc_4_5(1:length(binary_stream)), s);
    
    % Punctured 2/3
    punctured_2_3 = padded(logical(repmat(patterns.rate_2_3, 1, ceil(length(encoded_1_2)/16))));
    noisy_punc_2_3 = bsc(punctured_2_3, p);
    decoded_punc_2_3 = vitdec(noisy_punc_2_3, trellises.rate_1_2, 34, 'trunc', 'hard', patterns.rate_2_3);
    punc_2_3_frame = reconstruct_frame(decoded_punc_2_3(1:length(binary_stream)), s);
end

%rame reconstruction after recieving the bits
function reconstructed_frame = reconstruct_frame(binary_stream, s)
    % Split back into color channels
    bits_per_channel = s(1)*s(2)*8;
    R_stream_received = binary_stream(1:bits_per_channel);
    G_stream_received = binary_stream(bits_per_channel+1:2*bits_per_channel);
    B_stream_received = binary_stream(2*bits_per_channel+1:3*bits_per_channel);
    
    % Reshape back to original dimensions
    R_bin_received = reshape(R_stream_received, [s(1), s(2), 8]);
    G_bin_received = reshape(G_stream_received, [s(1), s(2), 8]);
    B_bin_received = reshape(B_stream_received, [s(1), s(2), 8]);
    
    % Convert back to pixel values
    R_recovered = uint8(bi2de(reshape(R_bin_received, [], 8)));
    G_recovered = uint8(bi2de(reshape(G_bin_received, [], 8)));
    B_recovered = uint8(bi2de(reshape(B_bin_received, [], 8)));
    
    % Reconstruct frame
    reconstructed_frame = zeros(s(1), s(2), 3, 'uint8');
    reconstructed_frame(:,:,1) = reshape(R_recovered, [s(1), s(2)]);
    reconstructed_frame(:,:,2) = reshape(G_recovered, [s(1), s(2)]);
    reconstructed_frame(:,:,3) = reshape(B_recovered, [s(1), s(2)]);
end