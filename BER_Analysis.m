% Hamming Code Project COM2
% İlker KESER, Ufuk Kaan BAYRAM, Yahya EKİN
%
%
%
% Main Loop ------------------------------------------------------------

% Step 0- Trials and finding the mean with multiple trials and BER
% calculation ----------------------------------------------------------
numberoftrials = 200; 
Error_for_trial_with_hammng = zeros(1, numberoftrials);  %to store  error per trial
Error_for_trial_without_hamming = zeros(1, numberoftrials); %to store  error per trial
snrArray = -5:1:10; 
errorArray_with_hamming = []; %to store average errors
errorArray_without_hamming = []; %to store average errors

for k=1:1:length(snrArray)
    snr = snrArray(k);
    X = ['SNR: ', num2str(snr)];
    disp(X);
for trial = 1:numberoftrials
disp(['Trial ' num2str(trial) ':']);



% Step 1- SIGNAL GENERATION --------------------------------------------

datagenerated = randi([0, 1], 1, 1000); %generating the transmitted signal

disp("Generated Data");
disp(datagenerated);
datalength = 4; %number of data bits for hamming code.
numberofarrays = length(datagenerated) / datalength; %number of arrays to be hamming coded
transmitted_data_subarray = cell(1, numberofarrays); %createing sub arrays for hamming code. Arrays of 4 bits generated.

for i = 1:numberofarrays %generating the subarrays
    transmitted_data_subarray{i} = datagenerated((i - 1) * datalength + 1 : i * datalength);
end

% Step 2- HAMMING ENCODING THE SUBARRAYS ---------------------------------

hamming_encoded_data = zeros(1, 7); %generating the array to hold hamming coded data
ParityArray = zeros(1,3);
for i = 1:numberofarrays
    [hamming_encoded_data(i, :), ParityArray(i, :)] = hamming_encoding(transmitted_data_subarray{i}); %hamming_encoding returns the hamming encoded data (7bit)
end

hamming_result_single_row = reshape(hamming_encoded_data.', 1, []); 
disp("Hamming Encoded Data:");
disp(hamming_result_single_row);

% Step 3- QPSK MODULATION ------------------------------------------------

symbols_withhamming = qpsk_modulation(hamming_result_single_row);
symbols_without_hamming = qpsk_modulation(datagenerated);
% disp("QPSK Modulated Signal"); % disp(symbols); % IQ diyagramını çiz
% figure;
% plot(symbols, '.');
% title('QPSK Modulation Constellation before transmission');
% xlabel('I (In-Phase)');
% ylabel('Q (Quadrature)');
% grid on;

% Step 4- SNR ADDING ----------------------------------------------------

%snr = -2; %// in dB
received_data_with_hamming = awgn(symbols_withhamming, snr, 'measured');
received_data_without_hamming = awgn(symbols_without_hamming, snr, 'measured');
% fh2 = figure;
% plot(real(received_data), imag(received_data), '.');
% title(['QPSK constellation at an SNR of ' num2str(snr) ' dB']);
% xlabel('real part');
% ylabel('imaginary part');
% grid on;

% Step 5- QPSK Demodulation ---------------------------------------------

received_bits_with_hamming = qpsk_demodulation(received_data_with_hamming); 
received_bits_without_hamming = qpsk_demodulation(received_data_without_hamming);

%Step 6- Hamming Decoding and Error Correction --------------------------

for i = 1:numberofarrays  
    
    result = hammingdem(received_bits_with_hamming((i - 1) * 7 + 1 : i * 7));
    encoded_data{i} = result.databits; 
    encoded_parity{i} = result.parity_bits;
    
end

dataresult = reshape(encoded_data.', 1, []);
dataresult2= cell2mat(dataresult);
disp("Hamming decoded Data Bits in a Single Row before correction")
disp(dataresult2)

% The steps of error correction
% 1-Ontain encoded data and encoded_parity
% 2- Generate new parity based on encoded_data
% 3- Compare the new_parity and encoded_parity and generate another array
% with this logic. If the bits of encoded_parity and new_parity is
% different, place 1. If they are same, place 0. Then convert this binary
% data to decimal. The deicmal value show the error index of received_data.
% Complement the value of received data in given error index. 

% We use hamming_encoding function to generate new_parity based on
% encoded_data.
new_hamming_encoded_data = zeros(1,7);
new_parity = zeros(1,3);
for i = 1:numberofarrays
    [new_hamming_encoded_data(i, :), new_parity(i,:)] = hamming_encoding(encoded_data{i}); %hamming_encoding returns the hamming encoded data (7bit)
end
new_parity_cell = mat2cell(new_parity, ones(1, numberofarrays), 3);
new_parity_cell = reshape(new_parity_cell.', 1, []); 
parity_error_decimal = zeros(1,numberofarrays);

for i = 1:numberofarrays
    parity_error_decimal(i) = bi2de(xor(new_parity_cell{i},encoded_parity{i}));
end

firstData = encoded_data{1};
first_element = firstData(2);
rowVector = zeros(1,4);
corrected_data_withHamming = repmat({rowVector}, 1, numberofarrays);


for i=1:numberofarrays
    tempData = encoded_data{i};
    
    switch parity_error_decimal(i)
    case 3
        tempData(1) = ~tempData(1);
    case 5
        tempData(2) = ~tempData(2);
    case 6
        tempData(3) = ~tempData(3);
    case 7
        tempData(4) = ~tempData(4);
    otherwise
        tempData = tempData;
    end
    corrected_data_withHamming{i} = tempData;
end

% disp('Corrected Data With Hamming');
% disp(corrected_data_withHamming);
corrected_data_withHamming_single_row1 = reshape(corrected_data_withHamming.', 1, []);
corrected_data_withHamming_single_row= cell2mat(corrected_data_withHamming_single_row1);
disp("Corrected Hamming Data Bits in a Signle Row")
disp(corrected_data_withHamming_single_row)

% Step 7- Error Calculation with Hamming & without hamming-----------------
e_with_hamming=0;
for i=1:length(datagenerated)
    if datagenerated(i) ~= corrected_data_withHamming_single_row(i)
        e_with_hamming=e_with_hamming+1;
    end
end

% disp("BER with Hamming:"); disp(e_with_hamming);
disp("With Ham BER Normalized Error:"); p1 = (e_with_hamming/length(datagenerated)); disp(p1);

e_without_hamming=0;
for i=1:length(datagenerated)
    if datagenerated(i) ~= received_bits_without_hamming(i)
        e_without_hamming=e_with_hamming+1;
    end
end

% disp("BER without Hamming:"); disp(e_without_hamming);
disp("Without Ham BER Normalized Error:"); p2 = (e_without_hamming/length(datagenerated)); disp(p2);

% Store error percentage for this trial
Error_for_trial_with_hamming(trial) = p1;
Error_for_trial_without_hamming(trial) = p2;

end

% Calculate and display average error over all trials
disp('Average Error with Hamming:');
disp(mean(Error_for_trial_with_hamming));
errorArray_with_hamming(k) = mean(Error_for_trial_with_hamming);

disp('Average Error without Hamming:');
disp(mean(Error_for_trial_without_hamming));
errorArray_without_hamming(k) = mean(Error_for_trial_without_hamming);

end

figure(1)
semilogy(snrArray, errorArray_with_hamming, '-o', 'LineWidth', 1.0);
xlabel('SNR Values');
ylabel('Average Error');
xticks(snrArray);
hold on
semilogy(snrArray, errorArray_without_hamming, '-o', 'LineWidth', 1.0);
title('BER with different SNR');
legend('BER with Hamming','BER without Hamming');
grid on;


% Functions -----------------------------------------------------------

% HAMMING ENCODER 
function [HamArray, ParityArray] = hamming_encoding(D)
    P1_Check = [D(1) D(2) D(4)];
    if sum(P1_Check) == 1 || sum(P1_Check) == 3
        p1 = 1;
    else
        p1 = 0;
    end

    P2_Check = [D(1) D(3) D(4)];
    if sum(P2_Check) == 1 || sum(P2_Check) == 3
        p2 = 1;
    else
        p2 = 0;
    end

    P3_Check = [D(2) D(3) D(4)];
    if sum(P3_Check) == 1 || sum(P3_Check) == 3
        p3 = 1;
    else
        p3 = 0;
    end

    HamArray = [p1 p2 D(1) p3 D(2) D(3) D(4)];
    ParityArray = [p1 p2 p3];
end

% HAMMING DECODER
function result = hammingdem(decod)
paritylength= ceil(log2(length(decod)+1));

parity=zeros(1,paritylength);
data=zeros(1,length(decod)-paritylength);

p=1;
d=1;

for i=1:length(decod)

   if log2(i) == fix(log2(i))
        
        parity(p) = decod(i);
        p=p+1;
    else
      
        data(d) = decod(i);
        d=d+1;
    end
end
  result.databits=data;
  result.parity_bits=parity;
end
% QPSK MODULATOR
function symbols = qpsk_modulation(data_bits)
    
    data_pairs = reshape(data_bits, 2, length(data_bits)/2)';
   
    qpsk_symbols = bi2de(data_pairs, 'left-msb');
    qpsk_corrected = [];
    for i=1:length(qpsk_symbols)
        if (qpsk_symbols(i) == 3)
            qpsk_corrected(i) = 1;
        end
        if (qpsk_symbols(i) == 1)
            qpsk_corrected(i) = 3;
        end
        if (qpsk_symbols(i) == 0)
            qpsk_corrected(i) = 5;
        end
        if (qpsk_symbols(i) == 2)
            qpsk_corrected(i) = 7;
        end
    end
    modulation_table = exp(1j * pi/4 * qpsk_corrected);
    symbols = modulation_table(:).';
end
%% QPSK DEMODULATION FUNCTION
function received_bits = qpsk_demodulation(received_data)

%This function has 3 steps. 
% 1:Thresholding and assigning quadrants.
% 2:Assigning desired complex values to quadrants.
% 3:Assigning bits to related quadrants (Converting symbols to bits).

quadrant = zeros(size(received_data)); %generating quadrant array

% Thresholding
for i=1:1:length(received_data) %Gelen gürültülü sinyalin thresholding ile bölgelere ayrılması
    if (real(received_data(i))>0 && imag(received_data(i))>0)
        quadrant(i) = 1; %her bir sembol, eğer tanımlı quadrant içindeyse o quadrantın numarasını alıyoruz.
    end
    if (real(received_data(i))<0 && imag(received_data(i))>0)
        quadrant(i) = 2;
    end
    if (real(received_data(i))<0 && imag(received_data(i))<0)
        quadrant(i) = 3;
    end
    if (real(received_data(i))>0 && imag(received_data(i))<0)
        quadrant(i) = 4;
    end
end
% Assigning to quadrants
thresholded_sig = zeros(size(received_data));
for i=1:1:length(quadrant) %received sinyalin quadrantlarına göre olması gereken bölgeye atanması
    if (quadrant(i)==1)
        thresholded_sig(i) = 0.707106781186548 + 0.707106781186548i;
    end
    if (quadrant(i)==2)
        thresholded_sig(i) = -0.707106781186548 + 0.707106781186548i;
    end
    if (quadrant(i)==3)
        thresholded_sig(i) = -0.707106781186548 - 0.707106781186548i;
    end
    if (quadrant(i)==4)
        thresholded_sig(i) = 0.707106781186548 - 0.707106781186548i;
    end

end

% figure(5) %Threshold uygulandıktan sonra oluşan sinyalin constellation u
% plot(real(thresholded_sig), imag(thresholded_sig), '.');
% title(['QPSK constellation at an SNR of ' num2str(snr) ' dB, Signal after threshodling']);
% xlabel('real part');
% ylabel('imaginary part');
% grid on;

% quadrantlara göre bit atama.
% Get the number of symbols
num_symbols = length(quadrant);

% Initialize the received_bits array
received_bits=[];
% Loop through each symbol
   
for i=1:num_symbols
        
        if (quadrant(i)==1)
            received_bits = [received_bits , [1 1]];
        end
        if (quadrant(i)==2)
          received_bits = [received_bits , [0 1]];
        end
        if (quadrant(i)==3)
            received_bits = [received_bits , [0 0]];
        end
        if (quadrant(i)==4)
           received_bits = [received_bits , [1 0]];
        end
end

end
