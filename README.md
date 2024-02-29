# BER-Performance-Analysis-of-Hamming-Codes
BER Performance Analysis of Hamming Codes for Quadrature Phase Shift Keying modulation in AWGN channel. 
As a project in Communication System 2 lecture, me and my 2 friend developed QPSK modulation functions and Hamming Encoder and Decoder functions without using MATLAB's built-in functions. The steps of process is:
<br>1- Generating transmitted signal (randomly generated, consist of 0's and 1's.
<br>2- Hamming Encoding the signals. 7,4 Hamming Code is used.
<br>3- QPSK Modulation the Hamming Coded Signal.
<br>4- For simulating AWGN channel, AWGN noise is added.
<br>5- QPSK Demodulation.
<br>6- Hamming decoding and error correction of the received signal.
<br>7- Error calculation with&without hamming code.
<br> To achieve more accurate results, there is an SNR array generated. This array includes SNR values to be tested. All these steps are performed numberoftrial times for each SNR value. After that, the graph is generated to show the BER result for Hamming coded and uncoded signals.

