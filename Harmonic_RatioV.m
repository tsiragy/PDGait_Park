function [Signal_peaks, Signal_Loc] = HarmonicRatio( x, fs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%x = signal of interest (typically trunk/whole-body COM acceleration)
% fs = capture frequency 
ff = 1; %the fundamental frequency 
S_length = length(x); %signal length (i.e. number of samples) 
FFT_aCOG = fft(x);
mag_FFT_aCOG = abs(FFT_aCOG);
avg_FFT = mean(mag_FFT_aCOG);
FFT_norm = mag_FFT_aCOG/length(mag_FFT_aCOG);
P1 = FFT_norm(1:S_length/2+1);
P1(2:end-1) = 2*P1(2:end-1);
fr = fs*(0:(S_length/2))/S_length;
%[pks,frq] = findpeaks(P1, fr, 'MinPeakDistance',1, 'MinPeakHeight',0.019);
%hold on 
%figure(1) 




   


[Signal_peaks_ff,ff_loc] = max(P1); %Find the fundamental frequency 
FF_loc = fr(ff_loc);    %location of the fundamental frequency
        
    
    
nf_t = (FF_loc/2); %determine the nyqulist frequency or the "first harmonic coefficient" 
nf = fr(fr>(nf_t - 0.3) & fr<(nf_t + 0.3)); %search for surrounding peaks 
 [~, y1_nf] = ismember(nf,fr);
    Signal_Loc_nf = P1(y1_nf);
    Signal_peaks_nf = max(Signal_Loc_nf);
 
 
 %h3 = fr(fr > (Signal_Loc_nf.*6)-0.2 & fr< (Signal_Loc_nf.*6)+0.2);
 h3 = fr(fr > (nf_t.*3)-0.2 & fr< (nf_t.*3)+0.2);
 [~, y1_h3] = ismember(h3,fr);
    Signal_values_h3 = P1(y1_h3);
    Signal_peaks_h3 = max(Signal_values_h3);
 
 h4 = fr(fr > (nf_t.*4)-0.2 & fr< (nf_t.*4)+0.2);
 [~, y1_h4] = ismember(h4,fr);
    Signal_values_h4 = P1(y1_h4);
    Signal_peaks_h4 = max(Signal_values_h4);
 
 h5 = fr(fr > (nf_t.*5)-0.2 & fr< (nf_t.*5)+0.2);
 [~, y1_h5] = ismember(h5,fr);
    Signal_values_h5 = P1(y1_h5);
    Signal_peaks_h5 = max(Signal_values_h5);
 
 h6 = fr(fr > (nf_t.*6)-0.2 & fr< (nf_t.*6)+0.2);
 [~, y1_h6] = ismember(h6,fr);
    Signal_values_h6 = P1(y1_h6);
    Signal_peaks_h6 = max(Signal_values_h6);
 
 h7 = fr(fr > (nf_t.*7)-0.2 & fr< (nf_t.*7)+0.2);
 [~, y1_h7] = ismember(h7,fr);
    Signal_values_h7 = P1(y1_h7);
    Signal_peaks_h7 = max(Signal_values_h7);
 
 h8 = fr(fr > (nf_t.*8)-0.2 & fr< (nf_t.*8)+0.2);
 [~, y1_h8] = ismember(h8,fr);
    Signal_values_h8 = P1(y1_h8);
    Signal_peaks_h8 = max(Signal_values_h8);
 
 h9 = fr(fr > (nf_t.*9)-0.2 & fr< (nf_t.*9)+0.2);
 [~, y1_h9] = ismember(h9,fr);
    Signal_values_h9 = P1(y1_h9);
    Signal_peaks_h9 = max(Signal_values_h9);
 
 h10 = fr(fr > (nf_t.*10)-0.2 & fr< (nf_t.*10)+0.2);
[~, y1_h10] = ismember(h10,fr);
    Signal_values_h10 = P1(y1_h10);
    Signal_peaks_h10 = max(Signal_values_h10);
    

 h11 = fr(fr > (nf_t.*11)-0.2 & fr< (nf_t.*11)+0.2);
[~, y1_h11] = ismember(h11,fr);
    Signal_values_h11 = P1(y1_h11);
    Signal_peaks_h11 = max(Signal_values_h11);
    
 h12 = fr(fr > (nf_t.*12)-0.2 & fr< (nf_t.*12)+0.2);
[~, y1_h12] = ismember(h12,fr);
    Signal_values_h12 = P1(y1_h12);
    Signal_peaks_h12 = max(Signal_values_h12);
    
 h13 = fr(fr > (nf_t.*13)-0.2 & fr< (nf_t.*13)+0.2);
[~, y1_h13] = ismember(h13,fr);
    Signal_values_h13 = P1(y1_h13);
    Signal_peaks_h13 = max(Signal_values_h13);
    
 h14 = fr(fr > (nf_t.*14)-0.2 & fr< (nf_t.*14)+0.2);
[~, y1_h14] = ismember(h14,fr);
    Signal_values_h14 = P1(y1_h14);
    Signal_peaks_h14 = max(Signal_values_h14);

    
 h15 = fr(fr > (nf_t.*15)-0.2 & fr< (nf_t.*15)+0.2);
[~, y1_h15] = ismember(h15,fr);
    Signal_values_h15 = P1(y1_h15);
    Signal_peaks_h15 = max(Signal_values_h15);
    
 h16 = fr(fr > (nf_t.*16)-0.2 & fr< (nf_t.*16)+0.2);
[~, y1_h16] = ismember(h16,fr);
    Signal_values_h16 = P1(y1_h16);
    Signal_peaks_h16 = max(Signal_values_h16);
    
 h17 = fr(fr > (nf_t.*17)-0.2 & fr< (nf_t.*17)+0.2);
[~, y1_h17] = ismember(h17,fr);
    Signal_values_h17 = P1(y1_h17);
    Signal_peaks_h17 = max(Signal_values_h17);
    
 h18 = fr(fr > (nf_t.*18)-0.2 & fr< (nf_t.*18)+0.2);
[~, y1_h18] = ismember(h18,fr);
    Signal_values_h18 = P1(y1_h18);
    Signal_peaks_h18 = max(Signal_values_h18);
    
 h19 = fr(fr > (nf_t.*19)-0.2 & fr< (nf_t.*19)+0.2);
[~, y1_h19] = ismember(h19,fr);
    Signal_values_h19 = P1(y1_h19);
    Signal_peaks_h19 = max(Signal_values_h19);
    
 h20 = fr(fr > (nf_t.*20)-0.2 & fr< (nf_t.*20)+0.2);
[~, y1_h20] = ismember(h20,fr);
    Signal_values_h20 = P1(y1_h20);
    Signal_peaks_h20 = max(Signal_values_h20);
    


Signal_peaks = [Signal_peaks_nf Signal_peaks_ff Signal_peaks_h3 Signal_peaks_h4 Signal_peaks_h5 Signal_peaks_h6 Signal_peaks_h7 Signal_peaks_h8 Signal_peaks_h9 Signal_peaks_h10 Signal_peaks_h11 Signal_peaks_h12 Signal_peaks_h13 Signal_peaks_h14 Signal_peaks_h15 Signal_peaks_h16 Signal_peaks_h17 Signal_peaks_h18 Signal_peaks_h19 Signal_peaks_h20];

[~, tt34] = ismember(Signal_peaks, P1);
Signal_Loc = fr(tt34);

hold on 
figure (3) 
plot(fr,P1)
hold on
scatter(Signal_Loc, Signal_peaks)
%set(gca, 'YScale', 'log')
hold off
title('HR V Single-Sided Amplitude Spectrum')
xlabel('f(Hz)')
ylabel('P1(f)')

%saveas(figure(3), 'HR Vert.png')


end

