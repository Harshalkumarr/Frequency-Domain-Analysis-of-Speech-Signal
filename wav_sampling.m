function[]=wav_sampling()
                       % Reading .wav file.Function wavread('file_name') returns the number of
                       % samples,sampling frequency,number of bits

[n,fs,nbits]=wavread('va.wav');
sp_samples=n;
sam_freq=fs
num_bits=nbits;

                      %soundsc(sp_samples,sam_freq);
                    
                      % Storing generated samples in file
                      
                      %[status,msginfo]=xlswrite('sampled_val.csv',sp_samples);


time=(1/sam_freq)*length(sp_samples)

                     % generation of row vector t of length of samples linearly spaced between
                        % and including 0 and time(duration of recorded audio)

figure(1)
                     %t=linspace(0,time,length(sp_samples));
                     %subplot(2,2,1)
                     
plot(sp_samples)
xlabel('time(sec)')
ylabel('amplitude')
title('noisy input speech sample')


                      % Ploting magnitude spectrum of noisy signal
%figure(2)
%sp_samples_mags=abs(fft(sp_samples));
%subplot(2,1,1)
%plot(sp_samples_mags)
%xlabel('Dft bins')
%ylabel('magnitude')
%title('magnitude response of noisy signal')


                      %Plot first half of Dft (normalized frequency)

%num_bins=length(sp_samples_mags);
%subplot(2,1,2)
%plot([0:1/(num_bins/2-1):1],sp_samples_mags(1:num_bins/2))
%xlabel('normalized frequency ')
%ylabel('magnitude')





                      % Filtering-butterworth bandpass filter
figure(2)
Wp=[10,10000]/(sam_freq/2);
Ws=[8,11000]/(sam_freq/2);
Rp=3;
Rs=11.4;
[n,Wn]=buttord(Wp,Ws,Rp,Rs)
[b,a]=butter(n,Wn,'bandpass')


[z,p,k]=butter(n,Wn,'bandpass')
subplot(3,2,4)
freqz(b,a)
title('Butterworth bandpass filter');

figure(3)
zplane(z,p)
title('Pole-Zero Plot')

%[b,a] = zp2tf(z,p,k)

                       % filtering of signal using the b (coefficients of numerator) and
                       % a (coefficients of denominator) obtained from butterworth design function



flt_sp_samples=filter(b,a,sp_samples);
                       %soundsc(flt_sp_samples,sam_freq)





                       %xlswrite('butt_flt_sampled_val.csv',flt_sp_samples);
                       %figure(4)
                       %create_signal_flow(b,a)

                  
                       % Plot the filtered signal

figure(4)
plot(flt_sp_samples,'r')
xlabel('time(sec)')
ylabel('amplitude')
title('Butterworth Bandpass Filter')




                      % Filtering-->Moving Average Filter %
figure(5)                    
for n=2:length(sp_samples)-1
   mov_sp_samples(n)=(sp_samples(n-1)+sp_samples(n)+sp_samples(n+1))/3;
end
                      %hist(mov_sp_samples)
plot(mov_sp_samples)
title('moving average filter')
                      %xlswrite('mov_filt_sampled_val.csv',mov_sp_samples)
                     
                      %soundsc(mov_sp_samples,sam_freq)

                    

                      % Framing-breaking signal into frames of 0.1 seconds %
new_mov_sp_samples=mov_sp_samples;                     
frame_duration=0.1;
% gives number of samples of that particular frame
frame_len=frame_duration*sam_freq;
N=length(new_mov_sp_samples);
num_frames=floor(N/frame_len);

non_silence=zeros(N,1);

 % what the following code does
 %{ frame1=mov_sp_samples(1:frame_len)
 %  frame2=mov_sp_samples(frame_len+1:frame_len*2)
 %  frame3=mov_sp_samples(frame_len*2+1:frame_len*3) and so on (number of
 %  frames)
 
 count=0;
 
 for k=1:num_frames
     
     % extracting frame of speech
     
     frame=new_mov_sp_samples((k-1)*frame_len+1:frame_len*k);
     
     %  identify the non_silent frames by finding frames with max amplitude
     %  more than 0.2212
   
     max_val=max(frame);
 
    if (max_val > 0.18)
      count=count+1 ;
       non_silence((count-1)*frame_len+1 : frame_len*count)=frame;
     end
     
    
 end
 
 
 
figure(6)
     
     plot(non_silence)
     title('Silence Removed From Filtered Signal')
     xlabel('Time')
     ylabel('Amplitude Values')
     soundsc(non_silence,sam_freq)
    
     % new wave file don't contain silence
     
    % wavwrite(non_silence,sam_freq,'new_va.wav')
     
     % feature extration using 'db4' wavelet-Signal decomposition                
     figure(7)
     
     new_samples=wavread('new_va.wav');
     S=size(new_samples);
     l=wmaxlev(S,'db4');
     [C,L]=wavedec(new_samples,l,'db4');
     subplot(311)
     plot(new_samples)
     title('Silence Removed Signal')
     xlabel('Freq')
     ylabel('Amplitude Values')
     subplot(312)
     plot(C)
     title('Wavelet decomposition structure,level l')
     xlabel('Coefficients for approximation at level 17 & details at levels 17 to 1')
     rec_new_samples=waverec(C,L,'db4');
     subplot(313)
     plot(rec_new_samples)
     title('Reconstructed signal from decomposition')
     xlabel('Freq')
     ylabel('Amplitude values')
     
     
     % power spectrum
      figure(8)
      cA17=appcoef(C,L,'db4');          % appr. coefficient vals
      
      pow_cA17=abs(cA17).^2;
      max_pow_cA17=max(pow_cA17);
      norm_pow_cA17=pow_cA17/max_pow_cA17;
      
      subplot(321)
      plot(norm_pow_cA17)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of approximatios')
      % detail coefficients at levels 1-17
     [cd1,cd2,cd3,cd4,cd5,cd6,cd7,cd8,cd9,cd10,cd11,cd12,cd13,cd14,cd15,cd16,cd17]=detcoef(C,L,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])
      pow_cd1=abs(cd1).^2;
      max_pow_cd1=max(pow_cd1);
      norm_pow_cd1=pow_cd1/max_pow_cd1;
      subplot(322)
      plot(norm_pow_cd1)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 1 ')
    
      pow_cd2=abs(cd2).^2;
      max_pow_cd2=max(pow_cd2);
      norm_pow_cd2=pow_cd2/max_pow_cd2;
      
      subplot(323)
      plot(norm_pow_cd2)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 2 ')
      pow_cd3=abs(cd3).^2;
      max_pow_cd3=max(pow_cd3);
      norm_pow_cd3=pow_cd3/max_pow_cd3;
      subplot(324)
      plot(norm_pow_cd3)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 3 ')
      pow_cd4=abs(cd4).^2;
      max_pow_cd4=max(pow_cd4);
      norm_pow_cd4=pow_cd4/max_pow_cd4;
      subplot(325)
      plot(norm_pow_cd4)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 4 ')
      pow_cd5=abs(cd5).^2;
      max_pow_cd5=max(pow_cd5);
      norm_pow_cd5=pow_cd5/max_pow_cd5;
      subplot(326)
      plot(norm_pow_cd5)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 5 ')
      pow_cd6=abs(cd6).^2;
      max_pow_cd6=max(pow_cd6);
      norm_pow_cd6=pow_cd6/max_pow_cd6;
      figure(9)
      subplot(321)
      plot(norm_pow_cd6)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 6 ')
     
      pow_cd7=abs(cd7).^2;
      max_pow_cd7=max(pow_cd7);
      norm_pow_cd7=pow_cd7/max_pow_cd7;
     subplot(322)
      plot(norm_pow_cd7)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 7 ')
      pow_cd8=abs(cd8).^2;
      max_pow_cd8=max(pow_cd8);
      norm_pow_cd8=pow_cd8/max_pow_cd8;
      subplot(323)
      plot(norm_pow_cd8)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 8 ')
      pow_cd9=abs(cd9).^2;
      max_pow_cd9=max(pow_cd9);
      norm_pow_cd9=pow_cd9/max_pow_cd9;
      subplot(324)
      plot(norm_pow_cd9)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 9')
      pow_cd10=abs(cd10).^2;
      max_pow_cd10=max(pow_cd10);
      norm_pow_cd10=pow_cd10/max_pow_cd10;
      subplot(325)
      plot(norm_pow_cd10)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 10 ')
      pow_cd11=abs(cd11).^2;
      max_pow_cd11=max(pow_cd11);
      norm_pow_cd11=pow_cd11/max_pow_cd11;
      subplot(326)
      plot(norm_pow_cd11)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 11 ')
      pow_cd12=abs(cd12).^2;
      max_pow_cd12=max(pow_cd12);
      norm_pow_cd12=pow_cd12/max_pow_cd12;
       figure(10)
      subplot(321)
      plot(norm_pow_cd12)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 12 ')
      pow_cd13=abs(cd13).^2;
      max_pow_cd13=max(pow_cd13);
      norm_pow_cd13=pow_cd13/max_pow_cd13;
   
      subplot(322)
      plot(norm_pow_cd13)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 13 ')
      pow_cd14=abs(cd14).^2;
      max_pow_cd14=max(pow_cd14);
      norm_pow_cd14=pow_cd14/max_pow_cd14;
      subplot(323)
      plot(norm_pow_cd14)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 14 ')
      pow_cd15=abs(cd15).^2;
      max_pow_cd15=max(pow_cd15);
      norm_pow_cd15=pow_cd15/max_pow_cd15;
      subplot(324)
      plot(norm_pow_cd15)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 15 ')
      pow_cd16=abs(cd16).^2;
      max_pow_cd16=max(pow_cd16);
      norm_pow_cd16=pow_cd16/max_pow_cd16;
      subplot(325)
      plot(norm_pow_cd16)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 16 ')
      pow_cd17=abs(cd17).^2;
      max_pow_cd17=max(pow_cd17);
      norm_pow_cd17=pow_cd17/max_pow_cd17;
      subplot(326)
      plot(norm_pow_cd17)
      xlabel('Freq')
      ylabel('power')
      title('Power spectrum for coeff. of details at level 17 ')
     
      % Spectral features-spectral centroid,spectral flux,spectral
      % spread,spectral skewness
      
      mag_cd1=abs(cd1); mag_cd2=abs(cd2);
      mag_cd3=abs(cd3); mag_cd4=abs(cd4);
      mag_cd5=abs(cd5); mag_cd6=abs(cd6);
      
      mag_cd7=abs(cd7); mag_cd8=abs(cd8);
      mag_cd9=abs(cd9); mag_cd10=abs(cd10);
      mag_cd11=abs(cd11); mag_cd12=abs(cd12);
      
      mag_cd13=abs(cd13); mag_cd14=abs(cd14);
      mag_cd15=abs(cd15); mag_cd16=abs(cd16);
      mag_cd17=abs(cd17);
      
      % spectral centroid
     
      spec_cntroid_cd1=sum(norm_pow_cd1.*cd1)/sum(norm_pow_cd1);
      spec_cntroid_cd2=sum(norm_pow_cd2.*cd2)/sum(norm_pow_cd2);
      spec_cntroid_cd3=sum(norm_pow_cd3.*cd3)/sum(norm_pow_cd3);
      spec_cntroid_cd4=sum(norm_pow_cd4.*cd4)/sum(norm_pow_cd4);
      spec_cntroid_cd5=sum(norm_pow_cd5.*cd5)/sum(norm_pow_cd5);
      spec_cntroid_cd6=sum(norm_pow_cd6.*cd6)/sum(norm_pow_cd6);
      spec_cntroid_cd7=sum(norm_pow_cd7.*cd7)/sum(norm_pow_cd7);
      spec_cntroid_cd8=sum(norm_pow_cd8.*cd8)/sum(norm_pow_cd8);
      spec_cntroid_cd9=sum(norm_pow_cd9.*cd9)/sum(norm_pow_cd9);
      spec_cntroid_cd10=sum(norm_pow_cd10.*cd10)/sum(norm_pow_cd10);
      spec_cntroid_cd11=sum(norm_pow_cd11.*cd11)/sum(norm_pow_cd11);
      spec_cntroid_cd12=sum(norm_pow_cd12.*cd12)/sum(norm_pow_cd12);
      spec_cntroid_cd13=sum(norm_pow_cd13.*cd13)/sum(norm_pow_cd13);
      spec_cntroid_cd14=sum(norm_pow_cd14.*cd14)/sum(norm_pow_cd14);
      spec_cntroid_cd15=sum(norm_pow_cd15.*cd15)/sum(norm_pow_cd15);
      spec_cntroid_cd16=sum(norm_pow_cd16.*cd16)/sum(norm_pow_cd16);
      spec_cntroid_cd17=sum(norm_pow_cd17.*cd17)/sum(norm_pow_cd17);
        
      spec_centroid=[ spec_cntroid_cd1,spec_cntroid_cd2, spec_cntroid_cd3, spec_cntroid_cd4, spec_cntroid_cd5, spec_cntroid_cd6, spec_cntroid_cd7, spec_cntroid_cd8, spec_cntroid_cd9, spec_cntroid_cd10, spec_cntroid_cd11, spec_cntroid_cd12, spec_cntroid_cd13, spec_cntroid_cd14, spec_cntroid_cd15, spec_cntroid_cd16, spec_cntroid_cd17]
      figure(11)
       plot(1:17,spec_centroid,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10)
     % plot(spec_centroid)
     % scatter(1:17,spec_centroid,'MarkerEdgeColor','b','MarkerFaceColor','c')
      xlabel('spectral band number');
      ylabel('spectral centroid value');
      title('spectral band vs spectral centroid ');
      
      over_spec_centroid=(spec_cntroid_cd1+spec_cntroid_cd2+spec_cntroid_cd3+spec_cntroid_cd4+ spec_cntroid_cd5+ spec_cntroid_cd6+spec_cntroid_cd7+spec_cntroid_cd8+spec_cntroid_cd9+spec_cntroid_cd10+spec_cntroid_cd11+spec_cntroid_cd12+ spec_cntroid_cd13 + spec_cntroid_cd14+spec_cntroid_cd15+ spec_cntroid_cd16+spec_cntroid_cd17)/17
     
      
      % spectral flux
      
       spec_flux1=sqrt(sum(abs((norm_pow_cd2-norm_pow_cd1(1:length(norm_pow_cd2))).^2)));
       spec_flux2=sqrt(sum(abs((norm_pow_cd3-norm_pow_cd2(1:length(norm_pow_cd3))).^2)));
       spec_flux3=sqrt(sum(abs((norm_pow_cd4-norm_pow_cd3(1:length(norm_pow_cd4))).^2)));
       spec_flux4=sqrt(sum(abs((norm_pow_cd5-norm_pow_cd4(1:length(norm_pow_cd5))).^2)));
       spec_flux5=sqrt(sum(abs((norm_pow_cd6-norm_pow_cd5(1:length(norm_pow_cd6))).^2)));
       spec_flux6=sqrt(sum(abs((norm_pow_cd7-norm_pow_cd6(1:length(norm_pow_cd7))).^2)));
       spec_flux7=sqrt(sum(abs((norm_pow_cd8-norm_pow_cd7(1:length(norm_pow_cd8))).^2)));
       spec_flux8=sqrt(sum(abs((norm_pow_cd9-norm_pow_cd8(1:length(norm_pow_cd9))).^2)));
       spec_flux9=sqrt(sum(abs((norm_pow_cd10-norm_pow_cd9(1:length(norm_pow_cd10))).^2)));
       spec_flux10=sqrt(sum(abs((norm_pow_cd11-norm_pow_cd10(1:length(norm_pow_cd11))).^2)));
       spec_flux11=sqrt(sum(abs((norm_pow_cd12-norm_pow_cd11(1:length(norm_pow_cd12))).^2)));
       spec_flux12=sqrt(sum(abs((norm_pow_cd13-norm_pow_cd12(1:length(norm_pow_cd13))).^2))); 
       spec_flux13=sqrt(sum(abs((norm_pow_cd14-norm_pow_cd13(1:length(norm_pow_cd14))).^2)));
       spec_flux14=sqrt(sum(abs((norm_pow_cd15-norm_pow_cd14(1:length(norm_pow_cd15))).^2)));
       spec_flux15=sqrt(sum(abs((norm_pow_cd16-norm_pow_cd15(1:length(norm_pow_cd16))).^2)));
       spec_flux16=sqrt(sum(abs((norm_pow_cd17-norm_pow_cd16(1:length(norm_pow_cd17))).^2)));
     
       spec_flux=[spec_flux1,spec_flux2,spec_flux3,spec_flux4,spec_flux5,spec_flux6,spec_flux7,spec_flux8,spec_flux9,spec_flux10,spec_flux11,spec_flux12,spec_flux13,spec_flux14,spec_flux15,spec_flux16]
       figure(12)
        plot(1:16,spec_flux,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10)
      % plot(spec_flux)
     % scatter(1:16,spec_flux,'MarkerEdgeColor','b','MarkerFaceColor','c')
       xlabel('spectra band number');
       ylabel('spectral flux value');
       title('spectral band vs spectral flux');       
       
       over_spec_flux=(spec_flux1+spec_flux2+spec_flux3+spec_flux4+spec_flux5+spec_flux6+spec_flux7+spec_flux8+spec_flux9+spec_flux10+spec_flux11+spec_flux12+spec_flux13+spec_flux14+spec_flux15+spec_flux16)/16
     
       
       % spectral spread
       
       spec_spread1=sqrt(sum((norm_pow_cd1-spec_cntroid_cd1).^2))/sqrt(sum((norm_pow_cd1).^2));
       spec_spread2=sqrt(sum((norm_pow_cd2-spec_cntroid_cd2).^2))/sqrt(sum((norm_pow_cd2).^2));
       spec_spread3=sqrt(sum((norm_pow_cd3-spec_cntroid_cd3).^2))/sqrt(sum((norm_pow_cd3).^2));
       spec_spread4=sqrt(sum((norm_pow_cd4-spec_cntroid_cd4).^2))/sqrt(sum((norm_pow_cd4).^2));
       spec_spread5=sqrt(sum((norm_pow_cd5-spec_cntroid_cd5).^2))/sqrt(sum((norm_pow_cd5).^2));
       spec_spread6=sqrt(sum((norm_pow_cd6-spec_cntroid_cd6).^2))/sqrt(sum((norm_pow_cd6).^2));
       spec_spread7=sqrt(sum((norm_pow_cd7-spec_cntroid_cd7).^2))/sqrt(sum((norm_pow_cd7).^2));
       spec_spread8=sqrt(sum((norm_pow_cd8-spec_cntroid_cd8).^2))/sqrt(sum((norm_pow_cd8).^2));
       spec_spread9=sqrt(sum((norm_pow_cd9-spec_cntroid_cd9).^2))/sqrt(sum((norm_pow_cd9).^2));       
       spec_spread10=sqrt(sum((norm_pow_cd10-spec_cntroid_cd10).^2))/sqrt(sum((norm_pow_cd10).^2));
       spec_spread11=sqrt(sum((norm_pow_cd11-spec_cntroid_cd11).^2))/sqrt(sum((norm_pow_cd11).^2));
       spec_spread12=sqrt(sum((norm_pow_cd12-spec_cntroid_cd12).^2))/sqrt(sum((norm_pow_cd12).^2));
       spec_spread13=sqrt(sum((norm_pow_cd13-spec_cntroid_cd13).^2))/sqrt(sum((norm_pow_cd13).^2));
       spec_spread14=sqrt(sum((norm_pow_cd14-spec_cntroid_cd14).^2))/sqrt(sum((norm_pow_cd14).^2));
       spec_spread15=sqrt(sum((norm_pow_cd15-spec_cntroid_cd15).^2))/sqrt(sum((norm_pow_cd15).^2));
       spec_spread16=sqrt(sum((norm_pow_cd16-spec_cntroid_cd16).^2))/sqrt(sum((norm_pow_cd16).^2));
       spec_spread17=sqrt(sum((norm_pow_cd17-spec_cntroid_cd17).^2))/sqrt(sum((norm_pow_cd17).^2));
       
       spec_spread=[ spec_spread1, spec_spread2, spec_spread3, spec_spread4, spec_spread5, spec_spread6, spec_spread7, spec_spread8, spec_spread9, spec_spread10, spec_spread11, spec_spread12, spec_spread13, spec_spread14, spec_spread15, spec_spread16, spec_spread17]
       figure(13);
       plot(1:17,spec_spread,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10)
       %plot(spec_spread)
       %scatter(1:17,spec_spread,'MarkerEdgeColor','b','MarkerFaceColor','c')
       xlabel('spectral band number');
       ylabel('spectral spread value');
       title('spectral band vs spectral spread');
       
       over_spec_spread=sum(spec_spread)/17
       
       % spectral skewness
       spec_skew1=(sum(((cd1-spec_cntroid_cd1).^2).*mag_cd1)/(sum(mag_cd1)));
       spec_skew2=(sum(((cd2-spec_cntroid_cd2).^2).*mag_cd2)/(sum(mag_cd2)));
       spec_skew3=(sum(((cd3-spec_cntroid_cd3).^2).*mag_cd3)/(sum(mag_cd3)));
       spec_skew4=(sum(((cd4-spec_cntroid_cd4).^2).*mag_cd4)/(sum(mag_cd4)));
       spec_skew5=(sum(((cd5-spec_cntroid_cd5).^2).*mag_cd5)/(sum(mag_cd5)));
       spec_skew6=(sum(((cd6-spec_cntroid_cd6).^2).*mag_cd6)/(sum(mag_cd6)));
       spec_skew7=(sum(((cd7-spec_cntroid_cd7).^2).*mag_cd7)/(sum(mag_cd7)));
       spec_skew8=(sum(((cd8-spec_cntroid_cd8).^2).*mag_cd8)/(sum(mag_cd8)));
       spec_skew9=(sum(((cd9-spec_cntroid_cd9).^2).*mag_cd9)/(sum(mag_cd9)));
       spec_skew10=(sum(((cd10-spec_cntroid_cd10).^2).*mag_cd10)/(sum(mag_cd10)));
       spec_skew11=(sum(((cd11-spec_cntroid_cd11).^2).*mag_cd11)/(sum(mag_cd11)));
       spec_skew12=(sum(((cd12-spec_cntroid_cd12).^2).*mag_cd12)/(sum(mag_cd12)));
       spec_skew13=(sum(((cd13-spec_cntroid_cd13).^2).*mag_cd13)/(sum(mag_cd13)));
       spec_skew14=(sum(((cd14-spec_cntroid_cd14).^2).*mag_cd14)/(sum(mag_cd14)));
       spec_skew15=(sum(((cd15-spec_cntroid_cd15).^2).*mag_cd15)/(sum(mag_cd15)));
       spec_skew16=(sum(((cd16-spec_cntroid_cd16).^2).*mag_cd16)/(sum(mag_cd16)));
       spec_skew17=(sum(((cd17-spec_cntroid_cd17).^2).*mag_cd17)/(sum(mag_cd17)));
       
       spec_skew=[spec_skew1,spec_skew2,spec_skew3,spec_skew4,spec_skew5,spec_skew6,spec_skew7,spec_skew8,spec_skew9,spec_skew10,spec_skew11,spec_skew12,spec_skew13,spec_skew14,spec_skew15,spec_skew16,spec_skew17]
       figure(14)
        plot(1:17,spec_skew,'--rs','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10)
       %plot(spec_skew)
       %scatter(1:17,spec_skew,'MarkerEdgeColor','b','MarkerFaceColor','c');
       xlabel('spectral band number');
       ylabel('spectral skew value');
       title('spectral band vs spectral skew');
       over_spec_skew=sum(spec_skew)/17
       
       % feature vector
       
       feature_vector=[over_spec_centroid,over_spec_flux,over_spec_spread,over_spec_skew]
       
       
       
       
end