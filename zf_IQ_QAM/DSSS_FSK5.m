% ****************************************************
%     Performance Simulation of DSSS-FSK
%        Error Probability of DSSS-FSK 
%              and Rake Receiving
%     Demodulation of spreaded signal is first
%
% ****************************************************
% clear all;
% close all;

NM=1000;
NC=31;
ncyc=1;                                       % Statistic number

tc1=0.005;                                    % chip宽度
td=NC*tc1;                                    % FSK的码元宽度

fs1=10e3;                                     % sample frequency
fc1=1e3;                                      % tone frequency  +1
fc2=2e3;                                      % tone frequency  -1
bs1=fs1*tc1;                                
bs2=fs1*td;
eps1=0.38;

flag_ds=1;
flag_rake=1;

err_fsk=0;err_fsk_en=0;
err_ds=0;err_rake=0;
err_ds_en=0;err_rake_en=0;

% ****** 载波信号 ******
tt=0:1/fs1:(td-1/fs1);
fsk1_cos=cos(2*pi*fc1*tt);
fsk2_cos=cos(2*pi*fc2*tt);

tt1=0:1/fs1:(tc1-1/fs1);
ds1_cos=cos(2*pi*fc1*tt1);
ds2_cos=cos(2*pi*fc2*tt1);

% ****** PN码产生 ******
a=[1 1 1 1 1];
for i=1:NC
	pn(i)=a(5);
	a(5)=a(4);  
    a(4)=a(3);
    a(3)=a(2);
    a(2)=a(1);
	 if((a(1)+a(2)+a(4)+pn(i))==1 || (a(1)+a(2)+a(4)+pn(i))==3)
		a(1)=1;
     else
        a(1)=0;
     end
end
for i=1:NC
    if (pn(i)==0)
        pn(i)=-1;
    end
end

% ****** Coherant Sequence ******
seqq_ds1=0;seqq_rake=0;
seqq_f1=0;seqq_f2=0;
ze=zeros(1,length(ds1_cos));
for j=1:NC                                       % DSSS
    if pn(j)==1
        seqq_f1=[seqq_f1 ds1_cos];
        seqq_f2=[seqq_f2 ds2_cos];
    else
        seqq_f1=[seqq_f1 ze];   
        seqq_f2=[seqq_f2 ze];  
    end
end

seqq_f1=seqq_f1(2:length(seqq_f1));
seqq_f2=seqq_f2(2:length(seqq_f2));
Lpn_ds=length(seqq_f1);

% ****** 信道参数 ******

pl=0;tp=0;
rpt =[134.227   65.569
      134.879   66.757
      137.526   66.396
      142.250   65.923];
  
% rpt =[134.2152   65.867
%   140.5465   66.061
%   135.7882   65.906
%   134.5818   65.971
%   135.0426   65.951
%   136.2212   66.079
%   141.2072   66.315];

% rpt = [121.2657   46.1790
%   121.1627   46.3129
%   126.2940   46.1332
%   128.4948   46.2669
%   143.3866   46.1812
%   128.2482   46.4922
%   132.3180   46.4071
%   128.5326   46.2748
%   129.2960   46.2072
%   128.8532   46.3496
%   131.8432   46.3626
%   128.9838   46.3262
%   133.0679   46.3668
%   141.5868   46.2857];

% rpt =[134.227   65.569];

pl=rpt(:,1);
tp=rpt(:,2);
pl=min(pl)-pl;
for j=1:length(pl)
    pl(j)=10^(pl(j)/20);
end
tp1=(tp-min(tp))*fs1;
tp1=round(tp1);
tpm=max(tp1);

% ****** 滤波器 ******

N2=256;
Wn1=2*[950,1050]/fs1;
be1=fir1(N2,Wn1 ,'bandpass');                    %滤波550－1050

Wn2=2*[1950,2050]/fs1;
be2=fir1(N2,Wn2 ,'bandpass');                    %滤波1950－2050

be3 =fir1(256,40/fs1);                           %包络低通

% % ****** 测试滤波器延时 ******
t_LFM=1.5;
lzero=2*fs1*td;
t=0:1/fs1:(t_LFM-1/fs1);   
x=chirp(t,1500,t_LFM,1900);
K=length(x);

% ee=xcorr(qq1,x);
% ee=ee/max(ee);
% [nmax,num1]=max(ee);
% num1=num1-length(qq1)+K+lzero-32;
% qq2=qq1(num1+1:length(qq1));                       
% 
% dd1800= filter(be2,1,qq2);                       % 滤波1780－1820
% dd_1= dd1800.^2;
% ddd_1= filter(be3,1,dd_1);                       % 包络检波
% 
% dd1600= filter(be21,1,qq2);                      % 滤波1580－1620
% dd_0= dd1600.^2;
% ddd_0= filter(be3,1,dd_0);                       % 包络检波
% 
% ddss=ddd_1-ddd_0;                                % 信号包络
% Lm=length(fsk1_cos);                             % 码元长度
% % tqq2=1:length(ddss);
% % figure(2),plot(tqq2,ddss);
% 
% for i=1:Lm
%     qq3=ddss(i:Lm+i-1);
%     qq4(i)=sum(qq3);
% end
% 
% [Y,NN]=max(qq4);
% nint_fsk=fix(bs2/nsam_fsk);
% nint_ds=fix(bs1/nsam_ds);                        % sample interval of DSSS signal
% nint_rake=fix(bs1/nsam_rake);
% MM_ds=fix(NN/nint_ds);
% MM_rake=fix(NN/nint_rake);
% 
% ****** 误码分析 ******
msg=0;smag=0;
SNR=-10:5:10;
% SNR=20;
en_fsk=0;en_fsk_en=0;
en_ds=0;en_ds_en=0;
en_rake=0;en_rake_en=0;

MM0=0.5;
MM=round(MM0*Lpn_ds);
MM1=10;
for jj=1:ncyc
    for ii=1:length(SNR)
        snr=SNR(ii);
        smag=0;
    
        msg=randi([0,1],1,NM);                      % 随机信息序列
        smag=msg;
        nt=length(smag);

   % ****** FSK调制 ******
        if flag_fsk==1                       
            probt=0;
            for j=1:nt                                      
                if (smag(j)==1)
                    probt=[probt fsk1_cos];
                else
                    probt=[probt fsk2_cos];
                end
            end
            probt=probt(2:length(probt));
            lst=length(probt);
            
            sr=0;srr=0;
            for j=1:length(pl)
                srr=pl(j)*[zeros(1,tp1(j)) probt zeros(1,tpm-tp1(j))]+srr;
            end
            sr=awgn(srr,snr,'measured'); 
            
            ssr1=0;
            ssr1=sr(1:lst);

            d_1=0;dd_1=0;ddd_1=0;
            d_0=0;dd_0=0;ddd_0=0;
            ddss=0;ddss1=0;
            
            d1000= filter(be1,1,ssr1);                     % 滤波1000－1050
            d_1= d1000.^2;
            dd_1= filter(be3,1,d_1);                                    % 包络检波

            d2000= filter(be2,1,ssr1);                                 % 滤波1580－1620
            d_0= d2000.^2; 
            dd_0= filter(be3,1,d_0);                                    % 包络检波
            
            rec=0;
            for j=1:nt
                mayuan1=0;jf1=0;
                mayuan2=0;jf2=0;
                mayuan1=dd_1((j-1)*bs2+1:j*bs2);
                jf1 = mean(mayuan1);
                mayuan2=dd_0((j-1)*bs2+1:j*bs2);
                jf2 = mean(mayuan2);
                if(jf1>jf2)
                    rec(j)=1;
                else
                    rec(j)=0;
                end
            end
        end   
        
        if flag_ds==1 | flag_rake==1
            probt=0;st=0;
            for i=1:nt
                if (smag(i)==1)
                    for j=1:NC                                       % DSSS
                        if pn(j)==1
                            probt=[probt ds1_cos];
                        else
                            probt=[probt ds2_cos];
                        end
                    end
                 else
                    for j=1:NC                                      
                        if pn(j)==-1
                            probt=[probt ds1_cos];
                        else
                            probt=[probt ds2_cos];
                        end
                    end
                end
            end
            probt=probt(2:length(probt));
            st=[x zeros(1,lzero) probt];                           % 加同步头信号
            lst=length(probt);
        
        % ****** 信道与多径 ************
            srr=0;sr=0;
            for j=1:length(pl)
                srr=pl(j)*[zeros(1,tp1(j)) st zeros(1,tpm-tp1(j))]+srr;
            end
            sr=awgn(srr,snr,'measured');                                                  % White noise 

        % ****** 同步,滤波器及信号解调  ******
            ee=0;
            ee=xcorr(sr,x);
            ee=ee/max(ee);
            
      % ****** DSSS Decode ******
            if flag_ds==1
               dd=0;
               [nmax,num]=max(ee);
               num=num-length(sr)+K+lzero-32;
               dd=sr(num+1:length(sr));                                  % 初同步
               
               d_1=0;d_2=0;
               d_f1=0;d_f2=0;d_1=0;d_2=0;

               d_f1= filter(be1,1,dd);                                  % 滤波950－1050
               d_1= [d_f1 zeros(1,MM1*MM)];
               
               d_f2= filter(be2,1,dd);                                 % 滤波1050－2050
               d_2= [d_f2 zeros(1,MM1*MM)];
               
               rec_ds=0;jff1=0;jff2=0;
               for j=1:nt
                   jf1=0;jf2=0;
                   for i=1:MM1                                            % 码元细同步
                        mayuan1=0;jiekuo1=0;dd_1=0;ddd_1=0;
                        mayuan2=0;jiekuo2=0;dd_2=0;ddd_2=0;
                        mayuan1=d_1((j-1)*Lpn_ds+(i-1)*MM+1:j*Lpn_ds+(i-1)*MM);
                        mayuan2=d_2((j-1)*Lpn_ds+(i-1)*MM+1:j*Lpn_ds+(i-1)*MM);
                        jiekuo1=mayuan1.*seqq_f1;
                        jiekuo2=mayuan2.*seqq_f2;
                        
                        for k=1:NC
                            jj1=00;jj2=0;
                            jj1=jiekuo1((k-1)*bs1+1:k*bs1);
                            jj2=jiekuo2((k-1)*bs1+1:k*bs1);
                            jff1(k)=sum(jj1);
                            jff2(k)=sum(jj2);

                            if jff1(k)>jff2(k)
                                rpn(k)=1;
                            else
                                rpn(k)=-1;
                            end
                        end
                        jf1(i)=sum(rpn.*pn);
                   end
                    [ZZ KK]=max(jf1.^2);
                    jf2=jf1(KK);
                    
                       if(jf2>0)
                          rec_ds(j)=1;
                       else
                          rec_ds(j)=0;
                       end
               end
           
               [en,er]=symerr(rec_ds,smag);
               en_ds(jj,ii)=en;
           end     
           
       % ****** RAKE接收 ******
           if flag_rake==1  
               jf=dsss_decode(NC,tc1,nt,pn,x,MM1,sr);
                   if(jf>0)
                       rec_rake(j)=1;
                   else
                       rec_rake(j)=0;
                   end
               
               [en,er]=symerr(rec_rake,smag);
               en_rake(jj,ii)=en;      
            end   
         end        
      end
  end
      if flag_ds==1
          err_ds=0;
         err_ds=sum(en_ds,1)/(ncyc*NM)
      end
     
      if flag_rake==1
          err_rake=0;
         err_rake=sum(en_rake,1)/(ncyc*NM)
      end
 
   figure(2),semilogy(SNR,err_ds,SNR,err_rake,'-*'),grid;
   xlabel('信噪比/dB');
   ylabel('误码率');
   legend('DSSS','DSSS-RAKE')%,sprintf('chip速率 = %0.5g chips/s',1/tc1),sprintf('chip长度 = %0.5g ',NC),3);    


