%{ 
   Multifractal Detrended Fluctuation Analysis (MFDFA)
   Reference: Kantelhardt et al., Physica A, 316, 1, 2002.

   Credit: IHLEN, E., Frontiers in Physiology, 3, 2012.
   Adapted by Neelakshi Joshi.
   
   Inputs: 
	data = timeseries (1-d array)
	m = detrending order (int)
        scale = scales (1-d array)
        q = moments (1-d array)
             needCorrection parameter decides how to define profile to implement the MFDFA, 
             according to the conditions described by Ihlen, 2012 in table 1.
	needCorrection = 0 then profile = cumsum(data - mean(data)
	needCorrection = 0 then profile = diff(data)	
	needCorrection = value other than 0 and 1 then profile = data
   
   Outputs:
	Fq = qth order Fluctuation function
	hq = the generalized Hurst exponents
	tq = the classical multifractal exponents
	alpha and falpha = multifractal exponents
%}  

 
function [Fq, hq, tq, alpha, falpha] = mfdfa1d(data,m,scale,q,needCorrection)
    if size(data,2)==1;
       data=transpose(data);
    end

    if needCorrection==1
       X = cumsum(data-mean(data));
    elseif needCorrection==0
       X = diff(data);
    else
       X = data;
    end
 
    for ns=1:length(scale),
        segments(ns)=floor(length(X)/scale(ns));
        for v=1:segments(ns),
            Index=((((v-1)*scale(ns))+1):(v*scale(ns)));
            C=polyfit(Index,X(Index),m);
            fit=polyval(C,Index);
            RMS_scale{ns}(v)=sqrt(mean((X(Index)-fit).^2));
        end
        for nq=1:length(q),
            qRMS{nq,ns}=RMS_scale{ns}.^q(nq);
            Fq(nq,ns)=mean(qRMS{nq,ns}).^(1/q(nq));
        end
        Fq(q==0,ns)=exp(0.5*mean(log(RMS_scale{ns}.^2)));
    end

    for nq=1:length(q),
        C = polyfit(log2(scale),log2(Fq(nq,:)),1);
        if needCorrection==1
           hq(nq) = C(1)+1;
        else
           hq(nq) = C(1);
        end
        qRegLine{nq} = polyval(C,log2(scale));
    end

    tq = hq.*q-1;
    alpha = diff(tq)./(q(2)-q(1));
    falpha = (q(1:end-1).*alpha)-tq(1:end-1);
