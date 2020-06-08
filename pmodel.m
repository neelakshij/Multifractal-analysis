%{ 
   Credit: Neelakshi Joshi 
   This code is written by Dr. Neelakshi Joshi for her thesis work at National Institute for Space Research (INPE), Brazil.
   Advisors: Dr. Reinaldo Rosa & Dr. Stephan Stephany.

   This code generates anaytical "P model singularity spectrum" based on the generalized two-scale weighted Cantor set.
   References: 1) Meneveau and Sreenivasan, Physical Review Letters, 59, 1987.
	       2) Halsey et al., Physical Review A, 33, 2, 1986.


   Inputs: 
	p1 = probability parameter, 0 < p1 < 1
	l1 = cascading length parameter, 0 < l1 < 1
	dp = loss of probability where p1+p2 < 1 and p1+p2+dp = 1
   Outputs:
	alpha and falpha = singularity spectrum

   Implementation: [alpha,falpha]=fit_pmodel_halsey(0.4,0.4,0.01)
   If no "dp" parameters is provided then default value is dp = 0.0.
   Value for l1 can be different from 0.5. 
%}

function [alpha,falpha]=pmodel(p1,l1,dp)

	if nargin < 3
	   dp = 0.0;
	end
 
	p2 = (1-p1) - dp; 
	l2 = 1-l1;
	q = -20:0.1:20;
	alpha = [];
	falpha = [];

	weight = p2/p1;
	nbym = 1 + (weight).^q;
	mbyn = 1./nbym;
	bt = (weight).^q;

 	for i=1:length(bt)
   	   alpha(i)=(log(p1) + (bt(i)*log(p2)))/(log(l1) + bt(i)*log(l2));
    	   falpha(i)=(bt(i)*log(bt(i)) - (bt(i)+1)*log(bt(i)+1))/(log(l1) + bt(i)*log(l2));
	end
 	
	if (log(p1)/log(l1)) < (log(p2)/log(l2))
     	   amin = log(p1)/log(l1);
           amax = log(p2)/log(l2);
        else
           amax = log(p1)/log(l1);
           amin = log(p2)/log(l2);
        end


	figure; plot(alpha,falpha,'-ok','linewidth',1,'markersize',5);
	title({'P model spectrum', ['p1 =' num2str(p1),' l1=',num2str(l1),' dp=' num2str(dp)]});
	xlabel('\alpha'); ylabel('f(\alpha)');
