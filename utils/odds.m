function odds(x,varargin)
% ODDS
% This function calculates the Risk Ratio and the Odds Ratio (OR) on a 2x2
% input matrix. Both ratios are computed with confidence intervals. If
% confidence interval of OR doesn't encompass the value OR=1, then the
% function computes the Bayesian Credibility Assessment of the test. If the
% test is credible, the function calculates the Association Parameter Phi.
% The association parameter Phi=sqrt(chisquare/N).
% The routine coumputes the Power and, if necessary, the sample sizes needed
% to achieve a power=0.80 using a modified asymptotic normal method with
% continuity correction as described by Hardeo Sahai and Anwer Khurshid in
% Statistics in Medicine, 1996, Vol. 15, Issue 1: 1-21.
% 
% Syntax: 	ODDS(X,ALPHA)
%      
%     Inputs:
%           X - 2x2 data matrix composed like this: 
% .............................................Cases...Controls
%                                              ___________
% Treated (or exposed to risk factor)          |  A  |  B  |
%                                             |_____|_____|
% Placebo (or not exposed to risk factors )    |  C  |  D  |
%                                             |_____|_____|
%                                                
%           ALPHA - Significance level (default=0.05).
% 
%     Outputs:
%           - Risk Ratio with Confidence interval.
%           - Absolute Risk Reduction.
%           - Relative Risk Reduction.
%           - Odds Ratio wirh Confidence interval
%           - Critical Odds Ratio (Bayesian Credibility Assessment)
%           - Phi (association parameter)
%           - Power and sample sizes calculation
% 
%      Example: 
% ..............................Cancer..Controls
%                                ___________
% Passive smoke exposed         |  25 |  21 |
%                               |_____|_____|
% Passive smoke not exposed     |  7  |  27 |
%                               |_____|_____|
% 
% Data matrix must be x=[25 21; 7 27];
% 
% Calling on Matlab the function: odds(x)
% answer is:
% 
% Significance level: 95%
%  
% Risk Ratio: 2.1021<2.6398<4.0597
% Absolute risk reduction: 33.8%
% Relative risk reduction: 62.1%
%  
% Odds Ratio: 1.6662<4.5918<12.6544
% Phi: 0.3149
% Moderate positive association (risk factor)
%  
% Bayesian Credibility Assessment
% Critical Odds Ratio: 2.4664
% OR>COR. Test is credible at the 95%
%  
% alpha = 0.0500  n1 = 46  n2 = 34
% Z1-b = 1.1924  Power (2-tails) = 0.8834
% 
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
% 
% To cite this file, this would be an appropriate format:
% Cardillo G. (2007) Odds: compute odds and risk ratio on a 2x2 matrix. 
% http://www.mathworks.com/matlabcentral/fileexchange/15347
%Input error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'real','finite','integer','nonnan','size',[2 2]}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
parse(p,x,varargin{:});
x=p.Results.x; alpha=p.Results.alpha;
clear p
x(x==0)=0.5;
fprintf('Significance level: %d%%\n', (1-alpha)*100)
disp(' ')
Za=-realsqrt(2)*erfcinv(2-alpha);
R=sum(x); %sum of the columns
p=x(1,:)./R; 
rr=p(1)/p(2); 
rrse=realsqrt(sum(1./x(1,:)-1./R)); %standard error of log(RR)
rrci=exp(reallog(rr)+([-1 1].*(Za*rrse))); %RR confidence interval
d=abs(diff(p)); %absolute risk reduction
rrr=d/p(2); %relative risk reduction
fprintf('Risk Ratio: %0.4f<%0.4f<%0.4f\n',rrci(1),rr,rrci(2))
if(rrci(1)<=1 && rrci(2)>=1)
   disp('Confidence interval encompasses RR=1. None significative association.')
end
fprintf('Absolute risk reduction: %0.1f%%\n',d*100)
fprintf('Relative risk reduction: %0.1f%%\n',rrr*100)
fprintf('Number Needed to Treat (NNT): %0.2f\n',1/d);
fprintf('Around %i patients need to be tested to correctly detect 100 positive tests for the presence of disease\n',ceil(100/d)) 
disp(' ')
%odd ratio (OR)
or=prod(diag(x))/prod(diag(rot90(x))); 
orse=realsqrt(sum(1./x(:))); %standard error of log(OR)
orci=exp(reallog(or)+([-1 1].*(Za*orse))); %OR confidence interval
fprintf('Odds Ratio: %0.4f<%0.4f<%0.4f\n',orci(1),or,orci(2))
if(orci(1)<=1 && orci(2)>=1)
   disp('Confidence interval encompasses OR=1. None significative association.')
else
    N=sum(x(:));
    Phi=(det(x)-N/2)/realsqrt(prod(sum(x))*prod(sum(x,2)));
    phi_hat=max(0,Phi^2-1/(N-1));
    k_hat=2-1/(N-1);
    V=sqrt(phi_hat/(k_hat-1));
    fprintf('Cramer''s V: %0.4f\n',V)
    switch sign(Phi)
        case -1
            txt2='negative association (protective factor)';
        case 1
            txt2='positive association (risk factor)';
    end
    V=abs(V);
    if V<=0.3
        txt1='Weak ';
    elseif (V>0.3 && V<=0.7)
        txt1='Moderate ';
    else
        txt1='Strong ';
    end
    disp([txt1 txt2])
    disp(' ')
    disp('Bayesian Credibility Assessment')
    orci=reallog(orci); 
    cor=exp(-diff(orci)^2/(4*realsqrt(prod(orci)))); %Critical odds ratio (COR)
    if or<1
        fprintf('Critical Odds Ratio: %0.4f\n',cor)
        if or<cor
            fprintf('OR<COR. Test is credible at the %d%%\n',(1-alpha)*100)
        else
            fprintf('OR>=COR. Test isn''t credible at the %d%%\n',(1-alpha)*100)
        end
    else
        cor=1/cor; %correct cor
        fprintf('Critical Odds Ratio: %0.4f\n',cor)
        if or>cor
            fprintf('OR>COR. Test is credible at the %d%%\n',(1-alpha)*100)
        else
            fprintf('OR<=COR. Test isn''t credible at the %d%%\n',(1-alpha)*100)
        end
    end
end
disp(' ')
%power (Asymptotic normal method)
k=R(2)/R(1);
q=1-p;
pm=(p(1)+k*p(2))/(k+1);
qm=1-pm;
Z1_b=(realsqrt(R(1)*d^2)-Za*realsqrt((1+1/k)*pm*qm))/realsqrt(p(1)*q(1)+p(2)*q(2)/k);
pwr=0.5*erfc(-Z1_b/realsqrt(2));
fprintf('alpha = %0.4f  n1 = %d  n2 = %d\n',alpha,R)
fprintf('Z1-b = %0.4f  Power (2-tails) = %0.4f\n',Z1_b,pwr)
if pwr<0.8
    %sample size (Modified Asymptotic normal method with continuity correction)
    nstar=(Za*realsqrt(pm*qm*(1+1/k))-realsqrt(2)*erfcinv(1.6)*realsqrt(p(1)*q(1)+p(2)*q(2)/k))^2/d^2;
    n1=round(nstar/4*(1+realsqrt(1+2*(k+1)/(k*d*nstar)))^2);
    n2=round(k*n1);
    disp(' ')
    disp('To achieve a recommended Power=0.80')
    fprintf('n1 = %d (add %d subjects to exposed row)\n',n1,n1-R(1))
    fprintf('n2 = %d (add %d subjects to not exposed row)\n',n2,n2-R(2))
end
