function [pvalue,KSSTAT]=kstestanaDDA(framenr,parameters,D, input,frametimerange,locerrorparam)
maxDfree = max(parameters(:,4),parameters(:,5));
maxrangeD = -log(1e-22)*maxDfree;
rangeD =maxrangeD/(input.precision*2):maxrangeD/input.precision:maxrangeD;
framescombined = zeros(length(rangeD),framenr);
for i = 1:numel(input.frametimerange)
input.frametime = input.frametimerange(i);
fractionframetimerange = length(D(:,D(3,:)==input.frametime))./length(D);

if input.confinement == true
[fx,fy] = Generateconfinedfunction(0:0.05:5,input);
fx = fx';
fy = fy';
else
    fx = NaN;
    fy = NaN;
end
if input.trackingwindow < 100
      maxD = (input.trackingwindow*input.pixelsize)^2/(4*input.frametime);
      maxDindtracking = round(maxD./(rangeD(2)-rangeD(1)));
else
    maxDindtracking = 0;
end
for ii = 1:size(parameters)
koff = parameters(ii,2);
kon = parameters(ii,3);
Dfree = max(parameters(ii,4),parameters(ii,5));
c = parameters(ii,1);
D1 = min(parameters(ii,4),parameters(ii,5));
if input.fitlocerror == 0
    locerror = input.sigmaerror.^2/input.frametime;
else
    locerror = locerrorparam.^2/input.frametime;
end
framescombinedtemp = DDistributiongenerator(koff,kon,Dfree,D1,rangeD,locerror,fx,fy,maxDindtracking,input,1);
framescombinedtemp = framescombinedtemp./sum(framescombinedtemp);
framescombinedtemp = fractionframetimerange*c*framescombinedtemp(:,framenr);
framescombined = framescombined + framescombinedtemp;
end
end
D = D(1,D(2,:)==framenr);
%framescombined = [framescombined(1)/2; framescombined];
func = @(x) interp1(framescombined,x,'spline');
D = sort(D);
cumtrapzframescombined = cumtrapz(framescombined(:,framenr));
cumtrapzframescombined = cumtrapzframescombined + 1-cumtrapzframescombined(end);
Dconverted = (D-rangeD(1))*framenr./(rangeD(3)-rangeD(2))+1;

number =  interp1(cumtrapzframescombined,Dconverted,'spline');
% for i = 1:numel(Dconverted)
% number(i) = integral(func,0.5,Dconverted(i));
% end
number(Dconverted<0.5) =0;
number = sort(number);
try
[ans,pvalue,KSSTAT]=kstest(D,'CDF',[D' number']);
catch
    keyboard
end
% figure
% hold on
% plot(D,number)
% plot(D,1/numel(number):1/numel(number):1)
% figure
% plot(diff([number; 1/numel(number):1/numel(number):1]))

end