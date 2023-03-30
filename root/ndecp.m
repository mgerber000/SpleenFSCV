function b = ndecp(a,n)
% Converts numbers in a vector to numbers with a set number of decimal places.
% B = ndecp(A,n)
%
% example:
% b = ndecp([pi,2*pi,3*pi],3)
% b =
%     3.1420    6.2830    9.4250

b = nan(size(a));
for i = 1:length(b)
    b(i) = round(a(i)*10^n)/10^n;
end
