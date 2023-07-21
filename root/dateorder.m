function order = dateorder(list)

numlist = str2double(list);
order = [];
while ~(size(list,1) == size(order,2))
   ind = double(numlist(:,3) == min(numlist(:,3)));
   ind(ind == 0) = NaN;
   i = numlist(:,2) .* ind;
   
   ind = double(i == min(i));
   ind(ind == 0) = NaN;
   i = numlist(:,1) .* ind;
   
   ind = double(i == min(i));
   ind(ind == 0) = NaN;
   i = numlist(:,4) .* ind;
   
   ind = double(i == min(i));
   ind(ind == 0) = NaN;
   i = numlist(:,5) .* ind;
   
   ind = double(i == min(i));
   ind(ind == 0) = NaN;
   i = numlist(:,6) .* ind;
   
   lowestsecond = find(i == min(i));
   
   order = [order,lowestsecond'];
   numlist(lowestsecond',:) = NaN;
end