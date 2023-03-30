%--------------------------------------------------------------------------
function idx = findPeaksSeparatedByMoreThanMinPeakDistance(y,x,iPk,Pd)
% Start with the larger peaks to make sure we don't accidentally keep a
% small peak and remove a large peak in its neighborhood. 
if isempty(iPk) || Pd==0
  IONE = ones('like',getIZERO);
  idx = (IONE:length(iPk)).';
  return
end
% copy peak values and locations to a temporary place
pks = y(iPk);
locs = x(iPk);
% Order peaks from large to small
if coder.target('MATLAB')
    [~, sortIdx] = sort(pks,'descend');
else
    sortIdx = coder.internal.sortIdx(pks,'d');
end
locs_temp = locs(sortIdx);
idelete = ones(size(locs_temp))<0;
for i = 1:length(locs_temp)
  if ~idelete(i)
    % If the peak is not in the neighborhood of a larger peak, find
    % secondary peaks to eliminate.
    idelete = idelete | (locs_temp>=locs_temp(i)-Pd)&(locs_temp<=locs_temp(i)+Pd); 
    idelete(i) = 0; % Keep current peak
  end
end
% report back indices in consecutive order
idx = sort(sortIdx(~idelete));