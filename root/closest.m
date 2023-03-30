function indices = closest(s,V)

% Returns the indices of the elements in a vector V that are closest in
% numerical value to each of the elements in another vector s.
%
% indices = closest(s,V);

for i = 1:length(s)
    [~,indices(i)] = min(abs(V-s(i)));
end