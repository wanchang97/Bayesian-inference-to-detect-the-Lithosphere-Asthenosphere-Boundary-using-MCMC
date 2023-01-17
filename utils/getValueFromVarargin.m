function [val,left] = getValueFromVarargin(key,vars,defaultValue)
% [value,argsleft] = getValueFromVarargin(key,varargin,defaultValue)
if nargin<3
   defaultValue = [];
end
left = cell(0,0);
val = defaultValue;
k=1;
for I = 1:2:length(vars)
   if strcmpi( vars{I}, key )
      val = vars{I+1};
   else
      left{k} = vars{I}; 
      left{k+1} = vars{I+1};
      k = k + 2;
   end
end
