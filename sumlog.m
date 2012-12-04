function [v]=sumlog(loga,logb)

if(logb==0)
 v=loga;
elseif (loga==0)
 v=logb;   
elseif (loga < logb)
      v = logb+log(1 + exp(loga-logb));
else
      v = loga+log(1 + exp(logb-loga));
end




