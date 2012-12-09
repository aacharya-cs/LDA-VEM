function [A] = ss_feature_m(model)

for n=1:model.N
  smallphin = squeeze(model.smallphi(n,:,:));
  A(n,:) = model.sumzeta(n,:)*smallphin;
end

end