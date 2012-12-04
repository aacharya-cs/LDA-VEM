function value = h_inv_g_alpha(alpha, g, k, lambda)

  
  qdiag=-k*[psi(1,alpha)+lambda];
  z=k*psi(1,sum(alpha));
  b=[sum(g./qdiag)]/[1/z+sum(1./qdiag)];
  value=(g-b)./qdiag;
