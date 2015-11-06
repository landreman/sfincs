function y = sfincs_xWeight(xx)

global xGrid_k

y = (xx.^xGrid_k) .* exp(-xx.*xx);

end