function y = Heaviside(x, b, eta)
y = (tanh(b*eta) + tanh(b*(x - eta)))./(tanh(b*eta) + tanh(b*(1 - eta)));
end
