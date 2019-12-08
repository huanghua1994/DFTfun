function [x, w] = mycheb2(N)
	x = zeros(N, 1);
	w = zeros(N, 1);
	for i = 1 : N
		x(i) = cos(i * pi / (N + 1));
		w(i) = pi / (N + 1) * (1 - x(i)) * (1 + x(i));
	end
end