%Generate null distribution for principal angles

D = 15;
N = 5000;
cos_theta = zeros(N,1);

for i = 1:N
	%Generate a and b
	a = normrnd(1,1, 1,D);
	b = normrnd(1,1, 1,D);
	%Save dot product
	cos_theta(i) = sum(a*b')/sqrt(sum(a*a')*sum(b*b'));
end

%Take inverse to get to angle in radians
acos_theta = acos(cos_theta);

%For all angles greater than pi/2, flip them to be less than pi/2
gtpi2 = acos_theta > (pi/2);
acos_theta(gtpi2) = pi - acos_theta(gtpi2);
acos_theta = sort(acos_theta);
%Find cutoff for which less than 5% lie
thr_idx = int16(N*0.05);
thr = acos_theta(N-thr_idx)

%Gives a threshold of 1.38
