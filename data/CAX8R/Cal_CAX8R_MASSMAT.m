clc
clear
R = [0 100 100 0 50 100 50 0];
% Xg = [0];
% Xg = [-0.5773502691896,0.5773502691896];
% W = [1 1];
Xg = [-0.774597 0 0.774597];
W = [0.5555555555555555 0.888888888888888 0.5555555555555555];
% Xg = [-0.8611363115941 -.3399810435849 0.3399810435849 0.8611363115941];
NtN = { };
NtNsum =zeros(16,16);
for i = 1:length(Xg)
    g = Xg(i);
    w1 = W(i);
    for j = 1:length(Xg)
     h = Xg(j);
     w2 = W(j);
     gp1 = 1+g;
	 gm1 = 1-g;
	 hp1 = 1+h;
	 hm1 = 1-h;
	 gp1ph = 1+g+h;
	 gp1mh = 1+g-h;
	 gm1ph = 1-g+h;
	 gm1mh = 1-g-h;
	
	N(1,1) =  -0.25*gm1*hm1*gp1ph;
	N(1,3) =  -0.25*gp1*hm1*gm1ph;
	N(1,5) =  -0.25*gp1*hp1*gm1mh;
	N(1,7) =  -0.25*gm1*hp1*gp1mh;
	N(1,9) =  0.5*gm1*gp1*hm1;
	N(1,11) =  0.5*hm1*hp1*gp1;
	N(1,13) =  0.5*gm1*gp1*hp1;
	N(1,15) =  0.5*hm1*hp1*gm1;
    N(2,2) =  -0.25*gm1*hm1*gp1ph;
	N(2,4) =  -0.25*gp1*hm1*gm1ph;
	N(2,6) =  -0.25*gp1*hp1*gm1mh;
	N(2,8) =  -0.25*gm1*hp1*gp1mh;
	N(2,10) =  0.5*gm1*gp1*hm1;
	N(2,12) =  0.5*hm1*hp1*gp1;
	N(2,14) =  0.5*gm1*gp1*hp1;
	N(2,16) =  0.5*hm1*hp1*gm1;
    
    r = 0;
    for rid = 1:length(R)
        r = r + R(rid)*N(2,rid*2);
    end
    scale = 2*pi*r*w1*w2*2500*7900;
    NtN{length(Xg)*(i-1)+j} = scale*N'*N;
    NtNsum = NtNsum + NtN{length(Xg)*(i-1)+j};
    end
end
NtN_col = sum(NtNsum,2);
eigNtN = eig(NtNsum)
    