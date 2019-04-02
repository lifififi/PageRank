clear all;

fid = fopen('input_graph.txt');
EDGES = fscanf(fid, '%i -> %i', [2 inf]);
fclose(fid);

d=0.86;
zbieznosc = 0.0001;
N=max(max(EDGES));

B=zeros(N);

for i=1:1:N
    B(EDGES(2,i),EDGES(1,i))=1;
end

L=zeros(1,N);

for i=1:1:N
   L(EDGES(1,i))=L(EDGES(1,i))+1;
end

for i=1:1:N
    if(L(i)==0)
        L(i)=N-1;
    end
end
LAMBDA=diag(1./L);
E=ones(N);
M=(d*B*LAMBDA+((1-d)/N)*E);

R0=zeros(N,1);

for i=1:1:N
R0(i)=1/N;
end

lastR=ones(N,1)*inf;
R=R0;
i=0;
Rfinal= [    0.0205
    0.0142
    0.4445
    0.4085
    0.0142
    0.0268
    0.0142
    0.0142
    0.0142
    0.0142
    0.0142];
while (norm((R-lastR),2)>zbieznosc)
i=i+1;
lastR=R;
R=M*R;
R=R/norm(R,1);
lambda=R.'*M*R/(R.'*R);
power_err(i)=norm(R-Rfinal);
end
stem(R);
grid on;
xlabel('numer strony');
ylabel('prawdopodobienstwo');
title('PageRank');
figure;

power_err_uciety=power_err(1:35);

semilogy(power_err_uciety);
grid on;
title('zbieznosc');

polyfit(1:length(power_err_uciety),log(power_err_uciety), 1)



M_w_wl=sort(abs(eig(M)),'descend');
log(M_w_wl(2)/M_w_wl(1))
