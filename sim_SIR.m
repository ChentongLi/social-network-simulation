% this code is coded by Prof. Junling Ma in 2014
% modified by Chentong Li in 2015
function x=sim_sir()
% create random graph
n=1000; % #nodes
degrees = poissrnd(4,n,1);
edges=[];
for i=1:n
    edges=[edges;ones(degrees(i),1)*i];
end
L=length(edges);
net=sparse(n,n);
failed=0;
while L>1 && failed<100
    i=floor(rand()*(L-1))+1;
    j=floor(rand()*(L-1))+1;
    A=edges(i);
    B=edges(j);
    if i~=j && A~=B && net(A,B) == 0
        net(A,B)=1;
        net(B,A)=1;
        edges([i,j])=[];
        L=L-2;
        failed=0;
    else
        failed=failed+1;
    end
end
%net
beta=0.1;
gamma=0.2;
I0=10;
degrees =full(sum(net));
N=length(degrees);
time=0;
nodes=floor(rand(1,I0)*(N-1))+ones(1,I0);
susp=ones(1,N);
susp(nodes)=0;
Tinf=exprnd(1./(degrees(nodes)*beta));
Trec=exprnd(1/gamma,1,length(nodes));
x=[time,length(nodes)];
while time<1000 && ~isempty(nodes) 
    [next_inf , ind_inf]=min(Tinf);
    [next_rec , ind_rec]=min(Trec);
    node_inf = nodes(ind_inf);
    node_rec = nodes(ind_rec);
    if next_inf < next_rec
        time = next_inf;
        neighbors = find(net(node_inf,:));
        ind = floor(rand()*(length(neighbors)-1))+1;
        qinf = neighbors(ind); 
        if susp(qinf)
            nodes = [nodes, qinf];
            susp(qinf)=0;
            Tinf=[Tinf,exprnd(1/(degrees(qinf)*beta))+time];
            Trec=[Trec,exprnd(1/gamma)+time];
        end
        Tinf(ind_inf)=exprnd(1/(degrees(node_inf)*beta))+time;
    else
        nodes(ind_rec)=[];
        Tinf(ind_rec)=[];
        Trec(ind_rec)=[];
        time=next_rec;
    end
    x=[x;time, length(nodes)];
    [time, length(nodes)]
end


















