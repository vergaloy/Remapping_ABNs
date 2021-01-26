function [clus,Z,HM]=significant_linkage(datain,varargin)

%% INTIALIZE VARIABLES
inp = inputParser;
valid_v = @(x) isnumeric(x);
valid_c = @(x) ischar(x);
addRequired(inp,'datain',valid_v) 
addParameter(inp,'alpha',0.05,valid_v)  %alpha
addParameter(inp,'Cdist','correlation')  % Distance metric
addParameter(inp,'Cmethod','average',valid_c )  % Clustering method


inp.KeepUnmatched = true;
parse(inp,datain,varargin{:});
warning off
Cdist=inp.Results.Cdist;
Cmethod=inp.Results.Cmethod;
alpha=inp.Results.alpha;



%% Get cluster strength
HM=get_dist(datain,Cdist);
Z=linkage(datain,Cmethod,Cdist);
thr=get_thr(datain,Cmethod,Cdist,alpha);
S=Z(:,3)<thr;
c = cophenet(Z,squareform(1-HM,'tovector'));
fprintf('Cophenetic correlation coefficient: %1.2f \n',c);
 L=find_leaves_in_node(Z);
[clus,~]=get_clusters(S,L); 
end
%%

function out=get_dist(in,distance)
    out=squareform(1-pdist(in,distance))+diag(ones(1,size(in,1)));
    out(isnan(out))=0;
end


function thr=get_thr(datain,Cmethod,Cdist,alpha)
for s=1:10000
    temp=shuffle_row(datain);
    temp=linkage(temp,Cmethod,Cdist);
    S(:,s)=temp(:,3);
end
thr=prctile(S,alpha,2);
end

function L=find_leaves_in_node(Z)
M=size(Z,1)+1;
L = num2cell(Z(:,1:2),2);
for i=1:size(Z,1)
    temp=L{i, 1};
    while (ismember(1,temp>M))
        f=find(temp>M,1);
        k=temp(f);
        temp(f)=[];
        temp=[temp,L{k-M, 1}];
    end
    L{i, 1}=temp;
end
end

function [clus,C]=get_clusters(S,L)
C=zeros(1,size(L,1));
for i=length(S):-1:1
    if (S(i))
        for k=1:size(L,1)
            if (ismember(L{k},L{i}))
                if (k~=i)
                    L{k}=[];
                    C(k)=i;
                end
                C(k)=i;
            end
        end
    else
        L{i}=[];
    end
end
clus=L;
clus=clus(~cellfun('isempty',clus));
u=unique(C);
u(u==0)=[];
for i=1:length(u)
    C(u(i)==C)=i;
end
end

function out=shuffle_row(in)
out=zeros(size(in,1),size(in,2));
for i=1:size(in,1)
    out(i,:)=in(i,randperm(size(in,2),size(in,2)));   
end
end










