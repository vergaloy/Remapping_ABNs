function [increasing,decreasing]=ABNs_responding_to_preS_and_postS(D,win)

for i=1:size(D,2)
    d=D{1,i};
    S(:,:,i)=Block_bootstrap(d,win); %get 10,000 bootstrap replicate, using a window size = 10 seconds.
end

%% HC vs pre
    Dif=S(:,:,2)-S(:,:,1); %Get surrogate differences. 
    P(:,1)=get_P_value(Dif); %Get the probability of including 0
    h=fdr_bh(P(:,1)); % correct by false-discovery rate
    ix=(h.*-(mean(Dif,2)<0))+(h.*(mean(Dif,2)>0)); % 1: increasing activity, -1: decreasing activity. 0:ns.
    IX(:,1)=ix;
 %% pre vs post   
        Dif=S(:,:,3)-S(:,:,2); %Get surrogate differences. 
    P(:,2)=get_P_value(Dif); %Get the probability of including 0
    h=fdr_bh(P(:,2)); % correct by false-discovery rate
    ix=(h.*-(mean(Dif,2)<0))+(h.*(mean(Dif,2)>0)); % 1: increasing activity, -1: decreasing activity. 0:ns.
    IX(:,2)=ix;

    
    
    
    
fprintf('ABNs increasing their activity in HC vs preS, and preS vs postS: \n');
increasing=sum(IX==1,1)
fprintf('ABNs increasing their activity in HC vs preS, and preS vs postS: \n');
decreasing=sum(IX==-1,1)

% plot stuff
X=[D{1,1},D{1,2},D{1,3}];temp=X(IX(:,1)==1,:);
temp=bin_data(temp,1,10);
temp=temp/max(temp(:));
mid(1)=size(D{1,1},2)/10;
mid(2)=mid(1)+size(D{1,2},2)/10;
[~,I]=sort(mean(temp(:,mid(1):mid(2)),2),'descend');
figure; imagesc(temp(I,:));colormap('hot');caxis([0 0.2]);xline(mid(1),'--c','LineWidth',2);xline(mid(2),'--c','LineWidth',2)


end

function S=Block_bootstrap(data,win)
parfor s=1:10000
    S(:,s)=mean(get_surrogate(data,win),2);
end
end


function t=get_surrogate(data,win)
x=floor(linspace(0,size(data,2),round(size(data,2)/win)+1));
n=size(x,2)-1;
s=[];
for i=1:n
    r=randperm(n,1);
    temp=data(:,x(r)+1:x(r+1));
    s=[s,temp];
end
t=mean(s,2);
end

function P=get_P_value(D)

L=size(D,2)/2;

for i=1:size(D,1)   
    dist=abs(sort(D(i,:)));
    if (mean(dist==0)<0.05) %% this is to assure with 95% accuaracy that the differences are not 0.      
        minv=min(dist);
        if (median(dist)>=0)
            ix     = find(dist == minv,1,'last');
        else
            ix     = find(dist == minv,1,'last');
        end
        P(i)=(L-abs(ix-L)+1)/(length(dist)/2+1);
    else
        P(i)=1;
    end
end
end



