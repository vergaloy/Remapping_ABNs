function O=get_activity_vectors_consolidation(X,N,hyp,div)
O=[];
for i=1:max(N)
    temp=X(N==i,:);
    h=hyp(i,:);
    o_temp=get_sleep_binned(temp,h,div);
    O=cat(1,O,o_temp);
end
end




function B=means_binned_data(in,bsize)
lin=round(linspace(1,size(in,2),bsize+1));
for i=1:size(lin,2)-1
B(:,i)=nanmean(in(:,lin(i):lin(i+1)),2);   
end    
end



function O=get_sleep_binned(temp,h,div)
S=get_activity_v(temp,h<0.5,div);
V=get_activity_v(temp,h>=0.5,div);
A=get_activity_v(temp,h>3,div);
O=cat(3,S,V,A);
end

function T=get_activity_v(temp,h,div)
T=temp;
T(:,h)=nan;
T=means_binned_data(T,div);
end

