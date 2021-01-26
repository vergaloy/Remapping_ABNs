function clus=cluster_activity_vectors(a)
a=a(:,[1,2,3,5]); % we remove the consolidation period from the analysis. 
rng('default')
a(a>1)=1;
ix=~(sum(a,2)==0); 
a=a(ix,:);


[clus,Z,HM]=significant_linkage(a','Cdist','cosine','Cmethod','ward');

subplot(1,2,1);imagesc(HM);subplot(1,2,2);dendrogram(Z);

end



