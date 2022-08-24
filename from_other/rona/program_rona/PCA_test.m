% load hald
ingredients = sts;%[observations#, variables#] %aim to find the component of the variables
% figure(333);hold on;
% for i = 1:length(ingredients)
%    plot(ingredients(:,i))
% end
% [coeff,score] = pca(ingredients);

new_ingredients=ingredients-repmat(mean(ingredients),size(ingredients,1),1);
covariancematrix=cov(new_ingredients);
[V,D] = eig(covariancematrix);
V; %eigenvector
D = diag(D) %eigenvalue

maxeigval=V(:,length(D)-3:length(D));
finaldata=maxeigval'*(new_ingredients)';
figure;hold on;plot(finaldata(1,:));plot(finaldata(2,:));%plot(finaldata(3,:));plot(finaldata(4,:));
figure;hold on;plot(finaldata(2,:),finaldata(1,:),'*');
% figure;hold on
% for i=1:size(finaldata,2)
%     if  finaldata(i)>=0
%         plot(x(i),y(i),'o')
%         plot(x(i),y(i),'r*')
%         
%     else
%         plot(x(i),y(i),'o')
%         plot(x(i),y(i),'g*')
%     end
%     
% end
