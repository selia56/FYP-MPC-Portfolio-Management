function [p_estimate, q_estimate, p_est_min,diff_vec]=estimate_pq(e_Cest,N)

for kk=1:size(e_Cest,3)
    vec(:,kk)= sum(sum(e_Cest(:,:,kk)))/N;
end

diff_vec=abs(diff(vec));

q_estimate=1;
for ii=1:size(diff_vec,2)
    if diff_vec(:,ii)>0.59
        q_estimate=q_estimate+1;
    else
        break;
    end
end

[sortedValue_X , X_Ranked] = sort(sum(e_Cest(:,:,q_estimate)'),'ascend');
diff_sort=diff(sortedValue_X);

count=0;
for kk=1:size(diff_sort,2)
    if diff_sort(:,kk)<=diff_sort(:,1)
        count=count+1;
    end
end
count

[~,p_est_min]=min(sum(e_Cest(:,:,q_estimate)'));
p_est_method_1=X_Ranked(:,count+1)


if p_est_min < p_est_method_1
    p_estimate=p_est_min;
else p_estimate=p_est_method_1;
end 

end