function [] = class_topic_correlation(model, data, vocabList, tagList, maxtopic)

k1 = size(tagList);

for i=1:data.Y
 ind = find(data.classlabels==i);
 stopicassignment = (mean(model.phi(ind,:),[],1))'; 
 topicind = [1:model.K]';
 stopicmat = [stopicind topicassignment];
 [~, index] = sort(topicmat(:,2));
 topicmatsort = topicmatsort(index,:);
 domtopicind  = topicmatsort(1:maxtopic,1); 
 model.beta()
end

end
