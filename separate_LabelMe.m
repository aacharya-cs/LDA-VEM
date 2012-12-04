function [] =  separate_LabelMe()

%% separates the visual words into categories
classnames = {'coast', 'forest', 'highway', 'insidecity', 'mountain', 'opencountry', 'street', 'tallbuilding'};
A = load('LabelMe_vwords.mat');

filenames = A.filename;
vwordsindex = A.vwordsindex;
vwordscount = A.vwordscount;
count = zeros(8,1);

for i=1:max(size(filenames))
 if(~isempty(regexp(filenames{i},'coast')))
   count(1) = count(1) + 1;
   Coast.filename{count(1)} = filenames{i}; 
   Coast.vwordsindex{count(1)} = vwordsindex{i}; 
   Coast.vwordscount{count(1)} = vwordscount{i};
 elseif(~isempty(regexp(filenames{i},'forest')))
   count(2) = count(2) + 1;   
   Forest.filename{count(2)} = filenames{i}; 
   Forest.vwordsindex{count(2)} = vwordsindex{i}; 
   Forest.vwordscount{count(2)} = vwordscount{i};
 elseif(~isempty(regexp(filenames{i},'highway')))
   count(3) = count(3) + 1;
   Highway.filename{count(3)} = filenames{i}; 
   Highway.vwordsindex{count(3)} = vwordsindex{i}; 
   Highway.vwordscount{count(3)} = vwordscount{i};
 elseif(~isempty(regexp(filenames{i},'insidecity')))
   count(4) = count(4) + 1;
   Insidecity.filename{count(4)} = filenames{i}; 
   Insidecity.vwordsindex{count(4)} = vwordsindex{i}; 
   Insidecity.vwordscount{count(4)} = vwordscount{i};
 elseif(~isempty(regexp(filenames{i},'mountain')))
   count(5) = count(5) + 1;
   Mountain.filename{count(5)} = filenames{i}; 
   Mountain.vwordsindex{count(5)} = vwordsindex{i}; 
   Mountain.vwordscount{count(5)} = vwordscount{i};
 elseif(~isempty(regexp(filenames{i},'opencountry')))
   count(6) = count(6) + 1;
   Opencountry.filename{count(6)} = filenames{i}; 
   Opencountry.vwordsindex{count(6)} = vwordsindex{i}; 
   Opencountry.vwordscount{count(6)} = vwordscount{i};
 elseif(~isempty(regexp(filenames{i},'street')))
   count(7) = count(7) + 1;
   Street.filename{count(7)} = filenames{i}; 
   Street.vwordsindex{count(7)} = vwordsindex{i}; 
   Street.vwordscount{count(7)} = vwordscount{i};
 elseif(~isempty(regexp(filenames{i},'tallbuilding')))
   count(8) = count(8) + 1;
   Tallbuilding.filename{count(8)} = filenames{i}; 
   Tallbuilding.vwordsindex{count(8)} = vwordsindex{i}; 
   Tallbuilding.vwordscount{count(8)} = vwordscount{i};
 else
   error('something is wrong!!');
 end 
end

A = Coast;
save('Coast_LabelMe_Vwords.mat','A');
A = Forest;
save('Forest_LabelMe_Vwords.mat','A');
A = Highway;
save('Highway_LabelMe_Vwords.mat','A');
A = Insidecity;
save('Insidecity_LabelMe_Vwords.mat','A');
A = Mountain;
save('Mountain_LabelMe_Vwords.mat','A');
A = Opencountry;
save('Opencountry_LabelMe_Vwords.mat','A');
A = Street;
save('Street_LabelMe_Vwords.mat','A');
A = Tallbuilding;
save('Tallbuilding_LabelMe_Vwords.mat','A');

end