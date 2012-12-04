function [] = disp_accuracy()

addpath(genpath('/home/ayan/Dropbox/aaaa/myLDA/'));

%% add paths of the saved .mat files

% % savepaths{1} = ['/home/ayan/Dropbox/confdata_17_May_2012/'];
% % savepaths{2} = ['/home/ayan/Dropbox/confdata_18_May_2012/'];
% % savepaths{3} = ['/home/ayan/Dropbox/confdata_19_May_2012/'];
savepaths{1} = ['/v/filer4b/v16q001/data/Ayans/confdata_31_May_2012/'];
%%savepaths{2} = ['/home/ayan/Dropbox/aaaa/confdata_25_May_2012/'];
%%savepaths{3} = ['/home/ayan/Dropbox/aaaa/confdata_26_May_2012/'];

for i=1:size(savepaths,2)
    addpath(genpath(savepaths{i}));
end

p1mat = [0.1:0.1:1];
p2mat = [0.1:0.1:1];
k2mat = [20:10:150];
prefixfilename = '/v/filer4b/v16q001/data/Ayans/Results_finalbigrun1/';
otherindex = 'finalbigrun1';
minvtopic = 1;
classfilenameprefix = '/lusr/u/ayan/Documents/aaaa/myLDA/conferencenames';
classnum = 6;

option = 4;
troption = 0;
svmoptionval = 4;
svmcval = 0.1;

%savedstat    = zeros(length(p1mat), length(p2mat), length(k2mat), 5);
%savedstatstd = zeros(length(p1mat), length(p2mat), length(k2mat), 5);

savedstat    = zeros(length(p1mat), length(k2mat), 5);
savedstatstd = zeros(length(p1mat), length(k2mat), 5);
mdirstr = ['mkdir ' prefixfilename(1:end-1)];

if(~exist(prefixfilename(1:end-1)))
    system(mdirstr);
end

for arg1=1:length(p1mat)
    %for arg2=1:length(p2mat)
    for arg3=1:length(k2mat)
        
        filename = [prefixfilename num2str(p1mat(arg1)) '_' num2str(p2mat(arg1)) '_' num2str(k2mat(arg3)) '.txt'];
        fp = fopen(filename,'w');
        
        WFacc = zeros(5,classnum);
        count = 0;
        indcount = [];
        
        for i=1:classnum
            
            classfilename = [classfilenameprefix num2str(i) '.txt'];
            B = textread(classfilename, '%s');
            
            tempfile = zeros(classnum,1);
            for ss=1:size(savepaths,2)
                savefilename{ss} = [savepaths{ss} B{end} '_' num2str(minvtopic) '_' num2str(p1mat(arg1)) '_' num2str(p2mat(arg1)) '_' num2str(k2mat(arg3)) '_' num2str(option) '_' num2str(troption) '_' num2str(svmoptionval)];
                savefilename{ss} = [savefilename{ss} '_' num2str(svmcval) '_' otherindex '.mat'];
                tempfile(ss) = exist(savefilename{ss});
            end
            sum(tempfile);
            if(sum(tempfile)>0)
                
                indcount = [indcount i];
                sfilename = savefilename{find(tempfile>0)};
                A = load(sfilename);
                savefilename;
                count = count + 1;
                fprintf(fp,'%s\n','                SVMTr    SVMT    SVMP MLDABTr MLDABT MLDABP DSLDAWLTr DSLDAWLT DSLDAWLP MLDATr   MLDAT MLDAP DSLDATr DSLDAT DSLDAP');
                p = [];
                
                %% SVM
                temp = A.svmperf.F_score_train; temp(isnan(temp)) = 0; p = [p temp];
                temp = A.svmperf.F_score_test; temp(isnan(temp)) = 0; p = [p temp];
                p = [p diag(A.svmperf.confmat_test)./sum(A.svmperf.confmat_test,2)];
                temp = A.svmperf.F_score_test; temp(isnan(temp)) = 0;
                wfscore1 = sum(temp.*sum(A.svmperf.confmat_test,2))/sum(sum(A.svmperf.confmat_test,2));
                wfscore1acc = sum(diag(A.svmperf.confmat_test))/sum(sum(A.svmperf.confmat_test));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% one vs all (non-transfer)
                Facctr = zeros(classnum,1);
                Facctest = zeros(classnum,1);
                Pacctr = zeros(classnum,1);
                Pacctest = zeros(classnum,1);
                wacctr = 0;
                wfscore2 = 0;
                wfscore2acc = 0;
                for ii=1:classnum
                    temp1 = A.trperfmedLDAova{ii}.accval;
                    temp2 = A.testperfmedLDAova{ii}.accval;
                    temp1(isnan(temp1)) = 0;
                    temp2(isnan(temp2)) = 0;
                    Facctr(ii,1) = temp1;
                    Facctest(ii,1) = temp2;
                    Pacctr(ii,1) = A.trperfmedLDAova{ii}.confmat(1,1)/sum(A.trperfmedLDAova{ii}.confmat(1,:));
                    Pacctest(ii,1) = A.testperfmedLDAova{ii}.confmat(1,1)/sum(A.testperfmedLDAova{ii}.confmat(1,:));
                    wacctr = wacctr + temp1*sum(A.trperfmedLDAova{ii}.confmat(1,:));
                    wfscore2 = wfscore2 + temp2*sum(A.testperfmedLDAova{ii}.confmat(1,:));
                    wfscore2acc = wfscore2acc + sum(diag(A.testperfmedLDAova{ii}.confmat))/sum(sum(A.testperfmedLDAova{ii}.confmat));
                end
                wfscore2 = wfscore2/sum(sum(A.testperfmedLDAova{1}.confmat,2));
                wfscore2acc = wfscore2acc/classnum; 
                p = [p Facctr];
                p = [p Facctest];
                p = [p Pacctest];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% DSLDA with no shared latent topic
                temp = A.trperfLLDA.accval; temp(isnan(temp)) = 0; p = [p temp];
                temp = A.testperfLLDA.accval; temp(isnan(temp)) = 0; p = [p temp];
                p = [p diag(A.testperfLLDA.confmat)./sum(A.testperfLLDA.confmat,2)];
                temp = A.testperfLLDA.accval; temp(isnan(temp)) = 0;
                wfscore3 = sum(temp.*sum(A.testperfLLDA.confmat,2))/sum(sum(A.testperfLLDA.confmat,2));
                wfscore3acc = sum(diag(A.testperfLLDA.confmat))/sum(sum(A.testperfLLDA.confmat,2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% medLDA
                temp = A.trperfmedLDA.accval; temp(isnan(temp)) = 0; p = [p temp];
                temp = A.testperfmedLDA.accval; temp(isnan(temp)) = 0; p = [p temp];
                p = [p diag(A.testperfmedLDA.confmat)./sum(A.testperfmedLDA.confmat,2)];
                temp = A.testperfmedLDA.accval; temp(isnan(temp)) = 0;
                wfscore4 = sum(temp.*sum(A.testperfmedLDA.confmat,2))/sum(sum(A.testperfmedLDA.confmat,2));
                wfscore4acc = sum(diag(A.testperfmedLDA.confmat))/sum(sum(A.testperfmedLDA.confmat,2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% DSLDA
                temp = A.trperf.accval; temp(isnan(temp)) = 0; p = [p temp];
                temp = A.testperf.accval; temp(isnan(temp)) = 0; p = [p temp];
                p = [p diag(A.testperf.confmat)./sum(A.testperf.confmat,2)];
                temp = A.testperf.accval; temp(isnan(temp)) = 0;
                wfscore5 = sum(temp.*sum(A.testperf.confmat,2))/sum(sum(A.testperf.confmat,2));
                wfscore5acc = sum(diag(A.testperf.confmat))/sum(sum(A.testperf.confmat,2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for j=1:size(p,1)
                    fprintf(fp, '%9s\t', B{j});
                    for jj=1:size(p,2)
                        fprintf(fp, '%.5f\t', p(j,jj));
                    end
                    fprintf(fp, '\n');
                end
                
                fprintf(fp,'\n weighted F-score: SVMtest: %f nontransfermedlDAtest: %f DSLDA without latent shared: %f medLDAtest: %f DSLDATest: %f\n\n', wfscore1, wfscore2, wfscore3, wfscore4, wfscore5);
                fprintf(fp,'\n weighted accuracy: SVMtestacc: %f nontransfermedlDAtestacc: %f DSLDA without latent shared: %f medLDAtestacc: %f DSLDATestacc: %f\n\n', wfscore1acc, wfscore2acc, wfscore3acc, wfscore4acc, wfscore5acc);
                WFacc(:,i) = [wfscore1 wfscore2 wfscore3 wfscore4 wfscore5]';
                Wacc(:,i) = [wfscore1acc wfscore2acc wfscore3acc wfscore4acc wfscore5acc]';
            end
        end
        if(count>0)
            temp = mean(WFacc(:,indcount),2);
            stdtemp = std(WFacc(:,indcount),0,2);
            tempacc = mean(Wacc(:,indcount),2);
            stdtempacc = std(Wacc(:,indcount),0,2); 
            %%savedstat(arg1,arg2,arg3,:) = temp';
            %%savedstatstd(arg1,arg2,arg3,:) = stdtemp';
            savedstat(arg1,arg3,:) = tempacc';
            savedstatstd(arg1,arg3,:) = stdtempacc';
        end
        fclose(fp);
    end
    %end
end


save(otherindex,'savedstat');
save([otherindex '_std'],'savedstatstd');

end

