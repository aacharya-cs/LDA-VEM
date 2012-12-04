function [] = converttoSLDA(trdata, testdata)

Ntrain = max(size(trdata.windex));
Ntest  = max(size(testdata.windex));

fp = fopen('/home/ayan/Dropbox/DSLDA_SDM/onlineHDP/trdatasz.txt','w+');
fprintf(fp, '%d', Ntrain);
fclose(fp);

fp = fopen('/home/ayan/Dropbox/DSLDA_SDM/onlineHDP/testdatasz.txt','w+');
fprintf(fp, '%d', Ntest);
fclose(fp);

fp = fopen('/home/ayan/Dropbox/DSLDA_SDM/onlineHDP/trdataSLDA.txt','w+');
for i=1:Ntrain
    fprintf(fp, '%s\n', convertTosLDAformat(trdata.windex{i}, trdata.wcount{i}));
end
fclose(fp);

fp = fopen('/home/ayan/Dropbox/DSLDA_SDM/onlineHDP/trdataSLDALabels.txt','w+');
fprintf(fp, '%d\n', trdata.classlabels-1);
fclose(fp);

fp = fopen('/home/ayan/Dropbox/DSLDA_SDM/onlineHDP/testdataSLDA.txt','w+');
for i=1:Ntest
    fprintf(fp, '%s\n', convertTosLDAformat(testdata.windex{i}, testdata.wcount{i}));
end
fclose(fp);

fp = fopen('/home/ayan/Dropbox/DSLDA_SDM/onlineHDP/testdataSLDALabels.txt','w+');
fprintf(fp, '%d\n', testdata.classlabels-1);
fclose(fp);

end
