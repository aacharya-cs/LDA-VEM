function [medLDAperf] = getmedLDAperf(trdata, testdata)


%%trdata
%%testdata
medLDAperf = [];

disp('inside getmedLDAperf');

trfilename   = '/lusr/u/ayan/MLDisk/medlda_files/medLDAtrainfile';
testfilename = '/lusr/u/ayan/MLDisk/medlda_files/medLDAtestfile';

fp1 = fopen(trfilename, 'w+');
fp2 = fopen(testfilename, 'w+');

trfilenameslda   = '/lusr/u/ayan/MLDisk/medlda_files/sLDAtrainfile';
trclassesslda   = '/lusr/u/ayan/MLDisk/medlda_files/sLDAtrainlabels';
testfilenameslda = '/lusr/u/ayan/MLDisk/medlda_files/sLDAtestfile';
testclassesslda   = '/lusr/u/ayan/MLDisk/medlda_files/sLDAtestlabels';

fp3 = fopen(trfilenameslda, 'w+');
fp4 = fopen(trclassesslda, 'w+');
fp5 = fopen(testfilenameslda, 'w+');
fp6 = fopen(testclassesslda, 'w+');

N1 = length(trdata.windex);
N2 = length(testdata.windex);

for i=1:max(N1,N2)
    if(i<=N1)
        str1 = convertToMedLDAformat(trdata.windex{i}', trdata.wcount{i}, trdata.classlabels(i));
        str11 = convertTosLDAformat(trdata.windex{i}', trdata.wcount{i});
        str12 = num2str(trdata.classlabels(i)-1);
        fprintf(fp1, '%s\n', str1);
        fprintf(fp3, '%s\n', str11); 
        fprintf(fp4, '%s\n', str12); 
    end
    if(i<=N2)
        str2 = convertToMedLDAformat(testdata.windex{i}', testdata.wcount{i}, testdata.classlabels(i));
        str21 = convertTosLDAformat(testdata.windex{i}', testdata.wcount{i});
        str22 = num2str(testdata.classlabels(i)-1);
        fprintf(fp2, '%s\n', str2);
        fprintf(fp5, '%s\n', str21); 
        fprintf(fp6, '%s\n', str22); 
    end
end

fclose(fp1);
fclose(fp2);
fclose(fp3);
fclose(fp4);
fclose(fp5);
fclose(fp6);

disp('inside getmedLDAperf');

estinfstr = './medlda estinf 50 12 0 100 1.0 random'
%%pause
%%system(estinfstr);

eststrslda = ['/lusr/u/ayan/Documents/DSLDA_SDM/slda/./slda est ' trfilenameslda ' ' trclassesslda ' settings.txt' ' 1.0' ' 20' ' random' ' /lusr/u/ayan/MLDisk/medlda_files/results']
%%system(eststrslda);
infstrslda = ['/lusr/u/ayan/Documents/DSLDA_SDM/slda/./slda inf ' testfilenameslda ' ' testclassesslda ' settings.txt' ' /lusr/u/ayan/MLDisk/medlda_files/results/final.model' ' /lusr/u/ayan/MLDisk/medlda_files/results'] 
%%system(infstrslda);
%%pause;

end
