function [accval, model] = mainfile (MAXCOUNT, MaxFun, MAXESTEPITER, MAXMSTEPITER, N, V, K, Y, option)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file for simulating synthetic data and checking the performance of the model on the same
% call variational_EM() directly with appropriate arguments for any other dataset  
%%%%%%%%%%%
% Input:
% MAXCOUNT: maximum number of EM iteration
% MAxFun: maximum number of function evaluations for fmincon
% MAXESTEPITER: maximum number of E-step iteration
% MAXMSTEPITER: maximum number of M-step iteration
% N: number of instances
% V: vocabulary size
% K: total number of topics (visible+latent) -- in synthetic data, K/2 visible and latent topics are maintained
% Y: number of classes
% option: 1,2,3,4 -- use 4 for the model we have; 1 for unsupervised LDA, 2 for labeled LDA, 3 for Med-LDA
% phase: 1 for training, 0 for testing 
%%%%%%%%%%%
% Output:
% accval: accuracy of classification
% model: model with all the parameter (variational and model) values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mainfile (2, 2, 2, 2, 20, 40, 5, 4, 4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long;
warning('off');
clc;

% these C++ files should be compiled first
mex sum_phi.cpp;
mex update_beta_cpp.cpp;
mex update_phi_cpp.cpp;

addpath(genpath('/lusr/u/ayan/Documents/aaaa/liblinear-v0.1/'));
addpath(genpath('/lusr/u/ayan/Documents/aaaa/siftpp/'));

p    = 0.9; %(should be > 0) multiplier for Dirichlet Distribution; change this to vary initialization 
trdata = synthdata(N, V, K, Y, option);
trdata.windex{1}
trdata.wcount{1}
sasx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data is a structure containing the following fields:
% annotations: N*k1 binary matrix with 1 indicating presence of an annotation in a document, N indicates number of documents and k1 indicates maximum number of visible topics/annotations
% w: an N*1 cell array with w{n} containing the indicies of the words in the vocabulary appearing in nth document.
% classlabels: N*1 matrix with each element indicating the class label of corresponding document
% Y: number of classes (at least one document from each of the classes should be present in the training data)
% k1: number of visible topics/annotations
% k2: number of latent topics -- note that K = k1+K2
% nwordspdoc: N*1 matrix with each element indicating the number of words in corresponding document
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
[trmodel, traccval] = variational_EM(trdata, MAXCOUNT, MAXESTEPITER, MAXMSTEPITER, MaxFun, p, K, option, 1, []); 
[testmodel, testaccval] = InferenceOnTest(trmodel, trdata, MAXESTEPITER, option, p); 

traccval
testaccval

end

