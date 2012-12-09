import sys, os
from corpus import *
import hdp
import cPickle
import utils
from utils import *
import random, time
from numpy import cumsum, sum, shape
from itertools import izip
from optparse import OptionParser
from glob import glob
import scipy.special as sp
sys.path.insert(0, "/lusr/u/ayan/Documents/DSLDA_SDM/DSLDA/liblinear-1.92/python/")
from liblinearutil import *
np = hdp.np

def parse_args():
  parser = OptionParser()
  os.system('clear')
  parser.set_defaults(T=6, K=3, D=27, W=13412, eta=0.1, alpha=6.0, gamma=5.0,
                      max_time=5, max_iter=2, var_converge=0.0001, random_seed=999931111, 
                      corpus_name=None, data_path="trdataSLDA.txt", test_data_path="trdataSLDA.txt", 
                      test_data_path_in_folds=None, directory=None)

  parser.add_option("--T", type="int", dest="T",
                    help="top level truncation [300]")
  parser.add_option("--K", type="int", dest="K",
                    help="second level truncation [20]")
  parser.add_option("--D", type="int", dest="D",
                    help="number of documents [-1]")
  parser.add_option("--W", type="int", dest="W",
                    help="size of vocabulary [-1]")
  parser.add_option("--eta", type="float", dest="eta",
                    help="the topic Dirichlet [0.01]")
  parser.add_option("--alpha", type="float", dest="alpha",
                    help="alpha value [1.0]")
  parser.add_option("--gamma", type="float", dest="gamma",
                    help="gamma value [1.0]")
  parser.add_option("--max_time", type="int", dest="max_time",
                    help="max time to run training in seconds [100]")
  parser.add_option("--max_iter", type="int", dest="max_iter",
                    help="max iteration to run training [-1]")
  parser.add_option("--var_converge", type="float", dest="var_converge",
                    help="relative change on doc lower bound [0.0001]")
  parser.add_option("--random_seed", type="int", dest="random_seed",
                    help="the random seed [999931111]")
  parser.add_option("--corpus_name", type="string", dest="corpus_name",
                    help="the corpus name: nature, nyt or wiki [None]")
  parser.add_option("--data_path", type="string", dest="data_path",
                    help="training data path or pattern [None]")
  parser.add_option("--test_data_path", type="string", dest="test_data_path",
                    help="testing data path [None]")
  parser.add_option("--test_data_path_in_folds", type="string",
                    dest="test_data_path_in_folds",
                    help="testing data prefix for different folds [None]")
  parser.add_option("--directory", type="string", dest="directory",
                    help="output directory [None]")

  (options, args) = parser.parse_args()
  return options 

def run_hdp():
  options = parse_args()
  random.seed(options.random_seed)

  # Read the training data.
  c_train_filename = options.data_path
  c_train  = read_data(c_train_filename)
  trlabels = utils.read_1D_array('trdataSLDALabels.txt');
  trsz     = utils.read_1D_array('trdatasz.txt');
  trsz = trsz[0];

  # Read the test data.
  test_data_path = options.test_data_path
  c_test = read_data(test_data_path)
  c_test_word_count = sum([doc.total for doc in c_test.docs])
  testlabels = utils.read_1D_array('trdataSLDALabels.txt');
  testsz     = utils.read_1D_array('trdatasz.txt');
  testsz = testsz[0];

  f = file("localsticks_a.txt", "w") 
  f.close()
  f = file("localsticks_b.txt", "w") 
  f.close()

  result_directory = "%s/corpus-%s" % (options.directory, options.corpus_name)
  print "creating directory %s" % result_directory
  if not os.path.isdir(result_directory):
    os.makedirs(result_directory)

  options_file = file("%s/options.dat" % result_directory, "w")
  for opt, value in options.__dict__.items():
    options_file.write(str(opt) + " " + str(value) + "\n")
  options_file.close()

  print "creating hdp instance for training."
  bhdp_hp = hdp.hdp_hyperparameter(options.alpha, options.alpha, options.gamma, options.gamma, True)
  bhdp = hdp.hdp(options.T, options.K, options.D, options.W, options.eta, trsz, bhdp_hp)

  print "setting up counters and log files."
  iter = 0 
  total_time = 0.0
  total_doc_count = 0

  likelihood = 0.0
  old_likelihood = 0.0
  converge = 1.0

  while options.max_iter == -1 or iter < options.max_iter:
    t0 = time.clock()
    # Run one step iteration.
    print "EM iteration starts!" 
    likelihood, trmodel = bhdp.em(c_train, trlabels, options.var_converge, trsz, fresh=(iter==0))
    if iter > 0:
        converge = (likelihood - old_likelihood)/abs(old_likelihood)
    old_likelihood = likelihood
    print "iter = %d, likelihood = %f, converge = %f" % (iter, likelihood, converge)
    if converge < 0:
        print "warning, likelihood is decreasing!"
    total_time += time.clock() - t0
    iter += 1  # increase the iter counter
    total_doc_count += options.D

  raw_input()
  # prediction on the fixed test data
  print "\tworking on fixed test data."
  print "creating hdp instance for test."
  bhdp_hp_test = hdp.hdp_hyperparameter(options.alpha, options.alpha, options.gamma, options.gamma, False)
  bhdp_test = hdp.hdp(options.T, options.K, options.D, options.W, options.eta, testsz, bhdp_hp_test)
  # prediction on test data 
  Elogbeta = dirichlet_expectation(bhdp.m_beta) # the topics
  Elogsticks_1st = expect_log_sticks(bhdp.m_var_sticks) # global sticks
  m_var_zeta = np.zeros((testsz, options.K))
  docnum = 0
  for doc in c_test.docs:
     temp, m_var_zeta = bhdp_test.doc_inference(doc, docnum, Elogbeta, Elogsticks_1st, options.var_converge, m_var_zeta)
     likelihood += temp
     docnum +=1  

  x = utils.converttoliblinear(m_var_zeta)
  y = testlabels
  p_label, p_acc, p_val = predict(y, x, trmodel)
  ACC, MSE, SCC = evaluations(y, p_label)
  print "test accuracy: %f" %ACC

if __name__ == '__main__':
  run_hdp()
