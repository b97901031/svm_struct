/***********************************************************************/
/*                                                                     */
/*   svm_struct_api.c                                                  */
/*                                                                     */
/*   Definition of API for attaching implementing SVM learning of      */
/*   structures (e.g. parsing, multi-label classification, HMM)        */ 
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 03.07.04                                                    */
/*                                                                     */
/*   Copyright (c) 2004  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-c:ommercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/

#include <stdio.h>
//#include <string.h>
#include <stdbool.h>
//#include "svm_struct/svm_struct_common.h"
//#include "svm_struct_api.h"
#include "math.h"
#include "np_helper.h"

void        svm_struct_learn_api_init(int argc, char* argv[])
{
  /* Called in learning part before anything else is done to allow
     any initializations that might be necessary. */
}

void        svm_struct_learn_api_exit()
{
  /* Called in learning part at the very end to allow any clean-up
     that might be necessary. */
}

void        svm_struct_classify_api_init(int argc, char* argv[])
{
  /* Called in prediction part before anything else is done to allow
     any initializations that might be necessary. */
}

void        svm_struct_classify_api_exit()
{
  /* Called in prediction part at the very end to allow any clean-up
     that might be necessary. */
}

SAMPLE      read_struct_examples(char *trainfile, char *wordFeaturefile, char *predicateFeaturefile, char *structfile, STRUCT_LEARN_PARM *sparm)
{
  /* Reads struct examples and returns them in sample. The number of
     examples must be written into sample.n */
  /* fill in your code here */
  return read_struct_examples_helper(trainfile, wordFeaturefile, predicateFeaturefile, structfile, sparm);
}

void        init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm)
{
     printf("init_struct_model...\n");
  /* Initialize structmodel sm. The weight vector w does not need to be
     initialized, but you need to provide the maximum size of the
     feature space in sizePsi. This is the maximum number of different
     weights that can be learned. Later, the weight vector w will
     contain the learned weights for the model. */
  int feature_size = sample.examples[0].x.feature_size + sample.examples[0].x.word_feature_size + sample.examples[0].x.predicate_feature_size;
  int label_size = sample.examples[0].x.label_size;
  

  if (sparm->joint==2) {
    sm->sizePsi= label_size*(2*feature_size+2) + 2*pow(label_size,2); 
    /* replace by appropriate number of features */
  } else {
    sm->sizePsi= label_size*feature_size + 2*pow(label_size,2);
  }
  // sample.label_num*sample.feature_num: x-y
  // pow(sample.label,2): y-y
  // pow(sample.label,2): multi-predicates

}

CONSTSET    init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Initializes the optimization problem. Typically, you do not need
     to change this function, since you want to start with an empty
     set of constraints. However, if for example you have constraints
     that certain weights need to be positive, you might put that in
     here. The constraints are represented as lhs[i]*w >= rhs[i]. lhs
     is an array of feature vectors, rhs is an array of doubles. m is
     the number of constraints. The function returns the initial
     set of constraints. */
  //printf("init_struct_constraints\n");
  CONSTSET c;
  long     sizePsi=sm->sizePsi;
  long     i;
  WORD     words[2];

  if(1) { /* normal case: start with empty set of constraints */
    c.lhs=NULL;
    c.rhs=NULL;
    c.m=0;
  }
  else { /* add constraints so that all learned weights are
            positive. WARNING: Currently, they are positive only up to
            precision epsilon set by -e. */
    c.lhs=my_malloc(sizeof(DOC *)*sizePsi);
    c.rhs=my_malloc(sizeof(double)*sizePsi);
    for(i=0; i<sizePsi; i++) {
      words[0].wnum=i+1;
      words[0].weight=1.0;
      words[1].wnum=0;
      /* the following slackid is a hack. we will run into problems,
         if we have move than 1000000 slack sets (ie examples) */
      c.lhs[i]=create_example(i,0,1000000+i,1,create_svector(words,"",1.0));
      c.rhs[i]=0.0;
    }
  }
  return(c);
}

LABEL       classify_struct_example(PATTERN x, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label yhat for pattern x that scores the highest
     according to the linear evaluation function in sm, especially the
     weights sm.w. The returned label is taken as the prediction of sm
     for the pattern x. The weights correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. If the
     function cannot find a label, it shall return an empty label as
     recognized by the function empty_label(y). */

  //printf("classify_struct_example\n");

  int i;

  LABEL y;
  y.labels = (int*)malloc(sizeof(int)*x.length);
  y.arguments = (int*)malloc(sizeof(int)*x.length);
  y.length = x.length;
  y.label_size = x.label_size;


  /* insert your code for computing the predicted label y here */
  if (sparm->joint==0) { //use gold standard.
    initialization_helper(x,&y,sm,1);
    if (sparm->decode==0){ // beam.
      beam_search_helper(x,&y,y,sm,sparm,0);
    } else if (sparm->decode==1) { //gibbs.
      gibbs_sampling_helper(x,&y,y,sm,sparm,0);
    }
  } else if (sparm->joint==1) { //only classification 
    initialization_helper(x,&y,sm,0);
    if (sparm->decode==0){ // beam.
      beam_search_helper(x,&y,y,sm,sparm,0);
    } else if (sparm->decode==1) { //gibbs.
      gibbs_sampling_helper(x,&y,y,sm,sparm,0);
    }
  } else if (sparm->joint==2) { //identification+classification
    initialization_joint_helper(x,&y,sm,0);
    if (sparm->decode==0){ // beam.
      beam_search_joint_helper(x,&y,y,sm,sparm,0);
    } else if (sparm->decode==1) { //gibbs.
      gibbs_sampling_joint_helper(x,&y,y,sm,sparm,0);
    } 
  }

  return(y);
}

LABEL       find_most_violated_constraint_slackrescaling(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label ybar for pattern x that that is responsible for
     the most violated constraint for the slack rescaling
     formulation. For linear slack variables, this is that label ybar
     that maximizes

            argmax_{ybar} loss(y,ybar)*(1-psi(x,y)+psi(x,ybar)) 

     Note that ybar may be equal to y (i.e. the max is 0), which is
     different from the algorithms described in
     [Tschantaridis/05]. Note that this argmax has to take into
     account the scoring function in sm, especially the weights sm.w,
     as well as the loss function, and whether linear or quadratic
     slacks are used. The weights in sm.w correspond to the features
     defined by psi() and range from index 1 to index
     sm->sizePsi. Most simple is the case of the zero/one loss
     function. For the zero/one loss, this function should return the
     highest scoring label ybar (which may be equal to the correct
     label y), or the second highest scoring label ybar, if
     Psi(x,ybar)>Psi(x,y)-1. If the function cannot find a label, it
     shall return an empty label as recognized by the function
     empty_label(y). */
  //printf("find_most_violated_constraint_slackrescaling\n");

  LABEL ybar;

  /* insert your code for computing the label ybar here */

  return(ybar);
}

LABEL       find_most_violated_constraint_marginrescaling(PATTERN x, LABEL y, 
						     STRUCTMODEL *sm, 
						     STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label ybar for pattern x that that is responsible for
     the most violated constraint for the margin rescaling
     formulation. For linear slack variables, this is that label ybar
     that maximizes

            argmax_{ybar} loss(y,ybar)+psi(x,ybar)

     Note that ybar may be equal to y (i.e. the max is 0), which is
     different from the algorithms described in
     [Tschantaridis/05]. Note that this argmax has to take into
     account the scoring function in sm, especially the weights sm.w,
     as well as the loss function, and whether linear or quadratic
     slacks are used. The weights in sm.w correspond to the features
     defined by psi() and range from index 1 to index
     sm->sizePsi. Most simple is the case of the zero/one loss
     function. For the zero/one loss, this function should return the
     highest scoring label ybar (which may be equal to the correct
     label y), or the second highest scoring label ybar, if
     Psi(x,ybar)>Psi(x,y)-1. If the function cannot find a label, it
     shall return an empty label as recognized by the function
     empty_label(y). */
  //printf("find_most_violated_constraint_marginrescaling\n");

  LABEL ybar;
  ybar.labels = (int*)malloc(sizeof(int)*x.length);
  ybar.arguments = (int*)malloc(sizeof(int)*x.length);
  ybar.length = x.length;
  ybar.label_size = x.label_size;  

  /* insert your code for computing the label ybar here */
  int i;
  if (sparm->joint==0) { //use gold standard.
    initialization_helper(x,&ybar,sm,1);
    if (sparm->decode==0){ // beam.
      beam_search_helper(x,&ybar,y,sm,sparm,1);
    } else if (sparm->decode==1) { //gibbs.
      gibbs_sampling_helper(x,&ybar,y,sm,sparm,1);
    }
  } else if (sparm->joint==1) { //only classification 
    initialization_helper(x,&ybar,sm,0);
    if (sparm->decode==0){ // beam.
      beam_search_helper(x,&ybar,y,sm,sparm,1);
    } else if (sparm->decode==1) { //gibbs.
      gibbs_sampling_helper(x,&ybar,y,sm,sparm,1);
    }
  } else if (sparm->joint==2) { //identification+classification
    initialization_joint_helper(x,&ybar,sm,0);
    if (sparm->decode==0){ // beam.
      beam_search_joint_helper(x,&ybar,y,sm,sparm,1);
    } else if (sparm->decode==1) { //gibbs.
      gibbs_sampling_joint_helper(x,&ybar,y,sm,sparm,1);
    }
  } 
  return(ybar);
}

int         empty_label(LABEL y)
{
  /* Returns true, if y is an empty label. An empty label might be
     returned by find_most_violated_constraint_???(x, y, sm) if there
     is no incorrect label that can be found for x, or if it is unable
     to label x at all */
  //printf("empty_label\n");
  return(0);
}

SVECTOR     *psi(PATTERN x, LABEL y, STRUCTMODEL *sm,
		 STRUCT_LEARN_PARM *sparm)
{
  /* Returns a feature vector describing the match between pattern x
     and label y. The feature vector is returned as a list of
     SVECTOR's. Each SVECTOR is in a sparse representation of pairs
     <featurenumber:featurevalue>, where the last pair has
     featurenumber 0 as a terminator. Featurenumbers start with 1 and
     end with sizePsi. Featuresnumbers that are not specified default
     to value 0. As mentioned before, psi() actually returns a list of
     SVECTOR's. Each SVECTOR has a field 'factor' and 'next'. 'next'
     specifies the next element in the list, terminated by a NULL
     pointer. The list can be though of as a linear combination of
     vectors, where each vector is weighted by its 'factor'. This
     linear combination of feature vectors is multiplied with the
     learned (kernelized) weight vector to score label y for pattern
     x. Without kernels, there will be one weight in sm.w for each
     feature. Note that psi has to match
     find_most_violated_constraint_???(x, y, sm) and vice versa. In
     particular, find_most_violated_constraint_???(x, y, sm) finds
     that ybar!=y that maximizes psi(x,ybar,sm)*sm.w (where * is the
     inner vector product) and the appropriate function of the
     loss + margin/slack rescaling method. See that paper for details. */
 
  //printf("psi\n");

  SVECTOR *fvec=NULL;
  WORD *words;
  words = (WORD*)my_malloc(sizeof(WORD)*(sm->sizePsi+1));  

  int i;
  int k=0;
  double *feature_vector = feature_vector_helper(x,y,sm,sparm);
  for (i=1; i<sm->sizePsi+1; i++){
    words[k].wnum = i;
    words[k].weight = feature_vector[k];
    k++;
  }
  words[k].wnum=0;
  words[k].weight=0.0;

  fvec = create_svector(words,"",1);
  free(words);
  free(feature_vector);  

  return(fvec);
}

double      loss(LABEL y, LABEL ybar, STRUCT_LEARN_PARM *sparm)
{
  /* loss for correct label y and predicted label ybar. The loss for
     y==ybar has to be zero. sparm->loss_function is set with the -l option. */
  int i, count, count_id, count_cl;

  //printf("loss\n");


//  if(sparm->loss_function == 0) { /* type 0 loss: 0/1 loss */
//    return 0;                      /* return 0, if y==ybar. return 1 else */
//  }
//  else {
    /* Put your code for different loss functions here. But then
       find_most_violated_constraint_???(x, y, sm) has to return the
       highest scoring label with the largest loss. */
    if (sparm->joint==2) {
      count_id = 0;
      count_cl = 0;
      for (i=0; i<y.length; i++) {
        if (y.labels[i]==ybar.labels[i]){
          count_cl++;
        }
        if (y.arguments[i]==ybar.arguments[i]) {
          count_id++;
        }
      }
      return 1.0-((double)count_id/(double)y.length+(double)count_cl/(double)y.length)/2.0;
    }
    else {
      count = 0;
      for (i=0; i<y.length; i++) {
        if (y.labels[i]==ybar.labels[i]) {
          count++;
        }
      }
      //printf("count= %d, length= %d\n",count,y.length); 
      //printf("loss= %lf\n",1.0-(double)count/(double)y.length);
      return 1.0-(double)count/(double)y.length;
    }
//  }
}

int         finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
			       STRUCT_LEARN_PARM *sparm)
{
  /* This function is called just before the end of each cutting plane iteration. ceps is the amount by which the most violated constraint found in the current iteration was violated. cached_constraint is true if the added constraint was constructed from the cache. If the return value is FALSE, then the algorithm is allowed to terminate. If it is TRUE, the algorithm will keep iterating even if the desired precision sparm->epsilon is already reached. */
  return(0);
}

void        print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm)
{
  /* This function is called after training and allows final touches to
     the model sm. But primarly it allows computing and printing any
     kind of statistic (e.g. training error) you might want. */
}

void        print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm, 
				       STRUCT_TEST_STATS *teststats)
{
  /* This function is called after making all test predictions in
     svm_struct_classify and allows computing and printing any kind of
     evaluation (e.g. precision/recall) you might want. You can use
     the function eval_prediction to accumulate the necessary
     statistics for each prediction. */
}

void        eval_prediction(long exnum, EXAMPLE ex, LABEL ypred, 
			    STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, 
			    STRUCT_TEST_STATS *teststats)
{
  /* This function allows you to accumlate statistic for how well the
     predicition matches the labeled example. It is called from
     svm_struct_classify. See also the function
     print_struct_testing_stats. */
  if(exnum == 0) { /* this is the first time the function is
		      called. So initialize the teststats */
  }
}

void        write_struct_model(char *file, STRUCTMODEL *sm, 
			       STRUCT_LEARN_PARM *sparm)
{
  /* Writes structural model sm to file file. */
/*
  Writes the learned weight vector sm->w to file after training. 
*/

  //printf("write_struct_model\n");

  FILE *modelfl;
  int i;
  
  modelfl = fopen(file,"w");
  if (modelfl==NULL) {
    printf("Cannot open model file %s for output!", file);
	exit(1);
  }
  
  /* write model information */
  fprintf(modelfl, "# size of psi: %ld\n", sm->sizePsi);
  fprintf(modelfl, "# joint: %d\n", sparm->joint);
  fprintf(modelfl, "# decode: %d\n", sparm->decode);

  for (i=1;i<sm->sizePsi+1;i++) {
    fprintf(modelfl, "%d:%.16g\n", i, sm->w[i]);
  }
  fclose(modelfl);
}

STRUCTMODEL read_struct_model(char *file, STRUCT_LEARN_PARM *sparm)
{
  /* Reads structural model sm from file file. This function is used
     only in the prediction module, not in the learning module. */

  STRUCTMODEL sm;

  FILE *modelfl;
  int sizePsi,i, fnum, joint, decode;
  double fweight;
  //char line[1000];
  
  modelfl = fopen(file,"r");
  if (modelfl==NULL) {
    printf("Cannot open model file %s for input!", file);
	exit(1);
  }

  if (fscanf(modelfl, "# size of psi: %d\n", &sizePsi)!=1) {
    printf("Cannot read model: lacking sizePsi information\n");
    exit(1);
  }

  if (fscanf(modelfl, "# joint: %d\n", &joint)!=1) {
    printf("Cannot read model: lacking joint information\n");
    exit(1);
  }

  if (fscanf(modelfl, "# decode: %d\n", &decode)!=1) {
    printf("Cannot read model: lacking decode information\n");
    exit(1);
  }

  sparm->joint = joint;
  sparm->decode = decode;

  sm.sizePsi = sizePsi;
  sm.w = (double*)malloc((sizePsi+1)*sizeof(double));

  for (i=0;i<sizePsi+1;i++) {
    sm.w[i] = 0.0;
  }
  
  while (!feof(modelfl)) {
    fscanf(modelfl, "%d:%lf", &fnum, &fweight);
	sm.w[fnum] = fweight;
  }

  fclose(modelfl);
  return(sm);

}

void        write_label(FILE *fp, LABEL y, LABEL ybar, PATTERN x)
{
  /* Writes label y to file handle fp. */
  int i;
  fprintf(fp, "#%s\n", x.filename);
  for (i=0; i<y.length; i++) {
    fprintf(fp, "%d %d\n", y.labels[i],ybar.labels[i]);
  }
} 

void write_confusion(FILE *fp, int **confusion, int size) {
  int i, j;
  fprintf(fp,"\n#confusion matrix:\n");
  for (i=0; i<size; i++) {
    for (j=0; j<size; j++) {
      fprintf(fp, "%d ", confusion[i][j]);
    }
    fprintf(fp, "\n");
  }
}

void        free_pattern(PATTERN x) {
  /* Frees the memory of x. */
  int i;
  for (i=0; i<x.length; i++){
    free(x.predicates[i]);
    free(x.words[i]);
    free(x.features[i]);
  }
  for (i=0; i<x.predicate_num; i++) {
    free(x.predicates[i]);
  }
  for (i=0; i<x.word_num; i++) {
    free(x.words[i]);
  }

  for (i=0; i<x.identical_length; i++) {
    free(x.structures[i]);
    free(x.structToIndex[i]);
  }

  free(x.structures);
  free(x.predicates);
  free(x.words);
  free(x.features);
  free(x.structToIndexSize);
  free(x.structToIndex);
}

void        free_label(LABEL y) {
  /* Frees the memory of y. */
  free(y.labels);
  free(y.arguments);
}

void        free_struct_model(STRUCTMODEL sm) 
{
  /* Frees the memory of model. */
  /* if(sm.w) free(sm.w); */ /* this is free'd in free_model */
  if(sm.svm_model) free_model(sm.svm_model,1);
  /* add free calls for user defined data here */
}

void        free_struct_sample(SAMPLE s)
{
  /* Frees the memory of sample s. */
  int i;
  for(i=0;i<s.n;i++) { 
    free_pattern(s.examples[i].x);
    free_label(s.examples[i].y);
  }
  free(s.examples);
}

void        print_struct_help()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_learn. */
  printf("         --* string  -> custom parameters that can be adapted for struct\n");
  printf("                        learning. The * can be replaced by any character\n");
  printf("                        and there can be multiple options starting with --.\n");
}

void         parse_struct_parameters(STRUCT_LEARN_PARM *sparm)
{
  /* Parses the command line parameters that start with -- */
  int i;

  for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
    switch ((sparm->custom_argv[i])[2]) 
      { 
      case 'a': i++; /* strcpy(learn_parm->alphafile,argv[i]); */ break;
      case 'e': i++; /* sparm->epsilon=atof(sparm->custom_argv[i]); */ break;
      case 'k': i++; /* sparm->newconstretrain=atol(sparm->custom_argv[i]); */ break;
      default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
	       exit(0);
      }
  }
}

void        print_struct_help_classify()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_classify. */
  printf("         --* string -> custom parameters that can be adapted for struct\n");
  printf("                       learning. The * can be replaced by any character\n");
  printf("                       and there can be multiple options starting with --.\n");
}

void         parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm)
{
  /* Parses the command line parameters that start with -- for the
     classification module */
  int i;

  for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
    switch ((sparm->custom_argv[i])[2]) 
      { 
      /* case 'x': i++; strcpy(xvalue,sparm->custom_argv[i]); break; */
      default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
	       exit(0);
      }
  }
}

