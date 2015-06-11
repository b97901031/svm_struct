#ifdef __cplusplus
extern "C" {
#endif
#include "svm_struct_api_types.h"
#include "./svm_struct/svm_struct_common.h"
#include "./svm_light/svm_common.h"
//#include <vector>
//#include <string>

//vector<string> split(const string &s, char delim);

SAMPLE read_struct_examples_helper(char *labelfile, char *wordFeaturefile, char *predicateFeaturefile, char *structfile, STRUCT_LEARN_PARM *sparm);

double inner_product(PATTERN x, STRUCTMODEL *sm, int start, int feature_size);

int find_element(char **list, int list_length, char *target);

void initialization_helper(PATTERN x, LABEL *y, STRUCTMODEL *sm, int gold);

void initialization_joint_helper(PATTERN x, LABEL *y, STRUCTMODEL *sm, int gold);

void gibbs_sampling_helper(PATTERN x, LABEL *y, LABEL ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int useLoss);

void gibbs_sampling_joint_helper(PATTERN x, LABEL *y, LABEL ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int useLoss);

void beam_search_helper(PATTERN x, LABEL *y, LABEL ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int useLoss);

void beam_search_joint_helper(PATTERN x, LABEL *y, LABEL ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int useLoss);

double* feature_vector_helper(PATTERN x, LABEL y, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm);

int converge(int* y1, int* y2, int size, double constraint);


#ifdef __cplusplus
}
#endif
