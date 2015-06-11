#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream> 
#include <algorithm>
#include "np_helper.h"
#include "math.h"
#include <stdbool.h>
#include <stdio.h> 
#include <string>
//#include "./svm_struct/svm_struct_common.h"

extern "C" {
//#include "svm_struct_api_types.h"
#include "./svm_light/svm_common.h"
#include "./svm_struct/svm_struct_common.h"
}

#define DEBUG 0

using namespace std;
/*
vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}
*/

SAMPLE read_struct_examples_helper(char *labelfile, char *wordFeaturefile, char *predicateFeaturefile, char *structfile, STRUCT_LEARN_PARM *sparm) {

  // file - word - fea
  vector< vector< vector<double> > > features;
  vector< vector<double> > feature;
  vector<double> fea;

  vector< vector< vector<double> > > wordFeatures;
  vector< vector< vector<double> > > predicateFeatures;

  /* change -> identical_index to index. should read filemap to get the ind for indexMaps */ 
  vector< map<int, vector<int> > > indexMaps;
  map<int, vector<int> > indexMap;

  /* word id to identical_index */ 
  vector< map<string,int> > identicalIndexMaps;
  map<string,int> identicalIndexMap;

  /* word id, predicate, label*/
  vector< vector<int> > labels;
  vector<int> label;

  vector< vector<string> > identicalWords;
  vector< vector<string> > words;
  vector<string> ws;

  vector< vector<string> > identicalPredicates;
  vector< vector<string> > predicates;
  vector<string> predicate;

  vector<int> length;
  vector<int> identicalLength;
  map<string,int> filemap; 

  vector<string> filenames;

  int i,j,k;
  int key, count, identical_count, fea_count, word_fea_count, predicate_fea_count, prev_fea_count, lab;
  int filecount=0;
  int label_size=0;
  double value;
  bool inVector, flag; 
  string id, filename, pre;

  printf("Read label file...\n");


  prev_fea_count=-1;
  ifstream fin(labelfile);
  string line;
  while ( getline(fin, line) ) {
    vector<string> elems;
    stringstream ss(line);
    string item;
    while (getline(ss, item, ' ')) {
        elems.push_back(item);
    }

    if (elems.size()>3){
      id = elems[0];
      pre = elems[1];
      lab = stoi(elems[2]);
      label_size = max(label_size,lab+1);

      // if word id already in list.
      if (identicalIndexMap.find(id)==identicalIndexMap.end()) {
        identicalIndexMap[id]=identical_count;
        indexMap[identical_count]=vector<int> (1,count);
        identical_count++;
      } else {
        indexMap[identicalIndexMap[id]].push_back(count);
      }
      count++;

      // start to read features 
      fea_count = elems.size()-3;
      if (!fea.empty()) {printf("feature vector suffers from data corruption.\n");}
      for (i=3; i<elems.size(); i++){
        vector<string> fx;
        stringstream ss2(elems[i]);
        while(getline(ss2, item, ':')){
          fx.push_back(item);
        }
        fea.push_back(stod(fx[1]));
      }
  
      if (prev_fea_count!=-1 and fea_count!=prev_fea_count) printf("feature length error!\n");
      prev_fea_count = fea_count;

      feature.push_back(fea);
      label.push_back(lab);
      ws.push_back(id);
      predicate.push_back(pre);
      fea.clear();

    } else if (elems.size()==2) {
      // start of a sentence
      if (elems[0].compare("#")==0) {
        count = 0;
        identical_count = 0;
        filename = elems[1];
        filemap[filename]=filecount;
        filecount+=1;
        filenames.push_back(filename);
      }
    } else if (elems.size()==1) {
      // end of a sentence   
      features.push_back(feature);
      for (i=0; i<feature.size(); i++) {feature[i].clear();}
      feature.clear();
      fea.clear();

      indexMaps.push_back(indexMap);
      indexMap.clear();

      identicalIndexMaps.push_back(identicalIndexMap);
      identicalIndexMap.clear();

      labels.push_back(label);
      label.clear();

      predicates.push_back(predicate);
      predicate.clear();

      words.push_back(ws);
      ws.clear();

      length.push_back(count);
      identicalLength.push_back(identical_count);
    }
  }

  printf("Read word feature file...\n");

  ifstream finw(wordFeaturefile);
  while ( getline(finw, line) ) {  
    vector<string> elems;
    stringstream ss(line);
    string item;
    while (getline(ss, item, ' ')) {elems.push_back(item);}

    if (elems[0].compare("#")==0) {
      if (elems.size()>1) {
        // start of a sentence 
      } else {
        // end of a sentence
        identicalWords.push_back(ws);
        ws.clear();
        wordFeatures.push_back(feature);
      }
    } else {
      id = elems[0];
      ws.push_back(id);
      word_fea_count = elems.size()-1;
      if (!fea.empty()) {printf("feature vector suffers from data corruption.\n");}
      for (i=1; i<elems.size(); i++){
        vector<string> fx;
        stringstream ss2(elems[i]);
        while(getline(ss2, item, ':')){
          fx.push_back(item);
        }
        fea.push_back(stod(fx[1]));
      }
      feature.push_back(fea);
      fea.clear();
    }
  }

  printf("Read predicate file...\n");

  ifstream finp(predicateFeaturefile);
  while ( getline(finp, line) ) {  
    vector<string> elems;
    stringstream ss(line);
    string item;
    while (getline(ss, item, ' ')) {
        elems.push_back(item);
    }

    if (elems[0].compare("#")==0) {
      if (elems.size()>1) {
        // start of a sentence 
      } else {
        // end of a sentence
        identicalPredicates.push_back(predicate);
        predicate.clear();
        predicateFeatures.push_back(feature);
      }
    } else {
      id = elems[0];
      predicate.push_back(id);
      predicate_fea_count = elems.size()-1;
      if (!fea.empty()) {printf("feature vector suffers from data corruption.\n");}
      for (i=1; i<elems.size(); i++){
        vector<string> fx;
        stringstream ss2(elems[i]);
        while(getline(ss2, item, ':')){
          fx.push_back(item);
        }
        fea.push_back(stod(fx[1]));
      }
      feature.push_back(fea);
      fea.clear();
    }
  }

  printf("Read structure file...\n");

  vector< vector< vector<int> > > structure;
  
  vector< vector<int> > mm;
  
  int index, ind1, ind2, rel;
  string w1, w2; 

  cout<<"train: "<<labelfile<<endl;
  cout<<"xxx: "<<wordFeaturefile<<" "<<predicateFeaturefile<<endl;
  cout<<"struct: "<<structfile<<endl;
 
  ifstream fins(structfile);
  while ( getline(fins, line) ) {
    vector<string> elems;
    stringstream ss(line);
    string item;
    while ( getline(ss, item, ' ') ) {
        elems.push_back(item);
    }
    
    if (elems[0].compare("#")==0) {  
      if (elems.size()>1) {
        // start of a sentence
        filename = elems[1];
        index = filemap[filename];
        vector< vector<int> > mm (identicalLength[index], vector<int>(identicalLength[index],0));
      } else {
        // end of a sentence
        structure.push_back(mm);
        for (i=0; i<mm.size(); i++) {mm[i].clear();}
        mm.clear();
      }
    } else {
      w1 = elems[0];
      w2 = elems[1];
      ind1 = identicalIndexMaps[index][w1];
      ind2 = identicalIndexMaps[index][w2];
      mm[ind1][ind2]=1;
      mm[ind2][ind1]=-1;
    }
  }
 
  printf("Set examples...\n");

  SAMPLE sample;
  sample.n = features.size();
  sample.examples = (EXAMPLE*)malloc(sizeof(EXAMPLE)*sample.n);

  if (features.size()!=labels.size() || predicates.size()!=features.size()){
    printf("Length mismatch!\nfeatures:%d, labels:%d, predicates:%d ",(int)features.size(), (int)labels.size(), (int)predicates.size());
  }

  for (i =0; i<sample.n; i++) {
    // x part
    int node_size = length[i];
    int identical_node_size = identicalLength[i];
    
    sample.examples[i].x.label_size = label_size;
    sample.examples[i].y.label_size = label_size;

    char *cstr = new char[filenames[i].length()+1];
    strcpy(cstr, filenames[i].c_str());
    sample.examples[i].x.filename = cstr;
    //delete [] cstr;
 
    /* features */
    sample.examples[i].x.feature_size = fea_count;
    sample.examples[i].x.features= (double**)malloc(sizeof(double*)*node_size);
    for (j=0; j<node_size; j++){
      sample.examples[i].x.features[j]=(double*)malloc(sizeof(double)*fea_count);
      for (k=0; k<features[i][j].size(); k++) {
        sample.examples[i].x.features[j][k] = features[i][j][k];
      }
    }

    printf("a\n");
    /* wordFeatures */
    int word_num = identicalWords[i].size();
    sample.examples[i].x.word_num = word_num;
    sample.examples[i].x.word_feature_size = word_fea_count;
    //sample.examples[i].x.words = (char**)malloc(sizeof(char*)*node_size);
    sample.examples[i].x.identicalWords = (char**)malloc(sizeof(char*)*word_num);
    sample.examples[i].x.wordFeatures = (double**)malloc(sizeof(double*)*word_num); 
    for (j=0; j<node_size; j++) {
      char *cstr = new char[identicalWords[i][j].length()+1];
      strcpy(cstr, identicalWords[i][j].c_str());
      sample.examples[i].x.identicalWords[j] = cstr;
      sample.examples[i].x.wordFeatures[j] = (double*)malloc(sizeof(double)*word_fea_count);
      for (k=0; k<word_fea_count; k++) {
        sample.examples[i].x.wordFeatures[j][k] = wordFeatures[i][j][k]; 
      }
    }

    printf("b\n");
    /* predicateFeatures */
    int predicate_num =  identicalPredicates[i].size(); //predicateAppears[i].size();
    sample.examples[i].x.predicate_num = predicate_num;
    sample.examples[i].x.predicate_feature_size = predicate_fea_count;
    //sample.examples[i].x.predicate = (char**)malloc(sizeof(char*)*node_size);
    sample.examples[i].x.identicalPredicates = (char**)malloc(sizeof(char*)*predicate_num);
    sample.examples[i].x.predicateFeatures = (double**)malloc(sizeof(double*)*predicate_num); 
    for (j=0; j<node_size; j++) {
      char *cstr = new char[identicalPredicates[i][j].length()+1];
      strcpy(cstr, identicalPredicates[i][j].c_str());
      sample.examples[i].x.identicalPredicates[j] = cstr;
      sample.examples[i].x.predicateFeatures[j] = (double*)malloc(sizeof(double)*predicate_fea_count);
      for (k=0; k<predicate_fea_count; k++) {
        sample.examples[i].x.predicateFeatures[j][k] = predicateFeatures[i][j][k]; 
      }
    }


    printf("c\n");
    //flag=false;
    sample.examples[i].x.length = node_size;
    sample.examples[i].x.identical_length = identical_node_size;
    sample.examples[i].x.structures = (int**)malloc(sizeof(int*)*identical_node_size);
    for (j=0; j<identical_node_size; j++) {
      flag=false;
      sample.examples[i].x.structures[j]=(int*)malloc(sizeof(int)*identical_node_size);
      for (k=0; k<identical_node_size; k++) {
        sample.examples[i].x.structures[j][k]=structure[i][j][k];
      }
    }

    printf("d\n");
    sample.examples[i].x.structToIndexSize = (int*)malloc(sizeof(int)*identical_node_size);
    sample.examples[i].x.structToIndex = (int**)malloc(sizeof(int*)*identical_node_size);
    for (j=0; j<identical_node_size; j++) {
      sample.examples[i].x.structToIndex[j] = (int*)malloc(sizeof(int)*indexMaps[i][j].size());
      sample.examples[i].x.structToIndexSize[j] = indexMaps[i][j].size();
      for (k=0; k<indexMaps[i][j].size(); k++) {
        sample.examples[i].x.structToIndex[j][k]=indexMaps[i][j][k];
      }
    }
    //}

    printf("e\n");
    // y part
    sample.examples[i].y.length = node_size;
    sample.examples[i].y.labels = (int*)malloc(sizeof(int)*node_size);
    sample.examples[i].y.arguments = (int*) malloc(sizeof(int)*node_size);
    sample.examples[i].x.predicates = (char**)malloc(sizeof(char*)*node_size);
    sample.examples[i].x.words = (char**)malloc(sizeof(char*)*node_size);
    for (j=0; j<node_size; j++) {
      char *cstr = new char[predicates[i][j].length()+1];
      strcpy(cstr, predicates[i][j].c_str());
      char *dstr = new char[words[i][j].length()+1];
      strcpy(dstr, words[i][j].c_str());

      sample.examples[i].x.words[j] = dstr;
      sample.examples[i].x.predicates[j] = cstr; 
      sample.examples[i].y.labels[j] = labels[i][j];
      if (labels[i][j]==0) {
        sample.examples[i].y.arguments[j]=0;
      } else {
        sample.examples[i].y.arguments[j]=1;
      }
    }
  }

  return sample;
}

double inner_product(PATTERN x, STRUCTMODEL *sm, int start, int feature_size) {
  double score = 0.0;
  for (int i=0; i<x.length; i++) {
    for (int j=start; j<start+feature_size; j++) {
      score += sm->w[j]*x.features[i][j]; // check dummy = head or end?
    }
  }
  return score;
}

void initialization_helper(PATTERN x, LABEL *y, STRUCTMODEL *sm, int gold){

  //cout<<"initialization_helper"<<endl;

  int max_label, index, i, j;
  double score, max_score;
  max_score = -1*pow(10,6);
  for (i=0; i<x.length; i++) { 
    max_label = 0;
    for (j=0; j<x.label_size; j++) {
      score = inner_product(x, sm, j*x.feature_size, x.feature_size);
      index = x.label_size*x.feature_size;
      score += inner_product(x, sm, index+j*x.word_feature_size, x.word_feature_size);
      index += x.label_size*x.word_feature_size;
      score += inner_product(x, sm, index+j*x.predicate_feature_size, x.predicate_feature_size);      

      if (score>max_score) {
        max_score = score;
        max_label = j;
      }
    }
    y->labels[i]=max_label;
  }
  //return y;
}

void initialization_joint_helper(PATTERN x, LABEL *y, STRUCTMODEL *sm, int gold) {

}

void gibbs_sampling_helper(PATTERN x, LABEL *y, LABEL ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int useLoss) {

  //cout<<"gibbs_sampling_helper"<<endl;

  int i, j, k, l, m, n, index, p, q;
  double scores[x.label_size];
  //vector<double> scores (x.label_size, 0.0);
  bool first=true;
  int prev_labels[x.length];
  double max_score;
  int max_index; 
 
  for (i=0;i<x.length;i++){
    prev_labels[i]=y->labels[i]; 
  }

  // iteratively update 
  while (converge(y->labels,prev_labels,x.length,0.01)==0||first==true) {

  for (i=0; i<x.length; i++) {prev_labels[i]=y->labels[i];}

  first = false;
  for (i=0; i<x.identical_length; i++) {
    for (j=0; j<x.structToIndexSize[i];j++){
      for (k=0; k<x.label_size; k++) {scores[k]=0.0;}
      p = x.structToIndex[i][j];
      // for each node in data:
      for (k=0; k<x.label_size; k++) {
	index = 1;
        // x-y label
        for (l=0; l<x.feature_size; l++) {
          scores[k] += sm->w[index+k*x.label_size+l]*x.features[p][l];
        }
        index += x.label_size*x.feature_size;
        for (l=0; l<x.word_feature_size; l++) {
          scores[k] += sm->w[index+k*x.label_size+l]*x.wordFeatures[p][l];
        } 
        index += x.label_size*x.word_feature_size;
        for (l=0; l<x.word_feature_size; l++) {
          scores[k] += sm->w[index+k*x.label_size+l]*x.predicateFeatures[p][l];
        } 

        // yt label
	index+=x.label_size*x.feature_size;
        for (m=0; m<x.identical_length; m++) {
	  if (x.structures[i][m]==1) {
	    for (n=0; n<x.structToIndexSize[m];n++) {
              q = x.structToIndex[m][n];
              // check that they should be in same parse tree. (same predicate)
              if (x.predicates[p]==x.predicates[q]){
                scores[k] += sm->w[index+k*x.label_size+y->labels[q]];
              }
            }
          }
        }

        // yl label
	index+=pow(x.label_size,2);
        for (m=0; m<x.structToIndexSize[i];m++) {
          q = x.structToIndex[i][m];
          // check that they are the same word but in different parse tree. (different predicate)
          if (strcmp(x.predicates[p],x.predicates[q])!=0) {
            scores[k] += sm->w[index+k*x.label_size+y->labels[q]];
          }
        }

        if (useLoss==1) {
          if (ybar.labels[p]!=y->labels[p]) {
            scores[k]+=(1.0/double(x.length)); //lossFunction(y,ybar);
          }
        }

      }
      // find max from scores, update y
      max_index = 0;
      max_score = scores[0];
      for (k=0;k<x.label_size;k++) {
        if (scores[k] >= max_score) {
          max_score = scores[k];
          max_index = k;
        }
      }
      y->labels[p]=max_index;
    } // j 
  } // i
  } // while
}

void gibbs_sampling_joint_helper(PATTERN x, LABEL *y, LABEL ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int useLoss) {
}

void beam_search_helper(PATTERN x, LABEL *y, LABEL ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int useLoss) {
  
}

void beam_search_joint_helper(PATTERN x, LABEL *y, LABEL ybar, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, int useLoss) {
  
}

int find_element(char **list, int list_length, char *target) {
  int i;
  int ind = -1;
  for (i=0; i<list_length; i++){
    if (strcmp(target, list[i])==0) {
      ind = i;
    }
  }
  return ind;
}

double* feature_vector_helper(PATTERN x, LABEL y, STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm){

  //cout<<"feature_vector_helper"<<endl;

  int i, j, k, l, index, p, q, ind;
  double* feature_vector = (double*)malloc(sizeof(double)*sm->sizePsi);
  //initialization
  for (i=0; i<sm->sizePsi; i++){
    feature_vector[i]=0.0;
  }

  if (sparm->joint==2){
    // todo

  } else {
    // x-y part joint features
    index = 0;
    for (i=0; i<x.length; i++) {
      for (j=0; j<x.feature_size; j++) {
        feature_vector[index+y.labels[i]*(x.feature_size)+j]+=x.features[i][j];
      }
    }
    index += x.label_size*x.feature_size;

    // x-y part word features
    for (i=0; i<x.length; i++){
      // TODO find word 
      ind = find_element(x.identicalWords, x.word_num, x.words[i]);
      if (ind==-1) {
        printf("fatal error for words: %s\n",x.words[i]);
        ind = 0;
      }

      for (j=0; j<x.word_feature_size; j++){
        feature_vector[index+y.labels[i]*(x.word_feature_size)+j]+=x.wordFeatures[ind][j];
      }
    }
    index += x.label_size*x.word_feature_size;
    // x-y part predicate features
    for (i=0; i<x.length; i++){
      // TODO find predicate
      ind = find_element(x.identicalPredicates, x.predicate_num, x.predicates[i]);
      if (ind==-1) {
        printf("fatal error for predicates: %s\n",x.predicates[i]);
      }
      for (j=0; j<x.predicate_feature_size; j++){
        feature_vector[index+y.labels[i]*(x.predicate_feature_size)+j]+=x.predicateFeatures[ind][j];
      }
    }
    index += x.label_size*x.predicate_feature_size;

    // yt part
    for (i=0; i<x.identical_length; i++) {
      for (j=0; j<x.identical_length; j++) {
        if (x.structures[i][j]==1){
          for (k=0; k<x.structToIndexSize[i]; k++) {
            for (l=0; l<x.structToIndexSize[j]; l++) {
              p = x.structToIndex[i][k];
              q = x.structToIndex[j][l];
              if (strcmp(x.predicates[p],x.predicates[q])==0) {
                feature_vector[index+y.labels[p]*x.label_size+y.labels[q]]+=1;
              }
            }
          }
        }
      }
    }

    // yl part
    index += pow(x.label_size,2);
    for (i=0; i<x.identical_length; i++) {
      for (k=0; k<x.structToIndexSize[i]; k++) {
        for (l=0; l<x.structToIndexSize[i]; l++) {
          p = x.structToIndex[i][k];
          q = x.structToIndex[i][l];
          if (x.predicates[p]!=x.predicates[q]) {
            // every link will count for twice.
            feature_vector[index+y.labels[p]*x.label_size+y.labels[q]]+=0.5;
          }
        }
      }
    }

  }
  return feature_vector;
} 

int converge(int* y1, int* y2, int size, double constraint) {
  int i, diff;
  diff=0;
  for (i=0; i<size; i++) {
    if (y1[i]!=y2[i]) {
      diff++;
    }
  }

  if (double(diff)/double(size) > constraint) {
    return 0;
  }
  return 1;
}
