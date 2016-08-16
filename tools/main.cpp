//
//  main.cpp
//  10xVarCall
//
//  Created by Will Richards  on 6/16/16.
//  Copyright © 2016 WillRichards. All rights reserved.
//

#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <utility>
#include <unordered_map>
#include <string>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <limits.h>
#include <tuple>


#include <shared/bamtools_global.h>
#include <api/api_global.h>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include <api/SamReadGroupDictionary.h>


using namespace std;
using namespace BamTools;


//typedef std::unordered_map<std::string,std::pair<int32_t, vector<pair< int32_t, int32_t> > > > hashmap_t;
typedef vector< vector<pair<int32_t,int32_t > > > cluster_t;
typedef vector< pair<string, cluster_t> > clusterList_t;
typedef vector< pair <BamAlignment, vector< tuple <int,char,char> > > > molecule_t;
typedef vector< tuple < int, char, char> > variants_t;
typedef unordered_map<string, vector<molecule_t> > hashmap_t;

//executes a command in the shell, and writes command output as funciton output
string exec(char const* cmd) {
  char buffer[128];
  string result = "";
  FILE* pipe = popen(cmd, "r");
  if (!pipe) throw runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
      if (fgets(buffer, 128, pipe) != NULL)
	result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  return result;
}



//Merges a set of intervals that overlap
vector<pair<int, int> >intervalMerge(vector<pair<int, int> > &ranges) {
  vector<pair<int, int> > result;
  sort(ranges.begin(),ranges.end());
  vector<pair<int, int> >::iterator it = ranges.begin();
  pair<int,int> current = *(it)++;
  while (it != ranges.end()){
    if (current.second > it->first){ 
      current.second = std::max(current.second, it->second); 
    } else {
      result.push_back(current);
      current = *(it);
    }
    it++;
  }
  result.push_back(current);
  return result;
}



//prints out the hashmap of barcodes and corresponding positions to cout
//void distance(hashmap_t &hashmap) {
//  hashmap_t::iterator hashIt;
//  vector<pair<int32_t, int32_t> >::iterator vecIt;
//  ofstream distance;
//  distance.open("distance.csv");
//  int i;
//  for (hashIt = hashmap.begin(); hashIt != hashmap.end(); hashIt++) {
//     for (vecIt = hashIt->second.second.begin(); vecIt != hashIt->second.second.end(); vecIt++) {
//     int dist = (vecIt+1)->first - vecIt->first;
//     if (dist != 0) {
//       distance  << dist << endl;
//  }
//  }
// }
//}


void printTokenizedCig(vector<pair<char, string> > &cigVec) {
  vector<pair<char,string> >::iterator cigIt;

  for (cigIt = cigVec.begin(); cigIt != cigVec.end(); cigIt++) {
    cout << "character: " << cigIt->first << "  number: " << cigIt->second << endl;
  }
}



vector<pair<char, string> > tokenizeCig(string &cigarString) {
 vector<pair<char, string> > cigVec;
 string::iterator cigIt;
 pair<char, string> cigPair;

  for (cigIt = cigarString.begin(); cigIt != cigarString.end(); cigIt++) {
       if (isalpha(*cigIt)) {
	 if (isalpha(cigPair.first)) {
      cigVec.push_back(cigPair);
	 }
      cigPair.first = *cigIt;
      //cout << "debug: " << cigPair.first << endl;
       cigPair.second = "";
    }

    if(isdigit(*cigIt)) {
      cigPair.second.push_back(*cigIt);
      // cout << "bebug2: " << cigPair.first << ", " << cigPair.second << endl;
    }
  }
  cigVec.push_back(cigPair);
  return cigVec; 
 
}





//calculates and prints out statistics of the meta-data of the barcode hash-map
//void barCodeStats(hashmap_t &hashmap) {
//  int i;
//  hashmap_t::iterator hashIt;
//  std::vector<int32_t> freqList;
//  vector<int32_t>::iterator vecIt;
//
//  for (hashIt = hashmap.begin(); hashIt != hashmap.end(); hashIt++) {
//    freqList.push_back(hashIt->second.first); 
//  }

//  int averageFreq = accumulate(freqList.begin(), freqList.end(), 0.0 )/ freqList.size();
//  vecIt = max_element(freqList.begin(), freqList.end());
//  cout << "averge freq is: " << averageFreq << endl;
//  cout << "max freq is:  "  << *vecIt << endl;
//}



//Clusters the alignments corresponding to a barcode based on position
//clusterList_t clusters(hashmap_t &hashmap) {
//  hashmap_t::iterator hashIt;
//  vector<pair<int32_t, int32_t> >::iterator vecIt;
//  int i;
//  vector< pair<string, cluster_t> > clusterList;
 
//   for (hashIt = hashmap.begin(); hashIt != hashmap.end(); hashIt++) {
//    pair<string, cluster_t> barCodeCluster;
//    vector<pair<int32_t, int32_t> > initCluster;
//     initCluster.push_back(make_pair(hashIt->second.second[0].first,hashIt->second.second[0].second));
//    barCodeCluster.first = hashIt->first;
//    barCodeCluster.second.push_back(initCluster);
   
    //        cout << hashIt->second.second[0].first << endl;

//    int distThresh = 100000;
//     for (vecIt = hashIt->second.second.begin(); vecIt != hashIt->second.second.end(); vecIt++) {
//         if ((*vecIt).first < (barCodeCluster.second.back().back().first + distThresh)) {
//         barCodeCluster.second.back().push_back(*vecIt);
//      	barCodeCluster.first = hashIt->first;
 
//          } else {
//      	vector<pair< int32_t, int32_t> > cluster;
//      	cluster.push_back(*vecIt);
//	//cout << vecIt->first << endl;
//      	barCodeCluster.second.push_back(cluster);
//        }
//       sort(barCodeCluster.second.begin(), barCodeCluster.second.end());

//       }
    
//       clusterList.push_back(barCodeCluster);       
//     }
//  return clusterList;
//}


void coverage(clusterList_t  &clusterList) {                                                                                                                                                 
   ofstream coverage;                                                                                                                                                                              
   coverage.open("coverage.csv");                                                                                                                                                                  
   clusterList_t::iterator vecIt;
   cluster_t::iterator clustIt;                                                                               
   vector<pair<int32_t, int32_t> >::iterator readIt;   
   float count;
   double minRead;
   double inf = numeric_limits<double>::infinity();
   float maxRead;
   for (vecIt = clusterList.begin(); vecIt != clusterList.end(); vecIt++) {
     for (clustIt = vecIt->second.begin(); clustIt!= vecIt->second.end(); clustIt++) {
	 for (readIt = clustIt->begin(); readIt!=clustIt->end()-1; readIt++) {
	   float first = readIt->first;
	   float second = readIt->second;

	   if ( second > maxRead) {
	     maxRead = second;
	   }


	   if (first < minRead) {
	     minRead = first;
	   }


	   //cout << "first is: " << first << "second is: " << second << endl;                                                                                                                                                                       
	   count = count + (second - first);
	  
	   }
	 float size = maxRead - minRead;
	 float cov = count / size;
	 coverage << cov << endl;
	 //cout << "size is " << size << endl;
	 //cout << "count is " << count << endl; 
	 }
     count = 0;
     minRead = inf;
     maxRead = 0;
     
       }
     }


double simpleMatch(string barcode1, string barcode2) {
  double smc = 0;
  for (string::size_type i = 0; i < barcode1.size(); i++){
    if (barcode1[i] == barcode2[i]) {
       smc++;
       // cout << smc << endl;
    }
  }
   smc = smc / barcode1.length();
   // cout << smc;
  return smc;
}


//hashmap_t filterBarcodes(hashmap_t &hashmap) {
//  hashmap_t::iterator hashIt;
//  hashmap_t filteredHashmap;
//   for(hashIt = hashmap.begin(); hashIt != hashmap.end(); hashIt++) {
//     if (hashIt->second.first == 1) {
       //cout << "found it!! \n";
//      filteredHashmap[hashIt->first] = hashIt->second;
//            }
//     }
//   return filteredHashmap;
//}

//void barCodeReAllign(hashmap_t &hashmap, hashmap_t &filteredHashmap) {
//  hashmap_t::iterator hashIt1;
//  hashmap_t::iterator hashIt2;

//  for (hashIt2 = filteredHashmap.begin(); hashIt2 != filteredHashmap.end(); hashIt2++) {
    //cout << hashIt2->second.second[0].first << endl;
    //cout << hashIt2->first << endl;
    //   cout << "found the rema"


//    for (hashIt1 = hashmap.begin(); hashIt1 != hashmap.end(); hashIt1++) {
      // cout << hashIt1->first << ", " << hashIt2->first << endl;
//      double smc = simpleMatch(hashIt2->first, hashIt1->first);
      //cout << smc << endl;

//      if ((smc >= 0.92) && (smc != 1)) {
//	cout << "smc: " << smc << endl;
//	  cout << "found the remap!: " << hashIt1->first << ", " << hashIt2->first << endl;
//	  cout << hashIt2->second.second[0].first << ", " << hashIt1->second.second[0].first << endl;
//	  break;
//      }	
//   }
//  }
//}



//generates and writes out statistics about the meta-data of the clusters
void clusterAnal(clusterList_t  &clusterList) {
  clusterList_t::iterator vecIt;
  cluster_t::iterator clustIt;
  vector<pair< int32_t, int32_t> >::iterator readIt;
  vector< vector<int> > totalClusterCount;
  vector< vector<int> >::iterator totalCCIt;
  vector<int>::iterator clusterCountIt;
  vector<vector<pair<int, int> > > totalIntList;
  vector<pair< int, int> > intervalList;
  vector<pair<int, int> > moleculeList;
  vector<vector<pair<int, int> > >totalMList;
  int size;
  ofstream clusterStats;
  ofstream moleculeHist2;
  clusterStats.open("clusterStats.csv");
  moleculeHist2.open("moleculeHistogram2.csv");
  ofstream readsHist2;
  readsHist2.open("readsHistogram2.csv");

  for (vecIt = clusterList.begin(); vecIt!= clusterList.end(); vecIt++) {
    int clustCount = 0;
    vector<int> clusterCountList;
       for (clustIt = vecIt->second.begin(); clustIt!= vecIt->second.end(); clustIt++) {
      clustCount++;
      int clustSize = 0;
      int maxRead = 0;
      int  minRead = INT_MAX;
      pair<int, int> readPair;
      pair<int, int> moleculePair;   
    
      //loops through clusters, and calculates the ;ength of each cluster generated  
    for (readIt = clustIt->begin(); readIt!=clustIt->end(); readIt++) {
      clustSize = clustSize + 1;

      if ((*readIt).second > maxRead) {
        maxRead = (*readIt).second;
      	}
      if ((*readIt).first < minRead) {
       	minRead = ((*readIt).first);
       	}
     }
       readsHist2 << clustSize << endl;
       // cout << " # of reads in clust: " << clustSize << endl;
       //out << "start of molecule: " << minRead << endl << "end of molecule: " << maxRead << endl;
       // cout << "molecule size: " << maxRead - minRead << endl;
       size = maxRead - minRead;
       moleculeHist2  << size << endl;
    }
  }
 
  readsHist2.close();
  moleculeHist2.close();

  // unordered_map<int, int>::iterator readsIt;
  //  ofstream readsHist;
  // readsHist.open("readsHistogram.csv");
  // for (readsIt = readsHistogram.begin(); readsIt != readsHistogram.end(); readsIt++) {
  //  readsHist << readsIt->first << "," << readsIt->second << endl;
  //  }

  //double clustPerBC = 0;
  // double readsPerCluster = 0;
  //   for (totalCCIt = totalClusterCount.begin(); totalCCIt!= totalClusterCount.end(); totalCCIt++) {
  //cout << "Clusters per barcode: " << totalCCIt->size() << endl;
  //  clustPerBC += totalCCIt->size();
  //   for (clusterCountIt = totalCCIt->begin(); clusterCountIt!= totalCCIt->end(); clusterCountIt++) {
  //  cout << "reads per cluster: " << *clusterCountIt << endl;
  //    readsPerCluster += *clusterCountIt;
  //    }
  //    }
  


  //  clusterStats << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  //  clusterStats << "total reads: " << readsPerCluster << endl;
  //  clusterStats << "average reads per molecule: "<< readsPerCluster / clustPerBC << endl;
  // clusterStats  << "total molecule: " << clustPerBC << endl; 
  // clustPerBC = clustPerBC / clusterList.size();
  // clusterStats << "avg molecules per barcode: " << clustPerBC << endl;
  //  clusterStats << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    
  //  clusterStats.close();

}

string getFastaPos(BamAlignment &al) {
  string  start = to_string(al.Position+1);
  string end = to_string(al.GetEndPosition());
  string pos = start + "-" + end;
  // cout << pos << endl;
  return pos;
}

string getCig(std::vector<CigarOp> & cigString) {
  std::vector< CigarOp >::iterator cigIt;
  string cigar;
    for (cigIt = cigString.begin(); cigIt != cigString.end(); cigIt++) {
      if ((*cigIt).Type == 'P' || (*cigIt).Type == 'H' || (*cigIt).Type == 'N') {
	
      }
     cigar += (*cigIt).Type + to_string((*cigIt).Length);
     
    }
    //cout << cigar << endl;
    return cigar;
}
  
	

vector<tuple<int, char, char> >  variantDetector(vector<pair<char, string> > &tokenizedCig, string &refSeq, string &readSeq, int &pos) {
  vector<tuple<int, char, char> > variants;
  vector<pair<char, string> >::iterator tokenIt;
  string::iterator refIt = refSeq.begin();
  string::iterator readIt = readSeq.begin();
 
  for(tokenIt = tokenizedCig.begin(); tokenIt != tokenizedCig.end(); tokenIt++) {
    if (tokenIt->first == 'S' || tokenIt->first == 'I') {
      advance(readIt, stoi(tokenIt->second));
    }

  if (tokenIt->first == 'D') {
    advance(refIt, stoi(tokenIt->second));
    pos += stoi(tokenIt->second);
  }

    if (tokenIt->first == 'M') {   
      for(int i = 0; i < stoi(tokenIt->second); i++) {
	pos += 1;
	if(*refIt == *readIt || *refIt == 'N') {
	  cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
	}
	else {
	  switch(*refIt) {
	  case 'R':
	    if (*readIt == 'A' || *readIt == 'G') {
	      cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
	    }
	    else {
	      cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
	      cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
	      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	      variants.push_back(make_tuple(pos, *refIt, *readIt));
	    }
	     break;
	 
	  case 'Y':
            if (*readIt== 'C' || *readIt == 'T') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;

	  case 'K':
	    if (*readIt== 'G' || *readIt == 'T') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;

	  case 'M':
	    if (*readIt== 'A' || *readIt == 'C') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;

	  case 'S':
	    if (*readIt== 'C' || *readIt == 'G') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;

	  case 'W':
	    if (*readIt== 'A' || *readIt == 'T') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;
 
	  case 'B':
	    if (*readIt != 'A') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;
	    
	  case 'D':
	    if (*readIt != 'C') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;
	    
	  case 'H':
	    if (*readIt != 'G') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;
	    
	  case 'V':
	    if (*readIt != 'T') {
              cout << "Reference: " << *refIt << " Read: " << *readIt << endl;
            }
            else {
              cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
              cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
              cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            }
	    break;

	  case 'A':
	  case 'C':
	  case 'G':
	  case 'T':
	    cout << "~~~~~~~~~~~~~~~~~~~Variant Detected~~~~~~~~~~~~~~~~~~~~~~" << endl;
	    cout << "Reference: " << *refIt << " Read: " << *readIt << " at pos: " << pos << endl;
	    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	    break;
	  }	  
	}
	refIt++;
	readIt++;
      }
    } 
  }
      
  return variants;
}

bool sortReads(pair<BamAlignment, tuple<int, char, char> > &read1, pair<BamAlignment, tuple<int, char, char> > &read2) {
  return read1.first.Position < read2.first.Position;
}




//main function where we loop through all allignments, and stores data of interest to a hashmap
//wes can then call other functions on our hashmap
int main() {
    
  const std::string &filename = "/uufs/chpc.utah.edu/common/home/u0401321/TenX/bams/HG00513_WGS_phased_possorted_bam.bam"; 
  string fasta1 = " /uufs/chpc.utah.edu/common/home/marth-d1/data/software/chpc/kingspeak.peaks/gkno_launcher/tools/vcflib/fastahack/fastahack -r chr1:";
  string fasta2 = " /uufs/chpc.utah.edu/common/home/marth-d1/data/software/chpc/kingspeak.peaks/gkno_launcher/resources/homo_sapiens_grch38/current/GRCh38_full_analysis_set_plus_decoy_hla.fa";
  string outputFilename;
  
  BamTools::BamReader reader;
  reader.Open(filename);
  cout << filename << "\n";
  if(!reader.IsOpen()) {
    cout<<"{\"status\":\"error\", \"message\":\"Cannot open the specified file\"}"<<endl;
    exit(1);
  }
    
  const BamTools::RefVector refVector = reader.GetReferenceData();
  //cout << refVector.size() << "\n";
  map<int32_t, string> chromIDNameMap;
  for(size_t i=0; i<refVector.size(); i++) {
    chromIDNameMap[reader.GetReferenceID(refVector[i].RefName)] = refVector[i].RefName;
    // cout << refVector[i].RefName  << "\n";
  }
    
  const SamHeader header = reader.GetHeader();
  const RefVector references = reader.GetReferenceData();

  BamAlignment al;
  int i = 0;
  hashmap_t hashMap;  
  string cigarString;
  std::vector< CigarOp >::iterator cigIt;
 
  // iterate though alignnets, only keeping one with a high-ish quality
  // Generates a multiset of all barcodes
  while ( reader.GetNextAlignment(al) && i < 100000) {
    if ( al.MapQuality >= 10) {
      string  bx;
      getFastaPos(al);
      //cout << al.AlignedBases << endl;
      cigarString = getCig(al.CigarData);
      vector<pair<char, string> > tokenizedCig = tokenizeCig(cigarString);
      // printTokenizedCig(tokenizedCig);
      string pos = getFastaPos(al);
      string cmd1 = fasta1 + pos + fasta2;
      const char* cmd2 = cmd1.c_str();
      string execute =  exec(cmd2);
      //cout << execute;
      //cout << al.AlignedBases << endl;
      //cout << al.QueryBases << endl << endl;  
      // variantDetector(tokenizedCig, execute, al.QueryBases, al.Position); 
      if (al.GetTag("BX", bx)) {
	//if key already exists, add to freq count and position list
	variants_t variants = variantDetector(tokenizedCig, execute, al.QueryBases, al.Position);
       	if (hashMap.count(bx) == 1) {
	  vector<molecule_t>::iterator molIt;
	  // molecules = hashMap.find(bx);
	  int check = 0;
       	  for (molIt = hashMap.find(bx)->second.begin(); molIt != hashMap.find(bx)->second.end(); molIt++) {
	    if (al.Position < (*molIt).back().first.Position + 175000 && al.Position > molIt->back().first.Position - 175000) {
	      (*molIt).push_back(make_pair(al, variants));
	      check = 1;
	    }
	  }
	  //sort(molIt->begin(), molIt->end(), sortReads);
	  if (!check) { 
	    molecule_t newMolecule;
	    newMolecule.push_back(make_pair(al, variants));
	    hashMap.find(bx)->second.push_back(newMolecule);
  }
}
	//if new key, create new pair of frequency count and position list
	else {
	  // std::vector<pair<int32_t,int32_t> > positionList;
	  // positionList.push_back(make_pair(al.Position, al.GetEndPosition()));
	  vector<molecule_t> molecules;
	  //pair<int, int> posPair = make_pair(al.Position, al.GetEndPosition());
	  molecule_t molecule;
	  molecule.push_back(make_pair(al, variants));
	  molecules.push_back(molecule);		      
	  hashMap[bx] = molecules;
	}
      }
      i++;
  }
 }
 
  //  distance(hashMap);
 // barCodeStats(hashMap);
  //  clusterList_t clustList =  clusters(hashMap);
  //  coverage(clustList);
 //printClusters(clustList);
 // distanceStats(clustList);
 // clusterAnal(clustList);
 // hashmap_t filteredHashmap = filterBarcodes(hashMap);
 // barCodeReAllign(hashMap, filteredHashmap);
  reader.Close();
  string f = "M3D1M32S116";
  vector<pair<char, string> > cig = tokenizeCig(f);
  //printTokenizedCig(cig);
 
  
  return 0;   

}
