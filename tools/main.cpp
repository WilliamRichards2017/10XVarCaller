//
//  main.cpp
//  10xVarCall
//
//  Created by Will Richards  on 6/16/16.
//  Copyright © 2016 WillRichards. All rights reserved.
//

#include <stdlib.h>
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


#include <shared/bamtools_global.h>
#include <api/api_global.h>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include <api/SamReadGroupDictionary.h>


using namespace std;
using namespace BamTools;


typedef std::unordered_map<std::string,std::pair<int32_t, vector<int32_t>>> hashmap_t;
typedef vector< vector<int32_t > > cluster_t;
typedef vector< pair<string, cluster_t> > clusterList_t;


//writes out the hashmap of barcodes and corresponding positions to a csv file
//first column is the unique barcode, all other colmns are positions corresponding to the barcodes
void hashMaptoFile(hashmap_t &hashmap) {
  hashmap_t::iterator hashIt;
  vector<int32_t>::iterator vecIt;
  int i;

  ofstream barCodeCSV;
  barCodeCSV.open("barCodeFreq.csv");

  for (hashIt = hashmap.begin(); hashIt != hashmap.end(); hashIt++) {
    barCodeCSV << hashIt->first << ","  << hashIt->second.first;
    for (vecIt = hashIt->second.second.begin(); vecIt != hashIt->second.second.end(); vecIt++) {
      barCodeCSV << ","  << vecIt[i];
    }
    barCodeCSV << "\n";
  }
}


//prints out the hashmap of barcodes and corresponding positions to cout
void printHashMap(hashmap_t &hashmap) {
  hashmap_t::iterator hashIt;
  vector<int32_t>::iterator vecIt;
  int i;
  for (hashIt = hashmap.begin(); hashIt != hashmap.end(); hashIt++) {
    cout << "barcode is: "  << hashIt->first << "  frequency is: "  << hashIt->second.first  << "  barcode positions are "; 
      for (vecIt = hashIt->second.second.begin(); vecIt != hashIt->second.second.end(); vecIt++) {
	cout << vecIt[i] << ", ";
       }
    cout << "\n";
  }
}

//calculates and prints out statistics of the meta-data of the barcode hash-map
void barCodeStats(hashmap_t &hashmap) {
  int i;
  hashmap_t::iterator hashIt;
  std::vector<int32_t> freqList;
  vector<int32_t>::iterator vecIt;

  for (hashIt = hashmap.begin(); hashIt != hashmap.end(); hashIt++) {
    freqList.push_back(hashIt->second.first); 
  }

 int averageFreq = accumulate(freqList.begin(), freqList.end(), 0.0 )/ freqList.size();
 vecIt = max_element(freqList.begin(), freqList.end());
 cout << "averge freq is: " << averageFreq << endl;
 cout << "max freq is:  "  << *vecIt << endl;
}



//Clusters the alignments corresponding to a barcode based on position
clusterList_t clusters(hashmap_t &hashmap) {
    hashmap_t::iterator hashIt;
    vector<int32_t>::iterator vecIt;
    int i;
    vector< pair<string, cluster_t> > clusterList;
 
    for (hashIt = hashmap.begin(); hashIt != hashmap.end(); hashIt++) {
      pair<string, cluster_t> barCodeCluster;
      vector<int32_t> initCluster;
      int32_t initInt = hashIt->second.second[0];

      initCluster.push_back(initInt);
      barCodeCluster.first = hashIt->first;
      barCodeCluster.second.push_back(initCluster);

      for (vecIt = hashIt->second.second.begin(); vecIt != hashIt->second.second.end(); vecIt++) {

	if ((vecIt[i] < (barCodeCluster.second.back().back() + 1000)) && (vecIt[i] > (barCodeCluster.second.back().back() - 1000))) {
	  barCodeCluster.second.back().push_back(vecIt[i]);
	  barCodeCluster.first = hashIt->first;
       }
	else {
	  vector<int32_t> cluster;
	  cluster.push_back(vecIt[i]);
	  barCodeCluster.second.push_back(cluster);
	}
      }
      clusterList.push_back(barCodeCluster);       
    }
   return clusterList;
}

//prints out generated clusters
//each barcode has a list of positions, where each cluster is seperated by '||'
void printClusters(vector< pair<string, cluster_t> >  &clusterList) {
  vector<pair<string, cluster_t> >::iterator vecIt;
  cluster_t::iterator clustIt;
  vector<int32_t>::iterator readIt;

  for (vecIt = clusterList.begin(); vecIt!= clusterList.end(); vecIt++) {
    cout << "barcode : " << vecIt->first << ", cluster(s): ";
    for (clustIt = vecIt->second.begin(); clustIt!= vecIt->second.end(); clustIt++) {
      cout << "||";
      for (readIt = clustIt->begin(); readIt!=clustIt->end(); readIt++) {
	cout  << *readIt << ", ";
      }
    }
    cout << "||" << endl;
  }
}

//Writes out the clusters assigned to each barcode
//Each line contains a barcode, as well as a list of clusters
//Each new cluster is marked with the '|' character
void clusterToFile(clusterList_t  &clusterList) {   
  ofstream clusters;
  clusters.open("clusters.csv");
   clusterList_t::iterator vecIt;
   cluster_t::iterator clustIt;                                                                                                         
   vector<int32_t>::iterator readIt;              

   for (vecIt = clusterList.begin(); vecIt!= clusterList.end(); vecIt++) {
     clusters << vecIt->first << ",";
     for (clustIt = vecIt->second.begin(); clustIt!= vecIt->second.end(); clustIt++) {
       clusters << "|";
       for (readIt = clustIt->begin(); readIt!=clustIt->end(); readIt++) {
	 clusters  << *readIt << ",";
       }
     }
     clusters << endl;
   }
}

//generates and writes out statistics about the meta-data of the clusters
void clusterAnal(clusterList_t  &clusterList) {
  clusterList_t::iterator vecIt;
  cluster_t::iterator clustIt;
  vector<int32_t>::iterator readIt;
  vector< vector<int> > totalClusterCount;
  vector< vector<int> >::iterator totalCCIt;
  vector<int>::iterator clusterCountIt;

  ofstream clusterStats;
  clusterStats.open("clusterStats.csv");


  for (vecIt = clusterList.begin(); vecIt!= clusterList.end(); vecIt++) {
    int clustCount = 0;
    vector<int> clusterCountList;
    for (clustIt = vecIt->second.begin(); clustIt!= vecIt->second.end(); clustIt++) {
      clustCount++;
      int clustSize = 0;
      for (readIt = clustIt->begin(); readIt!=clustIt->end(); readIt++) {
        clustSize++;
      }
        clusterCountList.push_back(clustSize);
    }
     totalClusterCount.push_back(clusterCountList);
   }

  double clustPerBC = 0;
  double readsPerCluster = 0;
  for (totalCCIt = totalClusterCount.begin(); totalCCIt!= totalClusterCount.end(); totalCCIt++) {
    //    cout << "Clusters per barcode: " << totalCCIt->size() << endl;
    clustPerBC += totalCCIt->size();
    for (clusterCountIt = totalCCIt->begin(); clusterCountIt!= totalCCIt->end(); clusterCountIt++) {
      //cout << "reads per cluster: " << *clusterCountIt << endl;
      readsPerCluster += *clusterCountIt;
    }
  }

  clusterStats << "total reads: " << readsPerCluster << endl;
  clusterStats << "average reads per cluster: "<< readsPerCluster / clustPerBC << endl;
  clusterStats << "total clusters: " << clustPerBC << endl; 
  clustPerBC = clustPerBC / clusterList.size();
  clusterStats << "avg clusters per barcode: " << clustPerBC << endl;

}



//main function where we loop through all allignments, and stores data of interest to a hashmap
//we can then call other functions on our hashmap
int main() {
    
  const std::string &filename = "/uufs/chpc.utah.edu/common/home/u0401321/TenX/bams/HG00512_WGS_phased_possorted_bam.bam"; 
    
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
    cout << refVector[i].RefName  << "\n";
  }
    
  const SamHeader header = reader.GetHeader();
  const RefVector references = reader.GetReferenceData();

  BamAlignment al;
  int i = 0;
  hashmap_t hashMap;  

  // iterate though alignnets, only keeping one with a high-ish quality
  // Generates a multiset of all barcodes
    while ( reader.GetNextAlignment(al)) {
     if ( al.MapQuality >= 10) {
        string  bx;
	
        if (al.GetTag("BX", bx)) {
	  //if key already exists, add to freq count and position list
          if (hashMap.count(bx) == 1) {
	    hashMap.find(bx)->second.second.push_back(al.Position);
       	    hashMap.find(bx)->second.first++;
	    }
	  //if new key, create new pair of frequency count and position list
	    else {
     	      std::vector<int32_t> positionList;
	      positionList.push_back(al.Position);
	       hashMap[bx] = make_pair(1, positionList);
	    }
	 }
        i++;
  }
 }

      // printHashMap(hashMap);
      // hashMaptoFile(hashMap);
      // barCodeStats(hashMap);
      clusterList_t clustList =  clusters(hashMap);
      // printClusters(clustList);
      // clusterToFile(clustList);
      clusterAnal(clustList);
 
  reader.Close();
  return 0;
    
}




