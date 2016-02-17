// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#define cimg_use_jpeg
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

void print_descriptor(const vector<SiftDescriptor> descriptors, CImg<double> &input) {
	for(int i=0; i<descriptors.size(); i++)
	  {
	    cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
	    for(int l=0; l<128; l++)
	      cout << descriptors[i].descriptor[l] << "," ;
	    cout << ")" << endl;

	    for(int j=0; j<5; j++)
	      for(int k=0; k<5; k++)
		if(j==2 || k==2)
		  for(int p=0; p<3; p++)
		  	if(descriptors[i].col+k < input.width() && descriptors[i].row+j < input.height())
		  		input(descriptors[i].col+k, descriptors[i].row+j, 0, p)=0;

	  }
}

int euclidean(const SiftDescriptor descriptor1, const SiftDescriptor descriptor2) {
	int sum = 0;
	for(int i=0;i<128;i++) {
		sum += pow((descriptor1.descriptor[i] - descriptor2.descriptor[i]), 2);
	}
	return sqrt(sum);
}

void calculate_distance_and_matches(vector<SiftDescriptor> &matches, vector<int> &distances, const vector<SiftDescriptor> descriptors, const vector<SiftDescriptor> descriptors2) {
	int dist, min_dist;
	for(int i=0;i<descriptors.size();i++) {
		min_dist = euclidean(descriptors[i], descriptors2[0]);
		matches[i] = descriptors2[0];
		for(int j=1;j<descriptors2.size();j++) {
			dist = euclidean(descriptors[i],descriptors2[j]);
			if(dist < min_dist) {
				min_dist = dist;
				matches[i] = descriptors2[j];
			}
		}
		distances[i] = min_dist;
	}
}

int quick_partition(int* max_matches, vector<string>& images, int low, int hi) {
	int pivot = max_matches[hi];
	int l = low;
	for (int j=low;j<hi;j++) {
		if(max_matches[j] >= pivot) {
			max_matches[j] = max_matches[j]+max_matches[l];
			max_matches[l] = max_matches[j]-max_matches[l];
			max_matches[j] = max_matches[j]-max_matches[l];

			string temp1 = images[j];
			images[j] = images[l];
			images[l] = temp1;

			l = l+1;
		}
	}
	max_matches[l] = max_matches[l]+max_matches[hi];
	max_matches[hi] = max_matches[l]-max_matches[hi];
	max_matches[l] = max_matches[l]-max_matches[hi];

	string temp1 = images[hi];
	images[hi] = images[l];
	images[l] = temp1;
	return l;
}

void quicksort(int* max_matches, vector<string>& images, int low, int hi) {
	int p;
	if(low < hi) {
		p = quick_partition(max_matches, images, low, hi);
		quicksort(max_matches, images, low, p-1);
		quicksort(max_matches, images, p+1, hi);
	}
}


int main(int argc, char **argv)
{
  try {

    if(argc < 4) {
		cout << "Insufficent number of arguments; correct usage:" << endl;
		cout << "    a2 part_id question_number ..." << endl;
		return -1;
    }

    string part = argv[1];
    string question = argv[2];

    if(part == "part1") {
    	string inputFile = argv[3];
    	CImg<double> input_image(inputFile.c_str());
    	int imgwidth = input_image.width();
    	CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
		vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
		print_descriptor(descriptors, input_image);
    	if(question == "q1") {
    		string infile2 = argv[4];
    		CImg<double> input_image2(infile2.c_str());
    		gray = input_image2.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
			print_descriptor(descriptors2, input_image2);
			vector<SiftDescriptor> matches(descriptors.size());
			vector<int> distances(descriptors.size());
			calculate_distance_and_matches(matches, distances, descriptors, descriptors2);
			CImg<double> output_image = input_image.append(input_image2);
			const unsigned char color[] = { 255,128,64 };
			for(int i=0;i<descriptors.size();i++) {
				if(distances[i] < 100) {
					output_image.draw_line(descriptors[i].col, descriptors[i].row, descriptors2[i].col+imgwidth, descriptors2[i].row, color);
				}
			}
			output_image.get_normalize(0,255).save("output1.png");
			input_image.get_normalize(0,255).save("sift.png");
			input_image2.get_normalize(0,255).save("sift2.png");
    	} else if(question == "q2") {
    		int max_matches[argc-4];
    		vector<string> images(argc-4);
    		for(int i=4;i<argc;i++) {
    			string infile2 = argv[i];
	    		CImg<double> input_image2(infile2.c_str());
	    		gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
				print_descriptor(descriptors2, input_image2);
				vector<SiftDescriptor> matches(descriptors.size());
				vector<int> distances(descriptors.size());
				calculate_distance_and_matches(matches, distances, descriptors, descriptors2);
				int count = 0;
				for (int j=0;j<descriptors.size();j++) {
					if(distances[j] < 50) {
						count++;
					}
				}
				max_matches[i-4] = count;
				images[i-4] = argv[i];
    		}
    		quicksort(max_matches, images, 0, argc-5);

    		cout << "The match status is : \n";

    		for(int i=0;i<argc-4;i++) {
    			cout << "image: " << images[i] << " || matches: " << max_matches[i] << "\n";
    		}
    	} else {
    		cout<<"Wrong option...Please enter a question number";
    	}
    }
    else if(part == "part2") {
	// do something here!
    }
    else
      throw std::string("unknown part!");

    // feel free to add more conditions for other parts (e.g. more specific)
    //  parts, for debugging, etc.
  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}








