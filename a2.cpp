// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.

#define cimg_use_jpeg
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>
#include <dirent.h>
#include <unistd.h>
#include <time.h>
#include <map>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

void print_descriptor(const vector<SiftDescriptor> descriptors, CImg<double> &input) {
	for(int i=0; i<descriptors.size(); i++)
	  {
	    // cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
	    // for(int l=0; l<128; l++)
	    //   cout << descriptors[i].descriptor[l] << "," ;
	    // cout << ")" << endl;

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

int quick_partition(int max_matches[], vector<string>& images, int low, int hi) {
	int pivot = max_matches[low];
	int left = low+1;
	int right = hi;
	bool complete = false;
	while(complete != true) {
		while(left <= right && max_matches[left] >= pivot)
			left++;
		while(max_matches[right] <= pivot && right >= left)
			right--;
		if(right < left) {
			complete = true;
		} else {
			swap(max_matches[left], max_matches[right]);
			iter_swap(images.begin()+left, images.begin()+right);
		}
	}
	swap(max_matches[low], max_matches[right]);
	iter_swap(images.begin()+low, images.begin()+right);
	return right;
}

void quicksort(int max_matches[], vector<string>& images, int low, int hi) {
	int p;
	if(low < hi) {
		p = quick_partition(max_matches, images, low, hi);
		quicksort(max_matches, images, low, p-1);
		quicksort(max_matches, images, p+1, hi);
	}
}

vector<vector<int> > quantize_vectors(const vector<SiftDescriptor> descriptors, float w, int k) {
	vector<vector<int> > x_vec;
	for(int i=0;i<descriptors.size();i++) {
		vector<int> temp;
		for(int l=0;l<k;l++) {
			float sum = 0.0;
			for(int j=0;j<128;j++) {
				sum += ((double)(rand() % 100) / 100.0) * descriptors[i].descriptor[j];
			}
			temp.push_back((int)(sum/w));
		}
		x_vec.push_back(temp);
	}
	return x_vec;
}


int main(int argc, char **argv)
{
  try {

  	DIR *dir;
	struct dirent *ent;
	char cwd[1024];
	char result[1024];
	vector<string> filelist;
	int filecount = 0;
	strcpy(result, getcwd(cwd, sizeof(cwd)));
	strcat(result, "/part1_images");
	if ((dir = opendir (result)) != NULL) {
  		while ((ent = readdir (dir)) != NULL) {
  			if (strstr(ent->d_name, ".jpg")) {
	    		filelist.push_back(ent->d_name);
	    		filecount++;
	    	}
  		}
  		closedir (dir);
	} else {
  		perror ("");
  		return EXIT_FAILURE;
	}

    string part = argv[1];
    string question = argv[2];
    srand (time(NULL));

    if(part == "part1") {
    	if (question == "q3") {
      		char temp_result[1024];
      		int dirmove = chdir(result);
    		cout<<"picking random...\n";
    		int imgindex = (rand() % filecount)+1;
    		cout<< "picked image : "<<filelist[imgindex]<<"\n";
    		string inputFile = strcpy(temp_result, filelist[imgindex].c_str());
    		string place = inputFile.substr(0, inputFile.find("_"));
    		CImg<double> input_image(inputFile.c_str());
			int imgwidth = input_image.width();
			CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
			print_descriptor(descriptors, input_image);
			int max_matches[filecount-1];
    		vector<string> images(filecount-1);
    		vector<SiftDescriptor> matches(descriptors.size());
			vector<int> distances(descriptors.size());
			for(int i=0;i<filecount;i++) {
				string infile2 = filelist[i];
	    		CImg<double> input_image2(infile2.c_str());
	    		gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
				print_descriptor(descriptors2, input_image2);
				calculate_distance_and_matches(matches, distances, descriptors, descriptors2);
				int count = 0;
				for (int j=0;j<descriptors.size();j++) {
					if(distances[j] < 150) {
						count++;
					}
				}
				cout<<filelist[i]<<"    "<<count<<" matches here \n";
				max_matches[i] = count;
				images[i] = filelist[i];
			}

			quicksort(max_matches, images, 0, filecount-1);

    		cout << "The match status is : \n";
    		int correct = 0;
    		for(int i=0;i<filecount;i++) {
    			cout << "image: " << images[i] << " \t|| matches: " << max_matches[i] << "\n";
    			if(i < 10 && (place.length() < images[i].length()) && (images[i].find(place) != std::string::npos)) {
    				correct++;
    			}
    		}

    		cout<<"Precision = "<<float(correct)/10.0<<"\n";

    	} else if(question == "q4") {
    		char temp_result[1024];
      		int dirmove = chdir(result);
    		cout<<"picking random...\n";
    		int imgindex = (rand() % filecount)+1;
    		cout<< "picked image : "<<filelist[imgindex]<<"\n";
    		string inputFile = strcpy(temp_result, filelist[imgindex].c_str());
    		string place = inputFile.substr(0, inputFile.find("_"));
    		CImg<double> input_image(inputFile.c_str());
			int imgwidth = input_image.width();
			CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
			print_descriptor(descriptors, input_image);
			vector<vector<int> > fv = quantize_vectors(descriptors, 500.0, 20);
			int max_matches[filecount-1];
    		vector<string> images(filecount-1);
    		vector<SiftDescriptor> matches(descriptors.size());
			vector<int> distances(descriptors.size());
			for(int i=0;i<filecount;i++) {
				string infile2 = filelist[i];
	    		CImg<double> input_image2(infile2.c_str());
	    		gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
				vector<vector<int> > f1v = quantize_vectors(descriptors2, 500.0, 20);
				int match_count = 0;
				print_descriptor(descriptors2, input_image2);
				for(int j=0;j<fv.size();j++) {
					int dist = 999999;
					for(int y=0;y<f1v.size();y++) {
						if(fv[j] == f1v[y]){
							// cout << "Matched....\n";
							int temp = euclidean(descriptors[j], descriptors2[y]);
							if(temp<dist) {
								// cout<<temp<<"\n";
								dist = temp;
							}
						}
					}
					if(dist < 400)
						match_count++;
				}
				cout<<"Matches: "<<match_count<<" for image: "<<infile2<<"\n";
			}
    	} else {
			string inputFile = argv[3];
			CImg<double> input_image(inputFile.c_str());
			int imgwidth = input_image.width();
			CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
			print_descriptor(descriptors, input_image);
	    	if(question == "q1") {
	    		if(argc < 4) {
					cout << "Insufficent number of arguments; correct usage:" << endl;
					cout << "    a2 part_id question_number ..." << endl;
					return -1;
	    		}
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
	    		if(argc < 4) {
					cout << "Insufficent number of arguments; correct usage:" << endl;
					cout << "    a2 part_id question_number ..." << endl;
					return -1;
	    		}
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
