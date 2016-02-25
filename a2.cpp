// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.

//2#define cimg_use_jpeg
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
#include <math.h>

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

void apply_transformation(CImg<double> &input) {
	double transformation[3][3] = {0.907, 0.258, -182, -0.153, 1.44, 58, -0.000306, 0.000731, 1};
	double transformation_inverse[3][3] = {1.12467, -0.314677,222.941,0.108839,0.685059,-19.9247,0.000264587,-0.000597069,1.08278};
	//double transformation_inverse[3][3] = {1.397602, -0.391042,277.044,0.13525199999999998,0.8513080000000001,-24.76,0.00032879700000000003,-0.000741965,1.345554};
	
	double I[3][3];
	
	CImg<double> warped(input.width(), input.height());
	double sum = 0;
	int new_coord[] = {0,0,1};
	double out[] = {0,0,0};

	warped.fill(255, 255, 255);
	cout << input.data(25,25,1,1) << endl;
	cout << input(25,25,1,1) << endl;

	//this is only to check whether the inverse is correct
	/*for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				sum += transformation[i][k] * transformation_inverse[k][j];
			}
			cout << sum << endl;
			sum = 0;
		}
	}*/
	cout << new_coord[0] << "," << new_coord[1] << "," << new_coord[2] << endl;
	
	//for (int i = 0; i < input.height(); i++) {
		//for (int j = 0; j < input.width(); j++){
	for (int i = 0; i < 1; i++) {
		for (int j = 0; j < 1; j++){
			new_coord[0] = i;
			new_coord[1] = j;			
			for (int m = 0; m < 3; m++) {
				sum = 0;
				for (int n = 0; n < 3; n++) {
					sum += transformation_inverse[m][n] * new_coord[n];
				}
				cout << "\n " << abs(sum);
				out[m] = abs(sum);
			}
			cout << endl << i << "," << j << "," << out[0]/out[2] << "," << out[1]/out[2];
			//cout << ",";// << input(out[0], out[1], 1, 0);
			//warped(i,j,1, 0) = input(out[0], out[1], 1, 0);
			
			//write rgb to warped
			//for(int rgb = 0; rgb<3;rgb++) {
				//warped(i,j,0,rgb) = input(out[0], out[1], 0, rgb);
			//}
			
		}
	}
	
	cout << endl << endl << input.height() << "," << input.width() << endl;
	warped.save("part2_q1_warped.png");
	
}

void inverseWarp(CImg<double> input_image)//,double arr[3][3])
{
  int sum=0;
  CImg<double> warped(input_image.width(),input_image.height(),input_image.depth(),input_image.spectrum());
  warped.fill(255,255,255);
  int w = 1;
  double coordinates[3]={0,0,w};
  double out[3];
  //double transformation_inverse[3][3] = {1.12467, -0.314677,222.941,0.108839,0.685059,-19.9247,0.000264587,-0.000597069,1.08278};
  double transformation_inverse[3][3] = {1.12467, -0.314677,222.941,0.108839,0.685059,-19.9247,0,0,1};
  
  double m[3][3] = {0.907, 0.258, -182, -0.153, 1.44, 58, -0.000306, 0.000731, 1};
  double minv[3][3];
  double det = 	m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
				m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
				m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

	double invdet = 1 / det;

	minv[0][0] = invdet * (m[1][1] * m[2][2] - m[2][1] * m[1][2]);
	minv[0][1] = invdet * (m[0][2] * m[2][1] - m[0][1] * m[2][2]);
	minv[0][2] = invdet * (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
	minv[1][0] = invdet * (m[1][2] * m[2][0] - m[1][0] * m[2][2]);
	minv[1][1] = invdet * (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
	minv[1][2] = invdet * (m[1][0] * m[0][2] - m[0][0] * m[1][2]);
	minv[2][0] = invdet * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
	minv[2][1] = invdet * (m[2][0] * m[0][1] - m[0][0] * m[2][1]);
	minv[2][2] = invdet * (m[0][0] * m[1][1] - m[1][0] * m[0][1]);

  //double transformation_inverse[3][3] = minv;
  
  for(int x=0;x<3;x++){
    for(int y=0;y<3;y++)
		cout << minv[x][y] << " ";
	cout << endl;
  }
  cout << endl << endl;
  
  
	for(int x=0;x<input_image.width();x++)
    for(int y=0;y<input_image.height();y++)
        {
          coordinates[0]=x;
          coordinates[1]=y;
          for(int i=0;i<3;i++)
           {
              sum=0;
              for(int k=0;k<3;k++)
                sum=sum+transformation_inverse[i][k]*coordinates[k];
                out[i]=abs(sum);
            }
			
        if(out[0]<input_image.width() && out[0]>=0 && out[1]<input_image.height() && out[1]>=0)
          {
            warped(x,y,0,0)=input_image(out[0],out[1],0,0);
            warped(x,y,0,1)=input_image(out[0],out[1],0,1);
            warped(x,y,0,2)=input_image(out[0],out[1],0,2);
			//warped(x,y,w,0)=input_image(out[0],out[1],w,0);
            //warped(x,y,w,1)=input_image(out[0],out[1],w,1);
            //warped(x,y,w,2)=input_image(out[0]out[1],w,2);
          }
      }
  warped.save("part2-q1.png");
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

	/*if(part == "part1") {
    	if (question == "q3") {

    		// randomly picking input image
    		
      		char temp_result[1024];
      		int dirmove = chdir(result);
    		cout<<"picking random...\n";
    		int imgindex = (rand() % filecount)+1;
    		cout<< "picked image : "<<filelist[imgindex]<<"\n";
    		string inputFile = strcpy(temp_result, filelist[imgindex].c_str());
    		string place = inputFile.substr(0, inputFile.find("_"));

    		// Getting image vector for input image and forming descriptors

    		CImg<double> input_image(inputFile.c_str());
			int imgwidth = input_image.width();
			CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
			print_descriptor(descriptors, input_image);

			// Initialization for matching process

			int max_matches[filecount-1];
    		vector<string> images(filecount-1);
    		vector<SiftDescriptor> matches(descriptors.size());
			vector<int> distances(descriptors.size());
			for(int i=0;i<filecount;i++) {

				// Getting image vector for query image and forming descriptors

				string infile2 = filelist[i];
	    		CImg<double> input_image2(infile2.c_str());
	    		gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
				print_descriptor(descriptors2, input_image2);

				// Calculating distances and finding matches

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

    		// randomly picking input image

    		char temp_result[1024];
      		int dirmove = chdir(result);
    		cout<<"picking random...\n";
    		int imgindex = (rand() % filecount)+1;
    		cout<< "picked image : "<<filelist[imgindex]<<"\n";
    		string inputFile = strcpy(temp_result, filelist[imgindex].c_str());
    		string place = inputFile.substr(0, inputFile.find("_"));

    		// Getting image vector for input image and forming descriptors

    		CImg<double> input_image(inputFile.c_str());
			int imgwidth = input_image.width();
			CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
			print_descriptor(descriptors, input_image);

			// Quantizing vectors with k dimensions

			vector<vector<int> > fv = quantize_vectors(descriptors, 400.0, 20);

			for(int i=0;i<filecount;i++) {

				// Getting image and forming the descriptors

				string infile2 = filelist[i];
	    		CImg<double> input_image2(infile2.c_str());
	    		gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);

				// Quantizing with k dimensions

				vector<vector<int> > f1v = quantize_vectors(descriptors2, 400.0, 20);
				int match_count = 0;
				print_descriptor(descriptors2, input_image2);

				// Running through all the vectors to find matches and finding the least distance match

				for(int j=0;j<fv.size();j++) {
					int dist = 999999;
					for(int y=0;y<f1v.size();y++) {
						if(fv[j] == f1v[y]){
							int temp = euclidean(descriptors[j], descriptors2[y]);
							if(temp<dist) {
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

    		// Getting image vector for input image and forming descriptors
    		
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

				// Getting image and forming the descriptors

	    		string infile2 = argv[4];
	    		CImg<double> input_image2(infile2.c_str());
	    		gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
				print_descriptor(descriptors2, input_image2);

				// Calculating distances and finding matches

				vector<SiftDescriptor> matches(descriptors.size());
				vector<int> distances(descriptors.size());
				calculate_distance_and_matches(matches, distances, descriptors, descriptors2);

				// Forming output image and drawing lines

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

					// Getting image and forming the descriptors

	    			string infile2 = argv[i];
		    		CImg<double> input_image2(infile2.c_str());
		    		gray = input_image2.get_RGBtoHSI().get_channel(2);
					vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
					print_descriptor(descriptors2, input_image2);

					// Calculating distances and finding matches

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
    else*/ 
		if(part == "part2") {
		if (question == "q1") {
			cout << "Part2 - Question 1:" << '\n';
			string inputFile = argv[3];
			CImg<double> input_image(inputFile.c_str());
			
			//apply_transformation(input_image);
			inverseWarp(input_image);
		}
    }
    else
      throw std::string("unknown part!");

  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}//0
