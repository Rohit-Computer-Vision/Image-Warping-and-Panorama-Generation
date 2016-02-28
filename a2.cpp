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

void calculate_distance_and_matches(vector<SiftDescriptor> &matches, vector<double> &distances, const vector<SiftDescriptor> descriptors, const vector<SiftDescriptor> descriptors2) {
	int dist, min_dist, next_min = 0;
	for(int i=0;i<descriptors.size();i++) {
		min_dist = euclidean(descriptors[i], descriptors2[0]);
		matches[i] = descriptors2[0];
		for(int j=1;j<descriptors2.size();j++) {
			dist = euclidean(descriptors[i],descriptors2[j]);
			if(dist < min_dist) {
				next_min = min_dist;
				min_dist = dist;
				matches[i] = descriptors2[j];
			}
		}
		distances[i] = (double)min_dist/(double)next_min;
	}
	sort(distances.begin(), distances.end());
}

int quick_partition(double max_matches[], vector<string>& images, int low, int hi) {
	double pivot = max_matches[low];
	int left = low+1;
	int right = hi;
	bool complete = false;
	while(complete != true) {
		while(left <= right && max_matches[left] <= pivot)
			left++;
		while(max_matches[right] >= pivot && right >= left)
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

void quicksort(double max_matches[], vector<string>& images, int low, int hi) {
	int p;
	if(low < hi) {
		p = quick_partition(max_matches, images, low, hi);
		quicksort(max_matches, images, low, p-1);
		quicksort(max_matches, images, p+1, hi);
	}
}

vector<vector<int> > quantize_vectors(const vector<SiftDescriptor> descriptors, const vector<vector<double> > x_vec, float w, int k) {
	vector<vector<int> > result;
	for(int i=0;i<descriptors.size();i++) {
		vector<int> temp;
		for(int l=0;l<k;l++) {
			double sum = 0.0;
			for(int j=0;j<128;j++) {
				sum += x_vec[l][j] * descriptors[i].descriptor[j];
			}
			temp.push_back((int)(sum/w));
		}
		result.push_back(temp);
	}
	return result;
}

void inverseWarp(CImg<double> input_image, CImg<double> warped,
		double transformation[3][3], int flag = 0) {
	double sum = 0.0;
	warped.fill(255, 255, 255);
	double coordinates[3] = { 0, 0, 1.0 };
	double out[3] = { 0, 0, 0 };

	double transformation1[3][3] = { 0.907, 0.258, -182.0, -0.153, 1.44, 58.0,
			-0.000306, 0.000731, 1.0 };
	if (flag == 1) {
		transformation = transformation1;
	}

	//cout << transformation[0][0] << "\n";
	double t_inv[3][3];
	double det = transformation[0][0]
			* (transformation[1][1] * transformation[2][2]
					- transformation[2][1] * transformation[1][2])
			- transformation[0][1]
					* (transformation[1][0] * transformation[2][2]
							- transformation[1][2] * transformation[2][0])
			+ transformation[0][2]
					* (transformation[1][0] * transformation[2][1]
							- transformation[1][1] * transformation[2][0]);

	cout << "\nDeterminant:" << det << "\n\n";
	double invdet = 1 / det;

	t_inv[0][0] = invdet
			* (transformation[1][1] * transformation[2][2]
					- transformation[2][1] * transformation[1][2]);
	t_inv[0][1] = invdet
			* (transformation[0][2] * transformation[2][1]
					- transformation[0][1] * transformation[2][2]);
	t_inv[0][2] = invdet
			* (transformation[0][1] * transformation[1][2]
					- transformation[0][2] * transformation[1][1]);
	t_inv[1][0] = invdet
			* (transformation[1][2] * transformation[2][0]
					- transformation[1][0] * transformation[2][2]);
	t_inv[1][1] = invdet
			* (transformation[0][0] * transformation[2][2]
					- transformation[0][2] * transformation[2][0]);
	t_inv[1][2] = invdet
			* (transformation[1][0] * transformation[0][2]
					- transformation[0][0] * transformation[1][2]);
	t_inv[2][0] = invdet
			* (transformation[1][0] * transformation[2][1]
					- transformation[2][0] * transformation[1][1]);
	t_inv[2][1] = invdet
			* (transformation[2][0] * transformation[0][1]
					- transformation[0][0] * transformation[2][1]);
	t_inv[2][2] = invdet
			* (transformation[0][0] * transformation[1][1]
					- transformation[1][0] * transformation[0][1]);

	for (int x = 0; x < input_image.width(); x++) {
		for (int y = 0; y < input_image.height(); y++) {
			coordinates[0] = x;
			coordinates[1] = y;
			for (int i = 0; i < 3; i++) {
				sum = 0.0;
				for (int k = 0; k < 3; k++) {
					sum += t_inv[i][k] * coordinates[k];
				}
				out[i] = sum;
			}
			out[1] = out[1] / out[2];
			out[0] = out[0] / out[2];

			if (out[0] < input_image.width() && out[0] >= 0
					&& out[1] < input_image.height() && out[1] >= 0) {
				warped(x, y, 0, 0) = input_image(floor(out[0]), floor(out[1]),
						0, 0);
				warped(x, y, 0, 1) = input_image(floor(out[0]), floor(out[1]),
						0, 1);
				warped(x, y, 0, 2) = input_image(floor(out[0]), floor(out[1]),
						0, 2);
			}
		}
	}
}

void convert_to_3x3(CImg<double> h, double homography[3][3]) {
	int k = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (k < 8) {
				homography[i][j] = h(0, k++);
			}
		}
	}
	homography[2][2] = 1;
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

			double max_matches[filecount-1];
    		vector<string> images(filecount-1);
    		vector<SiftDescriptor> matches(descriptors.size());
			vector<double> distances(descriptors.size());
			for(int i=0;i<filecount;i++) {

				// Getting image vector for query image and forming descriptors

				string infile2 = filelist[i];
	    		CImg<double> input_image2(infile2.c_str());
	    		gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
				print_descriptor(descriptors2, input_image2);

				// Calculating distances and finding matches

				calculate_distance_and_matches(matches, distances, descriptors, descriptors2);
				int ind = descriptors2.size()/40;
				max_matches[i] = distances[ind];
				images[i] = filelist[i];
				// cout<<"the match"<<filelist[i]<<"   "<<max_matches[i]<<"\n";
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
			vector<vector<double> > x_vec;
			int k_val = 10;
			for(int l=0;l<k_val;l++) {
				vector<double> temp_doub;
				for(int i=0;i<128;i++) {
					temp_doub.push_back((double)(rand() % 100) / 100.0);
				}
				x_vec.push_back(temp_doub);
			}

			vector<vector<int> > fv = quantize_vectors(descriptors, x_vec, 500.0, k_val);

			double max_matches[filecount];
    		vector<string> images(filecount);

			for(int i=0;i<filecount;i++) {

				// Getting image and forming the descriptors

				string infile2 = filelist[i];
	    		CImg<double> input_image2(infile2.c_str());
	    		gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);

				// Quantizing with k dimensions

				vector<vector<int> > f1v = quantize_vectors(descriptors2, x_vec, 500.0, k_val);
				double match_ratio = 0.0;
				print_descriptor(descriptors2, input_image2);

				// Running through all the vectors to find matches and finding the least distance match

				vector <double> temp_dists;

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
					temp_dists.push_back((double)dist);
				}
				sort(temp_dists.begin(), temp_dists.end());
				int ind = descriptors2.size()/40;
				match_ratio = temp_dists[ind];
				max_matches[i] = match_ratio;
				images[i] = filelist[i];
				// cout<<"Matches: "<<match_ratio<<" for image: "<<infile2<<"\n";
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
				vector<double> distances(descriptors.size());
				calculate_distance_and_matches(matches, distances, descriptors, descriptors2);

				// Forming output image and drawing lines

				CImg<double> output_image = input_image.append(input_image2);
				const unsigned char color[] = { 255,128,64 };
				for(int i=0;i<descriptors.size();i++) {
					if(distances[i] < 50) {
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
	    		double max_matches[argc-4];
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
					vector<double> distances(descriptors.size());
					calculate_distance_and_matches(matches, distances, descriptors, descriptors2);
					int ind = descriptors2.size()/40;
					max_matches[i-4] = distances[ind];
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
    } else if(part == "part2") {
		if (question == "q1") {
			cout << "Part2 - Question 1:" << '\n';
			string inputFile = argv[3];
			CImg<double> input_image(inputFile.c_str());
			CImg<double> warped = input_image;
			double transformation[3][3];
			inverseWarp(input_image, warped, transformation, 1);
			warped.save("part2-q1.png");
		} else if (question == "q2") {
			cout << "Part2 - Question 2:" << '\n';
				if (argc < 4) {
					cout << "Insufficent number of arguments; correct usage:"
							<< endl;
					cout << "    a2 part_id question_number ..." << endl;
					return -1;
				}

				// Getting image vector for input image and forming descriptors

				string inputFile = argv[3];
				CImg<double> input_image(inputFile.c_str());
				int imgwidth = input_image.width();
				CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);
				print_descriptor(descriptors, input_image);

				// Getting image and forming the descriptors

				string infile2 = argv[4];
				CImg<double> input_image2(infile2.c_str());
				gray = input_image2.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> descriptors2 = Sift::compute_sift(gray);
				print_descriptor(descriptors2, input_image2);

				// Calculating distances and finding matches

				vector<SiftDescriptor> matches(descriptors.size());
				vector<double> distances(descriptors.size());
				calculate_distance_and_matches(matches, distances, descriptors,
						descriptors2);

				//col is x-coord, row is y-coord
				//cout << descriptors[10].col << "," << descriptors[10].row << endl;

				int img_rows = input_image.width();
				int img_cols = input_image.height();
				//declaring these as CImg bcoz solver has these parameters and return value
				CImg<double> trans(8, 8);
				CImg<double> h(1, 8);
				CImg<double> global_h(1, 8);
				double homography[3][3]; //to pass to inverseWarp
				CImg<double> new_coords(1, 8);
				CImg<double> warped = input_image;

				int i, r, k, rand_x, rand_y;
				int x[4], y[4], x_dash[4], y_dash[4];
				double image_x, image_y, image_x_dash, image_y_dash;
				double new_image_x_dash, new_image_y_dash, new_image_w_dash;
				int max_inliers = 0, inliers = 0;

				//initializations for elements that are 1
				for (int i = 0; i < 8; i += 2) {
					trans(i, 2) = 1;
					trans(i + 1, 5) = 1;
				}

				for (int i = 0; i < 500; i++) {
					for (r = 0; r < 4; r++) {
						// gives a random number between 0 and image width
						rand_x = rand() % img_rows;
						rand_y = rand() % img_cols;
						x[r] = descriptors[rand_x].col;
						y[r] = descriptors[rand_y].row;
						x_dash[r] = matches[rand_x].col;
						y_dash[r] = matches[rand_y].row;
					}
					//calculating the 1st matrix trans and 3rd matrix new_coords (lec09_slide40)
					for (k = 0; k < 8; k += 2) {
						//har ek ko ulta
						trans(k, 0) = x[i];
						trans(k, 1) = y[i];
						trans(k, 6) = -x[i] * x_dash[i];
						trans(k, 7) = -y[i] * x_dash[i];
						new_coords(0, k) = x_dash[i];

						trans(k + 1, 3) = x[i];
						trans(k + 1, 4) = y[i];
						trans(k + 1, 6) = -x[i] * y_dash[i];
						trans(k + 1, 7) = -y[i] * y_dash[i];
						new_coords(0, k + 1) = y_dash[i];
					}
					h = new_coords.solve(trans);
					convert_to_3x3(h, homography);
					for (int m = 0; m < descriptors.size(); m++) {
						image_x = descriptors[m].col;
						image_y = descriptors[m].row;
						image_x_dash = matches[m].col;
						image_y_dash = matches[m].row;

						new_image_x_dash = homography[0][0] * image_x
								+ homography[0][1] * image_y + homography[0][2];
						new_image_y_dash = homography[1][0] * image_x
								+ homography[1][1] * image_y + homography[1][2];
						new_image_w_dash = homography[2][0] * image_x
								+ homography[2][1] * image_y + homography[2][2];

						new_image_x_dash /= new_image_w_dash;
						new_image_y_dash /= new_image_w_dash;

						//check number of inliers and save best h
						if (new_image_x_dash >= image_x_dash - 20
								&& new_image_x_dash <= image_x_dash + 20) {
							if (new_image_y_dash >= image_y_dash - 20
									&& new_image_y_dash <= image_y_dash + 20) {
								inliers++;
							}
						}
					}
					if (max_inliers < inliers) {
						max_inliers = inliers;
						global_h = h;
					}

					/*for(int m = 0; m<3;m++)
					 for(int n = 0; n<3;n++)
					 cout << homography[m][n] << " ";
					 */

				}
				//final homography
				inverseWarp(input_image, warped, homography);
				warped.save("part2-q2.png");

		}
    }
    else
      throw std::string("unknown part!");

  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}//0
