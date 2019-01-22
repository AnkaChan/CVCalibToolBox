#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "AC/AC.h"
#include "AC/cxxopts.hpp"
#include "CalibDef.h"

#include <direct.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <random>

#define NUM_THREADS 6
#define CB_SQUARE_SIZE 30
struct CalibrationConfigs
{
	std::string inImgPathL;
	std::string inImgPathR;

	std::string inImgSequenceFile;

	std::string inIntrinsicsFile;

	std::string outFile;

	bool stereoCalib = true;

	bool showCorners = false;

	bool outputSelectedFileSequence = false;
	std::string outImgSequenceFile;

	bool undistordImages = false;
	std::string outUndistortedImgpath;

	bool resampleInputs = false;
	int sampleNum = 100;
	enum ResampleMode
	{
		random,
		skip
	} resampleMode = random;

	int boardWidth = 9;
	int boardHeight = 6;

	std::string fileExt = "jpg";
	std::string cornersOutput = "";
	enum Mode
	{
		//Calibrate single camera's intrinsic parameters
		Single,
		//Calibrate two camera's intrinsic parameters
		Dual,
		//Calibrate two camera's intrinsic parameters, and apply stereo calibration on them
		DualStereo,
		//Only calibrate extrinsic parameters
		OnlyStereo
	};
	
	Mode mode = DualStereo;
};

void parseArgs(int argc, char ** argv, CalibrationConfigs & config) {
	cxxopts::Options options(argv[0], "App generate correspondences and triangles from label files.");
	options.add_options()
		("M,Mode", "[Optional]The calibration mode, should be: single, dual, dualStereo or onlyStereo. Defualt is dualStereo. For single camera only left camera will be used", cxxopts::value<std::string>())
		("S,Seq", "[Conditional Necessary]The .xml of input image sequences.", cxxopts::value<std::string>())
		("IL", "[Conditional Necessary]The path to folder containing the images of left camera.", cxxopts::value<std::string>())
		("RL", "[Conditional Necessary]The path to folder containing the images of right camera.", cxxopts::value<std::string>())
		("IP", "[Conditional Necessary]The file containing input intrinsic parameters of two camera. Only use this when you are in only Stereo mode.", cxxopts::value<std::string>())
		("O,Out", "[Necessary]The output file containing all the parameters you ased.", cxxopts::value<std::string>())
		("OutSeq", "[Optional]To output image sequences used to calibrate.", cxxopts::value<std::string>())
		("OutUndist", "[Optional]Wether to output undistorted images.Following the output path of undistorted.", cxxopts::value<std::string>(config.outUndistortedImgpath))
		("Sample", "[Optional]Wehter to sample the input set for calibration. Follows the number of samples.", cxxopts::value<int>(config.sampleNum))
		("ShowCorners", "[Optional]To show the corners.", cxxopts::value<bool>(config.showCorners))
		("bw", "[Optional]The number of columns of corners in calibration board", cxxopts::value<int>(config.boardWidth))
		("bh", "[Optional]The number of rows of corners in calibration board", cxxopts::value<int>(config.boardHeight))
		("ext", "[Optional]The extension name for the input images.", cxxopts::value<std::string>(config.fileExt))
		("cornersOutput", "[Optional]Output corners found in input images.", cxxopts::value<std::string>(config.cornersOutput))
		;
	try
	{
		auto result = options.parse(argc, argv);
		if (result.count("M"))
		{
			std::string modeStr = result["M"].as<std::string>();
			std::transform(modeStr.begin(), modeStr.end(), modeStr.begin(), ::tolower);
			if (modeStr == "single" ) {
				config.mode = CalibrationConfigs::Single;
			}
			else if (modeStr == "dual")
			{
				config.mode = CalibrationConfigs::Dual;
			}
			else if(modeStr == "dualstereo")
			{
				config.mode = CalibrationConfigs::DualStereo;
			}
			else
			{
				config.mode = CalibrationConfigs::OnlyStereo;
			}
		}
		if (config.mode == CalibrationConfigs::DualStereo || config.mode == CalibrationConfigs::Dual)
		{
			//For dual calibration
			if ((!result.count("S") && (!result.count("IL") || !result.count("RL")))
				|| !result.count("O"))
			{
				cout << "Input args not enough!\n";
				cout << options.help();
				exit(1);
			}
			if (result.count("S"))
			{
				config.inImgSequenceFile = result["S"].as<std::string>();
				config.inImgPathL = "";
				config.inImgPathR = "";
			}
			else
			{
				config.inImgPathL = result["IL"].as<std::string>();
				config.inImgPathR = result["RL"].as<std::string>();
			}
		}
		else if (config.mode == CalibrationConfigs::Single)
		{
			//For single calibration, only use left camera
			if ((!result.count("S") && (!result.count("IL")))
				|| !result.count("O"))
			{
				cout << "Input args not enough!\n";
				cout << options.help();
				exit(1);
			}
			if (result.count("S"))
			{
				config.inImgSequenceFile = result["S"].as<std::string>();
				config.inImgPathL = "";
				config.inImgPathR = "";
			}
			else
			{
				config.inImgPathL = result["IL"].as<std::string>();
			}
		}
		else if (config.mode == CalibrationConfigs::OnlyStereo)
		{
			if ((!result.count("S") && (!result.count("IL") || !result.count("RL")))
				|| !result.count("O") || !result.count("IP"))
			{
				cout << "Input args not enough!\n";
				cout << options.help();
				exit(1);
			}
			if (result.count("S"))
			{
				config.inImgSequenceFile = result["S"].as<std::string>();
				config.inImgPathL = "";
				config.inImgPathR = "";
			}
			else
			{
				config.inImgPathL = result["IL"].as<std::string>();
				config.inImgPathR = result["RL"].as<std::string>();
			}
			config.inIntrinsicsFile = result["IP"].as<std::string>();
		}
		
		config.outFile = result["O"].as<std::string>();
		if (result.count("Sample"))
		{
			config.resampleInputs = true;
		}

		if (result.count("OutUndist"))
		{
			config.undistordImages = true;
		}
		if (result.count("OutSeq"))
		{
			config.outputSelectedFileSequence = true;
			config.outImgSequenceFile = result["OutSeq"].as<std::string>();
		}
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		cout << options.help();
		exit(1);
	}

}

double calibrateSingleCam(AC::VecStr & files, cv::Size board_sz, cv::Mat & intrinsic, cv::Mat & distCoeffs, CalibrationConfigs & config) {
	int numSquares = board_sz.height * board_sz.width;
	int numCornersHor = board_sz.width;

	int numSuccess = 0;
	cv::Size imgSize;

	int successes = 0;
	VecP3f objP;
	for (int j = 0; j < numSquares; j++) {
		objP.push_back(CB_SQUARE_SIZE * cv::Point3f(j / numCornersHor, j%numCornersHor, 0.0f));
	}
	VecVecP2f imgPoints;
	VecVecP3f objPoints;
	
	int numThreads = config.showCorners ? 1 : NUM_THREADS;
#pragma omp parallel for num_threads(numThreads)
	for (int i = 0; i < files.size(); ++i)
	{
		std::string imgFile = files[i];
		VecP2f corners;
		

		cv::Mat img = cv::imread(imgFile);
		imgSize = img.size();
		cv::Mat imgGray;
		cv::cvtColor(img, imgGray, CV_BGR2GRAY);
		bool found = findChessboardCorners(img, board_sz, corners, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FILTER_QUADS | CV_CALIB_CB_FAST_CHECK);

#pragma omp critical
		{
			std::cout << "Processing: " << imgFile << std::endl;
		}
		if (found)
		{
			cornerSubPix(imgGray, corners, cv::Size(11, 11), cv::Size(-1, -1), cv::TermCriteria(CV_TERMCRIT_EPS | CV_TERMCRIT_ITER, 30, 0.001));
			if (config.showCorners)
			{
				drawChessboardCorners(img, board_sz, corners, found);
				cv::imshow("Corners", img);
				cv::waitKey(0);
			}
#pragma omp critical
			{
				imgPoints.push_back(corners);
				objPoints.push_back(objP);
				++numSuccess;
				std::cout << "Found in : " << imgFile << std::endl;
			}
		}
		else {
#pragma omp critical
			{
				std::cout << "Failed in : " << imgFile << std::endl;
			}
		}
	}
	std::vector<cv::Mat> rvecs;
	std::vector<cv::Mat> tvecs;

	std::cout << "Start calibration: " << std::endl;
	double rtval = calibrateCamera(objPoints, imgPoints, imgSize, intrinsic, distCoeffs, rvecs, tvecs);
	printf("Returning value: %f\n", rtval);
	std::cout << "Intrinsic matrix:\n" << intrinsic << std::endl;
	std::cout << "Dist Coeffs:\n" << distCoeffs << std::endl;

	return rtval;
};

double calibrateDualCam(AC::VecStr & imgsL, AC::VecStr & imgsR, AC::VecStr & imgsLSelected, AC::VecStr & imgsRSelected, cv::Size board_sz,
	cv::Mat & intrinsicL, cv::Mat & distCoeffsL, cv::Mat & intrinsicR, cv::Mat & distCoeffsR, bool stereoCalib, cv::Mat & R, cv::Mat & T, 
	CalibrationConfigs & config) {
	int numSquares = board_sz.height * board_sz.width;
	int numCornersHor = board_sz.width;

	int numSuccess = 0;
	cv::Size imgSize;

	int successes = 0;
	VecP3f objP;
	for (int j = 0; j < numSquares; j++) {
		objP.push_back(CB_SQUARE_SIZE * cv::Point3f(j / numCornersHor, j%numCornersHor, 0.0f));
	}

	VecVecP2f imgPointsL, imgPointsR;
	VecVecP3f objPointsL, objPointsR;

	assert(imgsL.size() == imgsR.size());

	int numThreads = config.showCorners ? 1 : NUM_THREADS;
#pragma omp parallel for num_threads(numThreads)
	for (int i  = 0; i < imgsL.size(); ++i)
	{
		VecP2f cornersL, cornersR;
		std::string imgFileL = imgsL[i];
		std::string imgFileR = imgsR[i];
		cv::Mat imgL = cv::imread(imgFileL);
		cv::Mat imgR = cv::imread(imgFileR);
		imgSize = imgL.size();
		cv::Mat imgGrayL, imgGrayR;
		cv::cvtColor(imgL, imgGrayL, CV_BGR2GRAY);
		cv::cvtColor(imgR, imgGrayR, CV_BGR2GRAY);
		bool found = findChessboardCorners(imgL, board_sz, cornersL, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FILTER_QUADS | CV_CALIB_CB_FAST_CHECK)
			&& findChessboardCorners(imgR, board_sz, cornersR, CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FILTER_QUADS | CV_CALIB_CB_FAST_CHECK);
#pragma omp critical
		{
			std::cout << "Processing: " << imgFileL << std::endl;
		}
		if (found)
		{
			cornerSubPix(imgGrayL, cornersL, cv::Size(11, 11), cv::Size(-1, -1), cv::TermCriteria(CV_TERMCRIT_EPS | CV_TERMCRIT_ITER, 30, 0.001));
			cornerSubPix(imgGrayR, cornersR, cv::Size(11, 11), cv::Size(-1, -1), cv::TermCriteria(CV_TERMCRIT_EPS | CV_TERMCRIT_ITER, 30, 0.001));
			drawChessboardCorners(imgL, board_sz, cornersL, found);
			drawChessboardCorners(imgR, board_sz, cornersR, found);

			if (config.cornersOutput != "")
			{
				AC::IO::FileParts fpL = AC::IO::fileparts(imgFileL), fpR = AC::IO::fileparts(imgFileR);
				cv::imwrite(config.cornersOutput + "/" + fpL.name + "L.png", imgL);
				cv::imwrite(config.cornersOutput + "/" + fpR.name + "R.png", imgR);
			}
			if (config.showCorners)
			{
				cv::Mat imgLShow, imgRShow;
				cv::resize(imgL, imgLShow, imgL.size() / 2);
				cv::resize(imgR, imgRShow, imgR.size() / 2);
				cv::imshow("CornersL", imgLShow);
				cv::imshow("CornersR", imgRShow);
				cv::waitKey(0);
			}
			

#pragma omp critical
			{
				imgPointsL.push_back(cornersL);
				imgPointsR.push_back(cornersR);
				objPointsL.push_back(objP);
				objPointsR.push_back(objP);
				++numSuccess;
				std::cout << "Found in : " << imgFileL << " and " << imgFileR << std::endl;
			}
		}
		else {
#pragma omp critical
			{
				std::cout << "Failed in : " << imgFileL << " or " << imgFileR << std::endl;
			}
		}
	}
	std::vector<cv::Mat> rvecsL, rvecsR;
	std::vector<cv::Mat> tvecsL, tvecsR;

	if (config.mode == CalibrationConfigs::Dual || config.mode == CalibrationConfigs::DualStereo)
	{
		std::cout << "Start calibration: " << std::endl;
		double rtvalL = calibrateCamera(objPointsL, imgPointsL, imgSize, intrinsicL, distCoeffsL, rvecsL, tvecsL);
		double rtvalR = calibrateCamera(objPointsR, imgPointsR, imgSize, intrinsicR, distCoeffsR, rvecsR, tvecsR);
		printf("Left Returning value: %f\n", rtvalL);
		printf("Right Returning value: %f\n", rtvalR);
		if (!stereoCalib)
		{
			return rtvalL + rtvalR;
		}
	}
	
	std::cout << "Left Intrinsic matrix:\n" << intrinsicL << std::endl;
	std::cout << "Right Intrinsic matrix:\n" << intrinsicR << std::endl;

	std::cout << "Left Dist Coeffs:\n" << distCoeffsL << std::endl;
	std::cout << "Right Dist Coeffs:\n" << distCoeffsR << std::endl;
	
	cv::Mat E, F;

	double rms = cv::stereoCalibrate(objPointsL, imgPointsL, imgPointsR,
		intrinsicL, distCoeffsL,
		intrinsicR, distCoeffsR,
		imgSize, R, T, E, F,
		cv::CALIB_FIX_ASPECT_RATIO +
		cv::CALIB_ZERO_TANGENT_DIST +
		cv::CALIB_FIX_INTRINSIC +
		//cv::CALIB_SAME_FOCAL_LENGTH +
		cv::CALIB_RATIONAL_MODEL +
		cv::CALIB_FIX_K1 + cv::CALIB_FIX_K2 + cv::CALIB_FIX_K3 + cv::CALIB_FIX_K4 + cv::CALIB_FIX_K5 + cv::CALIB_FIX_K6,
		cv::TermCriteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, 100, 1e-5));

	cout << "R:\n" << R << endl;
	cout << "T:\n" << T << endl;
	cout << "done with RMS error=" << rms << endl;

	return rms;
};

int myrandom(int i) { return std::rand() % i; }

void getRandomSampledIndices(int total, int numSample, AC::VecInt & indices) {
	indices.clear();
	indices.reserve(total);
	for (int i = 0; i < total; i++)
	{
		indices.push_back(i);
	}
	std::srand(unsigned(std::time(0)));
	//std::srand(12345);
	//std::srand(54321);
	std::random_shuffle(indices.begin(), indices.end(), myrandom);

	indices.resize(numSample);
}

void extractFilesFromIndices(AC::VecInt & indices, AC::VecStr & allFiles, AC::VecStr & extractedFiles) {
	extractedFiles.clear();
	for (int i = 0; i < indices.size(); i++)
	{
		extractedFiles.push_back(allFiles[indices[i]]);
	}
}


void undistortImgFiles(AC::VecStr & files, std::string undistortedPath, cv::Mat intrinsic, cv::Mat distCoeffs) {
	for (auto imgFile : files)
	{
		AC::IO::FileParts fp = AC::IO::fileparts(imgFile);
		cv::Mat img = cv::imread(imgFile);
		cv::Mat imgUndistorted;
		cv::undistort(img, imgUndistorted, intrinsic, distCoeffs);
	
		cout << "Writting to: " << imgFile << endl;
		cv::imwrite(undistortedPath + "\\" + fp.name + fp.ext, imgUndistorted);
		//cv::imshow("win1", img);
		//cv::imshow("win2", imgUndistorted);
		//cv::waitKey(0);
	}
}

static bool readimgList(const std::string& filename, std::vector<std::string>& imgsL, std::vector<std::string>& imgsR)
{
	imgsL.resize(0);
	imgsR.resize(0);

	cv::FileStorage fs(filename, cv::FileStorage::READ);
	if (!fs.isOpened())
		return false;
	cv::FileNode nL = fs["imgListLeft"];
	if (nL.type() != cv::FileNode::SEQ)
		return false;
	cv::FileNodeIterator itL = nL.begin(), itL_end = nL.end();
	for (; itL != itL_end; ++itL)
		imgsL.push_back((std::string)*itL);

	cv::FileNode nR = fs["imgListRight"];
	if (nR.type() != cv::FileNode::SEQ)
		return false;
	cv::FileNodeIterator itR = nR.begin(), itR_end = nR.end();
	for (; itR != itR_end; ++itR)
		imgsR.push_back((std::string)*itR);

	return true;
}

int main(int argc, char ** argv) {
	CalibrationConfigs config;
	parseArgs(argc, argv, config);
	AC::VecStr imgsL, imgsR;


	if (config.inImgSequenceFile != "")
	{
		readimgList(config.inImgSequenceFile, imgsL, imgsR);
	}
	else
	{
		AC::IO::getFilesWithExt(config.inImgPathL, config.fileExt, imgsL);
		std::sort(imgsL.begin(), imgsL.end());

		if (config.mode != CalibrationConfigs::Single)
		{
			AC::IO::getFilesWithExt(config.inImgPathR, config.fileExt, imgsR);
			std::sort(imgsR.begin(), imgsR.end());
		}

	}

	cv::Size boardSize(config.boardWidth, config.boardHeight);

	AC::VecInt indices;
	AC::VecStr sampledImgsL, sampledImgsR;
	AC::VecStr imgsLSelected, imgsRSelected;

	if (config.resampleInputs)
	{
		switch (config.resampleMode)
		{
		case CalibrationConfigs::random:
			getRandomSampledIndices(imgsL.size(), config.sampleNum, indices);
			extractFilesFromIndices(indices, imgsL, sampledImgsL);
			if (config.mode != CalibrationConfigs::Single)
			{
				extractFilesFromIndices(indices, imgsR, sampledImgsR);
			}
			break;
		case CalibrationConfigs::skip:
			break;
		default:
			break;
		}
	}
	else
	{
		sampledImgsL = imgsL;
		sampledImgsR = imgsR;
	}
	cv::Mat intrinsicL, intrinsicR, distCoeffsL, distCoeffsR, R, t;
	if (config.mode == CalibrationConfigs::Dual || config.mode == CalibrationConfigs::DualStereo)
	{
		bool stereoCalib = config.mode == CalibrationConfigs::DualStereo;
		calibrateDualCam(sampledImgsL, sampledImgsR, imgsLSelected, imgsRSelected, boardSize,
			intrinsicL, distCoeffsL, intrinsicR, distCoeffsR, stereoCalib, R, t, config);

		cv::FileStorage fs(config.outFile, cv::FileStorage::WRITE);
		if (fs.isOpened())
		{
			fs << "M1" << intrinsicL << "D1" << distCoeffsL <<
				"M2" << intrinsicR << "D2" << distCoeffsR << "R" << R << "T" << t;
			fs.release();
		}
	}
	else if (config.mode == CalibrationConfigs::Single)
	{
		calibrateSingleCam(sampledImgsL, boardSize, intrinsicL, distCoeffsL, config);
		cv::FileStorage fs(config.outFile, cv::FileStorage::WRITE);
		if (fs.isOpened())
		{
			fs << "M1" << intrinsicL << "D1" << distCoeffsL;
			fs.release();
		}
		if (config.undistordImages)
		{
			mkdir(config.outUndistortedImgpath.c_str());
			undistortImgFiles(imgsL, config.outUndistortedImgpath, intrinsicL, distCoeffsL);
		}
	}
	else if (config.mode == CalibrationConfigs::OnlyStereo)
	{
		cv::FileStorage fsIn(config.inIntrinsicsFile, cv::FileStorage::READ);
		fsIn["M1"] >> intrinsicL;
		fsIn["D1"] >> distCoeffsL;
		fsIn["M2"] >> intrinsicR;
		fsIn["D2"] >> distCoeffsR;

		bool stereoCalib = config.stereoCalib == CalibrationConfigs::DualStereo;
		calibrateDualCam(sampledImgsL, sampledImgsR, imgsLSelected, imgsRSelected, boardSize,
			intrinsicL, distCoeffsL, intrinsicR, distCoeffsR, stereoCalib, R, t, config);

		cv::FileStorage fs(config.outFile, cv::FileStorage::WRITE);
		if (fs.isOpened())
		{
			fs << "M1" << intrinsicL << "D1" << distCoeffsL <<
				"M2" << intrinsicR << "D2" << distCoeffsR << "R" << R << "T" << t;
			fs.release();
		}
	}
	
	//Generate and save the undistored images
	//std::string undistortedPath = std::string(argv[1]) + "/Undistorted/";
	//mkdir(undistortedPath.c_str());
	//cv::Mat imgUndistorted;
	//getchar();
}