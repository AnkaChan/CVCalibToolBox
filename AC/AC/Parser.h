#ifndef _AC_PARSER_H_
#define _AC_PARSER_H_
#include "cxxopts.hpp"

namespace AC {
	
	struct OutputSelection
	{
		bool toOutput = false;
		std::string path;
			
	};

};

//void parseArgs(int argc, char ** argv, WhitePatternExtractionConfig & config) {
//	cxxopts::Options options(argv[0], "App generate correspondences and triangles from label files.");
//	options.add_options()
//		("m, maxReproj", "[Optional]The max cut off reprojection error. Defualt is 6.0.", cxxopts::value<float>(config.maxReprojErr))
//		("r, reproj", "[Optional]Wether to output the reprojection error. Must be followed by the output Reproj file path.", cxxopts::value<std::string>());
//
//	try
//	{
//		auto result = options.parse(argc, argv);
//
//		if (result.count("r")) {
//			config.outputReprojError = true;
//			config.reprojErrorFile = result["r"].as<std::string>();
//		}
//
//	}
//	catch (const cxxopts::OptionException& e)
//	{
//		std::cout << "error parsing options: " << e.what() << std::endl;
//		cout << options.help();
//		exit(1);
//	}
//
//}

#endif //_PATTERN_OCR_H_