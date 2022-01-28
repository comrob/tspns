/*
 * File name: tgsoa.cc
 * Date:      2016/12/07 08:33
 * Author:    Jan Faigl
 */

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <crl/logging.h>
#include <crl/perf_timer.h>

#include <crl/gui/guifactory.h>
#include <crl/gui/win_adjust_size.h>

#include <crl/boost_args_config.h>
#include <crl/config.h>

#include "gsoa.h"
#include "glkh.h"

using crl::logger;
using namespace gsoa;
using namespace glkh;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

const std::string GSOA_VERSION = "1.0";

typedef crl::gui::CCanvasBase Canvas;

/// ----------------------------------------------------------------------------
/// Program options variables
/// ----------------------------------------------------------------------------
std::string guiType = "none";

crl::CConfig glkhConfig;
crl::CConfig gsoaConfig;
std::string canvasOutput = "";
std::string solutionFile = "";

/// ----------------------------------------------------------------------------
/// Global variable
/// ----------------------------------------------------------------------------
crl::gui::CGui *g = 0;
#define GUI(x) \
	if (gui)   \
	{          \
		x;     \
	}

/// ----------------------------------------------------------------------------
bool parseArgs(int argc, char *argv[])
{
	bool ret = true;
	std::string configFile;
	std::string glkhConfigFile;
	std::string guiConfigFile;
	std::string loggerCfg = "";

	po::options_description desc("General options");
	desc.add_options()("help,h", "produce help message")
	("config,c", po::value<std::string>(&configFile)->default_value(std::string(argv[0]) + ".cfg"), "configuration file")
	// ("config,c", po::value<std::string>(&glkhConfigFile)->default_value(std::string(argv[0]) + "-glkh.cfg"), "configuration file")
	("logger-config,l", po::value<std::string>(&loggerCfg)->default_value(loggerCfg), "logger configuration file")("solution-file", po::value<std::string>(&solutionFile)->default_value(""));
	try
	{
		// boost_args_add_options(guiConfig, "", guiOptions);
		// guiOptions.add_options()("canvas-output", po::value<std::string>(&canvasOutput), "result canvas outputfile");

		po::options_description gsoaOptions("GSOA options");
		boost_args_add_options(CGSOA::getConfig(gsoaConfig), "", gsoaOptions);

		po::options_description cmdline_options;
		cmdline_options.add(desc).add(gsoaOptions);

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
		po::notify(vm);

		std::ifstream ifs(configFile.c_str());
		store(parse_config_file(ifs, cmdline_options), vm);
		po::notify(vm);
		ifs.close();

		if (vm.count("help"))
		{
			std::cerr << std::endl;
			std::cerr << "TSPNS GSOA solver ver. " << GSOA_VERSION << std::endl;
			std::cerr << cmdline_options << std::endl;
			ret = false;
		}
		if (
			ret &&
			loggerCfg != "" &&
			fs::exists(fs::path(loggerCfg)))
		{
			crl::initLogger("gsoa", loggerCfg.c_str());
		}
		else
		{
			crl::initLogger("gsoa");
		}
		const std::string problemFile = gsoaConfig.get<std::string>("problem");
		if (!fs::exists(fs::path(problemFile)))
		{
			ERROR("Problem file '" + problemFile + "' does not exists");
			ret = false;
		}
	}
	catch (std::exception &e)
	{
		std::cerr << std::endl;
		std::cerr << "Error in parsing arguments: " << e.what() << std::endl;
		ret = false;
	}
	return ret;
}

/// - main ---------------------------------------------------------------------
int main(int argc, char *argv[])
{
	Canvas *canvas = 0;
	int ret = -1;
	if (parseArgs(argc, argv))
	{
		INFO("Start Logging");
		try
		{
			crl::CPerfTimer t("Load problem time real:");

			if (gsoaConfig.get<bool>("lkh") || gsoaConfig.get<bool>("glkh"))
			{
				CGLKH glkh(gsoaConfig);
				{
					t.start("TSP solve time: ");
					glkh.solve();
					t.stop();
				}
			}
			else
			{
				CGSOA gsoa(gsoaConfig);
				{
					t.start("TSP solve time: ");
					gsoa.solve();
					t.stop();
				}
			}

			INFO("End Logging");
		}
		catch (crl::exception &e)
		{
			ERROR("Exception " << e.what() << "!");
		}
		catch (std::exception &e)
		{
			ERROR("Runtime error " << e.what() << "!");
		}
		ret = EXIT_SUCCESS;
	}
	crl::shutdownLogger();
	return ret;
}

/* end of tgsoa.cc */
